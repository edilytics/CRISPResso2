"""Unit tests for CRISPResso2.storage — VariantStore Stage 1 (count_reads).

Parity strategy: each test builds the *reference* dedup dict using the exact
arithmetic the current pandas path uses in ``CRISPRessoCORE.process_fastq`` and
``process_paired_fastq`` (single-read: ``{seq: count}``; paired:
``{R1+'+'+RC(R2): [count, quals]}``), then asserts that
``VariantStore.count_reads`` reproduces it byte-for-byte on the eager path and
that the spill path reproduces the same counts (and first-occurrence quals).
"""
import gzip
import os

import numpy as np
import pytest

from CRISPResso2 import CRISPRessoShared
from CRISPResso2.storage import (
    AlignedShardWriter,
    DEFAULT_MEMORY_BUDGET_MB,
    ReadCounts,
    VariantStore,
    count_reads_from_fastq,
    iter_aligned_shard,
    payload_to_row,
    variant_parquet_generator_process,
)

_RC = CRISPRessoShared.reverse_complement


def _write_fastq(path, records):
    """records: list of (name, seq, qual)."""
    with open(path, "w") as f:
        for name, seq, qual in records:
            f.write(f"@{name}\n{seq}\n+\n{qual}\n")
    return path


def _reference_single(records):
    """Mirror process_fastq: {seq: count}."""
    d = {}
    for _, seq, _ in records:
        d[seq] = d.get(seq, 0) + 1
    return d


def _reference_paired(records1, records2):
    """Mirror process_paired_fastq: {R1+'+'+RC(R2): [count, quals]}.

    quals = R1_qual + ' ' + reversed(R2_qual); first-occurrence retained.
    """
    d = {}
    for (_, s1, q1), (_, s2, q2) in zip(records1, records2):
        key = s1 + "+" + _RC(s2)
        quals = q1 + " " + q2[::-1]
        if key in d:
            d[key][0] += 1
        else:
            d[key] = [1, quals]
    return d


def _gz(path, text):
    with gzip.open(path, "wt") as f:
        f.write(text)
    return path


# -- single-read eager parity -------------------------------------------------


def test_single_read_eager_parity(temp_dir):
    recs = [
        ("r1", "ATCGATCG", "IIIIIIII"),
        ("r2", "ATCGATCG", "IIIIIIII"),  # dup of r1
        ("r3", "GCTAGCTA", "IIIIIIII"),
        ("r4", "TTTTAAAA", "IIIIIIII"),
    ]
    fq = _write_fastq(os.path.join(temp_dir, "in.fastq"), recs)
    rc = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False)
    assert not rc.spilled
    assert rc.num_total == 4
    assert rc.num_unique == 3
    assert rc.to_dict() == _reference_single(recs)


def test_single_read_items_yields_count_and_none_quals(temp_dir):
    recs = [("r1", "ATCG", "IIII"), ("r2", "ATCG", "IIII"), ("r3", "GGGG", "IIII")]
    fq = _write_fastq(os.path.join(temp_dir, "in.fastq"), recs)
    rc = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False)
    items = sorted(rc.items())
    assert items == [("ATCG", 2, None), ("GGGG", 1, None)]


def test_single_read_eager_with_quals_parity(temp_dir):
    recs = [
        ("r1", "ATCGATCG", "IIIIIIII"),
        ("r2", "ATCGATCG", "HHHHHHHH"),  # dup seq, different quals -> first kept
        ("r3", "GCTAGCTA", "IIIIIIII"),
    ]
    fq = _write_fastq(os.path.join(temp_dir, "in.fastq"), recs)
    rc = count_reads_from_fastq(fq, None, temp_dir, keep_quals=True)
    assert not rc.spilled
    expected = _reference_single(recs)
    # to_dict_with_quals: {seq: [count, first_quals]}
    expected_q = {k: [v, recs[[r[1] for r in recs].index(k)][2]] for k, v in expected.items()}
    assert rc.to_dict_with_quals() == expected_q
    # count column correct
    assert rc.to_dict() == expected


# -- paired eager parity ------------------------------------------------------


def test_paired_eager_parity(temp_dir):
    r1 = [("a", "AAATTT", "IIIIII"), ("b", "AAATTT", "IIIIII"), ("c", "CCCGGG", "IIIIII")]
    r2 = [("a", "TTTAAA", "IIIIII"), ("b", "TTTAAA", "IIIIII"), ("c", "GGGCCC", "IIIIII")]
    fq1 = _write_fastq(os.path.join(temp_dir, "r1.fastq"), r1)
    fq2 = _write_fastq(os.path.join(temp_dir, "r2.fastq"), r2)
    rc = count_reads_from_fastq(fq1, fq2, temp_dir, keep_quals=True)
    assert not rc.spilled
    assert rc.num_total == 3
    assert rc.num_unique == 2
    assert rc.to_dict_with_quals() == _reference_paired(r1, r2)


def test_paired_key_uses_reverse_complement_of_r2(temp_dir):
    # R2 raw "AAATTT"; RC = "AAATTT" (palindrome) -> key = R1 + "+" + "AAATTT"
    r1 = [("a", "GGGCCC", "IIIIII")]
    r2 = [("a", "AAATTT", "IIIIII")]
    fq1 = _write_fastq(os.path.join(temp_dir, "r1.fastq"), r1)
    fq2 = _write_fastq(os.path.join(temp_dir, "r2.fastq"), r2)
    rc = count_reads_from_fastq(fq1, fq2, temp_dir, keep_quals=True)
    keys = list(rc.to_dict_with_quals().keys())
    assert keys == ["GGGCCC+AAATTT"]
    # quals = R1_qual + ' ' + reversed(R2_qual)
    assert rc.to_dict_with_quals()["GGGCCC+AAATTT"][1] == "IIIIII IIIIII"


def test_paired_quals_first_occurrence_retained(temp_dir):
    r1 = [("a", "GGGCCC", "IIIIII"), ("b", "GGGCCC", "HHHHHH")]
    r2 = [("a", "AAATTT", "JJJJJJ"), ("b", "AAATTT", "KKKKKK")]
    fq1 = _write_fastq(os.path.join(temp_dir, "r1.fastq"), r1)
    fq2 = _write_fastq(os.path.join(temp_dir, "r2.fastq"), r2)
    rc = count_reads_from_fastq(fq1, fq2, temp_dir, keep_quals=True)
    d = rc.to_dict_with_quals()
    key = "GGGCCC+AAATTT"
    assert d[key][0] == 2  # both reads same key
    assert d[key][1] == "IIIIII JJJJJJ"  # first occurrence quals


# -- spill path ---------------------------------------------------------------


def test_spill_forced_reproduces_eager_counts_single(temp_dir):
    recs = [
        ("r%d" % i, "ATCG" if i % 2 == 0 else "GCTA", "IIII") for i in range(200)
    ]
    fq = _write_fastq(os.path.join(temp_dir, "in.fastq"), recs)
    eager = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False)
    spill = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False, force_spill=True)
    assert spill.spilled
    assert spill.parquet_path and os.path.exists(spill.parquet_path)
    assert spill.num_total == eager.num_total
    assert spill.num_unique == eager.num_unique
    assert dict((k, c) for k, c, _ in spill.items()) == eager.to_dict()


def test_spill_forced_reproduces_eager_with_quals(temp_dir):
    recs = [("r%d" % i, "ATCG" if i % 2 == 0 else "GCTA", "I%dII" % (i % 4)) for i in range(200)]
    fq = _write_fastq(os.path.join(temp_dir, "in.fastq"), recs)
    eager = count_reads_from_fastq(fq, None, temp_dir, keep_quals=True)
    spill = count_reads_from_fastq(fq, None, temp_dir, keep_quals=True, force_spill=True)
    assert spill.spilled
    # counts match
    assert dict((k, c) for k, c, _ in spill.items()) == eager.to_dict()
    # first-occurrence quals match (stable sort preserves file order within group)
    eager_q = eager.to_dict_with_quals()
    spill_q = {k: [c, q] for k, c, q in spill.items()}
    assert spill_q == eager_q


def test_spill_paired_reproduces_eager(temp_dir):
    r1 = [("r%d" % i, "AAATTT" if i % 3 == 0 else "CCCGGG", "I%dII" % (i % 5)) for i in range(150)]
    r2 = [("r%d" % i, "TTTAAA" if i % 3 == 0 else "GGGCCC", "I%dII" % (i % 5)) for i in range(150)]
    fq1 = _write_fastq(os.path.join(temp_dir, "r1.fastq"), r1)
    fq2 = _write_fastq(os.path.join(temp_dir, "r2.fastq"), r2)
    eager = count_reads_from_fastq(fq1, fq2, temp_dir, keep_quals=True)
    spill = count_reads_from_fastq(fq1, fq2, temp_dir, keep_quals=True, force_spill=True)
    assert spill.num_unique == eager.num_unique
    assert {k: [c, q] for k, c, q in spill.items()} == eager.to_dict_with_quals()


def test_eager_threshold_triggers_spill_on_large_input(temp_dir):
    # 5kb reads with ~all-unique keys; tiny memory budget forces spill.
    import random
    rng = random.Random(0)
    recs = []
    for i in range(2000):
        seq = "".join(rng.choice("ACGT") for _ in range(5000))
        recs.append(("r%d" % i, seq, "I" * 5000))
    fq = _write_fastq(os.path.join(temp_dir, "big.fastq"), recs)
    rc = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False, memory_budget_mb=1)
    assert rc.spilled
    assert rc.num_total == 2000
    assert rc.num_unique == 2000
    # all counts == 1
    assert all(c == 1 for _, c, _ in rc.items())


def test_eager_stays_in_memory_under_budget(temp_dir):
    recs = [("r%d" % i, "ATCG%d" % (i % 4), "IIIII") for i in range(100)]
    fq = _write_fastq(os.path.join(temp_dir, "small.fastq"), recs)
    rc = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False, memory_budget_mb=DEFAULT_MEMORY_BUDGET_MB)
    assert not rc.spilled
    assert rc._counts is not None


def test_budget_triggered_spill_with_quals_parity(temp_dir):
    # Non-forced spill: tiny budget forces the eager dict over the threshold.
    import random
    rng = random.Random(42)
    recs = [("r%d" % i, "".join(rng.choice("ACGT") for _ in range(2000)), "I" * 2000) for i in range(500)]
    fq = _write_fastq(os.path.join(temp_dir, "in.fastq"), recs)
    eager = count_reads_from_fastq(fq, None, temp_dir, keep_quals=True, memory_budget_mb=DEFAULT_MEMORY_BUDGET_MB)
    spill = count_reads_from_fastq(fq, None, temp_dir, keep_quals=True, memory_budget_mb=1)
    assert spill.spilled, "budget should have forced spill"
    # First-occurrence quals must match the eager (file-order) dict.
    assert {k: [c, q] for k, c, q in spill.items()} == eager.to_dict_with_quals()


def test_spill_count_vector_matches_eager_fuzz(temp_dir):
    # Fuzz: many duplicates + unique reads; spill and eager must agree on every count.
    import random
    rng = random.Random(7)
    recs = []
    for i in range(1000):
        # ~70% draw from a small alphabet (heavy dups), ~30% unique.
        if rng.random() < 0.7:
            seq = "SEQ%d" % rng.randint(0, 19)
        else:
            seq = "UNIQ%d" % i
        recs.append(("r%d" % i, seq, "I" * len(seq)))
    fq = _write_fastq(os.path.join(temp_dir, "fuzz.fastq"), recs)
    eager = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False)
    spill = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False, force_spill=True)
    eager_d = eager.to_dict()
    spill_d = dict((k, c) for k, c, _ in spill.items())
    assert spill_d == eager_d
    assert spill.num_total == eager.num_total == len(recs)
    assert spill.num_unique == eager.num_unique == len(eager_d)


# -- edge cases ---------------------------------------------------------------


def test_empty_fastq(temp_dir):
    fq = _write_fastq(os.path.join(temp_dir, "empty.fastq"), [])
    rc = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False)
    assert rc.num_total == 0
    assert rc.num_unique == 0
    assert list(rc.items()) == []


def test_empty_fastq_spill(temp_dir):
    fq = _write_fastq(os.path.join(temp_dir, "empty.fastq"), [])
    rc = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False, force_spill=True)
    assert rc.spilled
    assert rc.num_total == 0
    assert rc.num_unique == 0
    assert list(rc.items()) == []


def test_single_read(temp_dir):
    fq = _write_fastq(os.path.join(temp_dir, "one.fastq"), [("r1", "ATCG", "IIII")])
    rc = count_reads_from_fastq(fq, None, temp_dir, keep_quals=True)
    assert rc.num_total == 1
    assert rc.to_dict_with_quals() == {"ATCG": [1, "IIII"]}


def test_gzip_fastq(temp_dir):
    recs = [("r%d" % i, "ATCG" if i % 2 == 0 else "GCTA", "IIII") for i in range(20)]
    text = "".join(f"@{n}\n{s}\n+\n{q}\n" for n, s, q in recs)
    fq = _gz(os.path.join(temp_dir, "in.fastq.gz"), text)
    eager = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False)
    spill = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False, force_spill=True)
    assert eager.to_dict() == spill.to_dict() if not spill.spilled else (
        eager.to_dict() == dict((k, c) for k, c, _ in spill.items())
    )
    assert eager.num_unique == 2


def test_keep_quals_false_drops_quals(temp_dir):
    recs = [("r1", "ATCG", "IIII"), ("r2", "ATCG", "HHHH")]
    fq = _write_fastq(os.path.join(temp_dir, "in.fastq"), recs)
    rc = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False)
    with pytest.raises(RuntimeError):
        rc.to_dict_with_quals()
    for _, _, q in rc.items():
        assert q is None


def test_to_dict_raises_when_spilled(temp_dir):
    fq = _write_fastq(os.path.join(temp_dir, "in.fastq"), [("r1", "ATCG", "IIII")])
    rc = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False, force_spill=True)
    assert rc.spilled
    with pytest.raises(RuntimeError):
        rc.to_dict()
    with pytest.raises(RuntimeError):
        rc.to_dict_with_quals()


def test_paired_mismatched_lengths_raises(temp_dir):
    fq1 = _write_fastq(os.path.join(temp_dir, "r1.fastq"), [("a", "ATCG", "IIII")])
    fq2 = _write_fastq(os.path.join(temp_dir, "r2.fastq"), [
        ("a", "ATCG", "IIII"), ("b", "GGGG", "IIII")
    ])
    with pytest.raises(CRISPRessoShared.InputFileFormatException):
        count_reads_from_fastq(fq1, fq2, temp_dir, keep_quals=True)


def test_store_constructor_creates_directory(tmp_path):
    out = tmp_path / "nested" / "store"
    store = VariantStore(str(out), keep_quals=False)
    assert os.path.isdir(str(out))
    assert store.keep_quals is False
    assert store.memory_budget_bytes == DEFAULT_MEMORY_BUDGET_MB * 1024 * 1024


def test_read_counts_dataclass_fields():
    rc = ReadCounts(num_total=5, num_unique=3, spilled=False, keep_quals=True, _counts={"a": [2, "q"]})
    assert rc.num_total == 5
    assert rc.num_unique == 3
    assert rc.spilled is False


def test_items_eager_keep_quals_false_yields_none(temp_dir):
    fq = _write_fastq(os.path.join(temp_dir, "in.fastq"), [("r1", "ATCG", "IIII")])
    rc = count_reads_from_fastq(fq, None, temp_dir, keep_quals=False)
    assert list(rc.items()) == [("ATCG", 1, None)]


# -- Stage 2: schema round-trip (Spike S2) -----------------------------------

def _make_realistic_payload(ref_name="ref1", aln_len=20):
    """Build a payload with every field type from find_indels_substitutions.

    Includes numpy arrays, Python lists, tuples (coordinates), and all the
    scalar fields added by get_new_variant_object. This is the S2 spike —
    if every field round-trips through parquet, the schema is correct.
    """
    sub_payload = {
        "ref_name": ref_name,
        "aln_seq": "ATCGATCG--ATCGATCGAT",
        "aln_ref": "ATCGATCGAAATCGATCGAT",
        "aln_strand": "+",
        "classification": "MODIFIED",
        "irregular_ends": False,
        "insertions_outside_window": 1,
        "deletions_outside_window": 0,
        "substitutions_outside_window": 1,
        "total_mods": 3,
        "mods_in_window": 2,
        "mods_outside_window": 1,
        "insertion_n": 0,
        "deletion_n": 2,
        "substitution_n": 1,
        "ref_positions": [0, 1, 2, 3, 4, 5, 6, 7, -8, -8, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
        "all_insertion_positions": [7, 8],
        "all_insertion_left_positions": [7],
        "insertion_positions": [],
        "insertion_coordinates": [],
        "insertion_sizes": [],
        "all_deletion_positions": [8, 9],
        "deletion_positions": [8, 9],
        "deletion_coordinates": [(8, 10)],
        "deletion_sizes": [2],
        "all_substitution_positions": [4],
        "substitution_positions": [4],
        "substitution_values": np.array(["G"]),
    }
    return {
        "count": 1,
        "aln_ref_names": [ref_name],
        "aln_scores": [95.123, 12.0],
        "best_match_score": 95.123,
        "class_name": ref_name + "_MODIFIED",
        "best_match_name": ref_name,
        "caching_is_ok": True,
        "variant_" + ref_name: sub_payload,
    }


def test_s2_round_trip_all_fields(temp_dir):
    """Spike S2: every payload field round-trips through parquet with correct type."""
    payload = _make_realistic_payload()
    row = payload_to_row("ATCG+CGAT", 5, payload)
    shard = os.path.join(temp_dir, "test.parquet")
    with AlignedShardWriter(shard) as w:
        w.write_row(row)

    # read back
    items = list(iter_aligned_shard(shard))
    assert len(items) == 1
    read_key, count, back = items[0]
    assert read_key == "ATCG+CGAT"
    assert count == 5
    assert back["best_match_score"] == 95.123
    assert back["class_name"] == "ref1_MODIFIED"
    assert back["aln_ref_names"] == ["ref1"]
    assert back["aln_scores"] == [95.123, 12.0]
    assert back["caching_is_ok"] is True

    sub = back["variant_ref1"]
    assert sub["aln_seq"] == "ATCGATCG--ATCGATCGAT"
    assert sub["aln_ref"] == "ATCGATCGAAATCGATCGAT"
    assert sub["classification"] == "MODIFIED"
    assert sub["irregular_ends"] is False
    assert sub["insertion_n"] == 0
    assert sub["deletion_n"] == 2
    assert sub["substitution_n"] == 1
    assert sub["ref_positions"] == [0, 1, 2, 3, 4, 5, 6, 7, -8, -8, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
    assert sub["all_insertion_positions"] == [7, 8]
    assert sub["all_insertion_left_positions"] == [7]
    assert sub["insertion_coordinates"] == []
    assert sub["deletion_coordinates"] == [(8, 10)]
    assert sub["deletion_sizes"] == [2]
    assert sub["substitution_positions"] == [4]
    # substitution_values comes back as numpy array
    assert isinstance(sub["substitution_values"], np.ndarray)
    assert sub["substitution_values"].tolist() == ["G"]


def test_s2_round_trip_multi_ref(temp_dir):
    """A read aligning to 2 refs produces a list-of-structs with 2 entries."""
    p1 = _make_realistic_payload("ref1")
    p2 = _make_realistic_payload("ref2")
    payload = {
        "count": 1,
        "aln_ref_names": ["ref1", "ref2"],
        "aln_scores": [95.0, 95.0],
        "best_match_score": 95.0,
        "class_name": "ref1_MODIFIED&ref2_MODIFIED",
        "best_match_name": "ref1",
        "caching_is_ok": True,
        "variant_ref1": p1["variant_ref1"],
        "variant_ref2": p2["variant_ref2"],
    }
    row = payload_to_row("ATCG", 3, payload)
    shard = os.path.join(temp_dir, "multi.parquet")
    with AlignedShardWriter(shard) as w:
        w.write_row(row)
    items = list(iter_aligned_shard(shard))
    assert len(items) == 1
    _, _, back = items[0]
    assert back["aln_ref_names"] == ["ref1", "ref2"]
    assert "variant_ref1" in back
    assert "variant_ref2" in back
    assert back["variant_ref1"]["ref_name"] == "ref1"
    assert back["variant_ref2"]["ref_name"] == "ref2"


def test_s2_round_trip_unaligned_read(temp_dir):
    """Unaligned reads (best_match_score <= 0) have null payloads."""
    payload = {
        "count": 1,
        "aln_scores": [5.0, 3.0],
        "best_match_score": 0,
    }
    row = payload_to_row("NNNN", 2, payload)
    shard = os.path.join(temp_dir, "unaln.parquet")
    with AlignedShardWriter(shard) as w:
        w.write_row(row)
    items = list(iter_aligned_shard(shard))
    assert len(items) == 1
    _, count, back = items[0]
    assert count == 2
    assert back["best_match_score"] == 0.0
    assert back["aln_ref_names"] == []
    assert "variant_ref1" not in back
    assert back["caching_is_ok"] is True  # default


def test_s2_round_trip_empty_arrays(temp_dir):
    """Reads with no indels/substitutions have empty arrays, not null."""
    payload = _make_realistic_payload()
    # zero out all indel/substitution arrays
    sub = payload["variant_ref1"]
    sub["insertion_coordinates"] = []
    sub["deletion_coordinates"] = []
    sub["deletion_sizes"] = []
    sub["insertion_sizes"] = []
    sub["substitution_values"] = np.array([])
    sub["deletion_n"] = 0
    sub["insertion_n"] = 0
    sub["substitution_n"] = 0
    sub["all_deletion_positions"] = []
    sub["deletion_positions"] = []
    row = payload_to_row("ATCG", 1, payload)
    shard = os.path.join(temp_dir, "empty_arrays.parquet")
    with AlignedShardWriter(shard) as w:
        w.write_row(row)
    _, _, back = next(iter(iter_aligned_shard(shard)))
    sub_back = back["variant_ref1"]
    assert sub_back["deletion_coordinates"] == []
    assert sub_back["insertion_coordinates"] == []
    assert sub_back["deletion_sizes"] == []
    assert isinstance(sub_back["substitution_values"], np.ndarray)
    assert sub_back["substitution_values"].tolist() == []


def test_s2_round_trip_numpy_int_arrays(temp_dir):
    """Numpy int arrays (from find_indels_substitutions) round-trip correctly."""
    payload = _make_realistic_payload()
    sub = payload["variant_ref1"]
    sub["ref_positions"] = np.array([0, 1, 2, 3], dtype=np.int32)
    sub["all_insertion_positions"] = np.array([2, 3], dtype=np.int64)
    sub["all_substitution_positions"] = np.array([1], dtype=np.int8)
    row = payload_to_row("ATCG", 1, payload)
    shard = os.path.join(temp_dir, "numpy.parquet")
    with AlignedShardWriter(shard) as w:
        w.write_row(row)
    _, _, back = next(iter(iter_aligned_shard(shard)))
    sub_back = back["variant_ref1"]
    assert sub_back["ref_positions"] == [0, 1, 2, 3]
    assert sub_back["all_insertion_positions"] == [2, 3]
    assert sub_back["all_substitution_positions"] == [1]


def test_empty_shard_tolerated(temp_dir):
    """A worker whose slice has 0 aligned reads writes a valid empty parquet (edge #6)."""
    shard = os.path.join(temp_dir, "empty.parquet")
    with AlignedShardWriter(shard) as w:
        pass  # no rows written
    assert os.path.exists(shard)
    items = list(iter_aligned_shard(shard))
    assert items == []


def test_shard_batched_write(temp_dir):
    """Multiple batches flush correctly (batch_size < num_rows)."""
    shard = os.path.join(temp_dir, "batched.parquet")
    payloads = []
    with AlignedShardWriter(shard, batch_size=5) as w:
        for i in range(12):
            p = _make_realistic_payload()
            p["best_match_score"] = float(i)
            row = payload_to_row(f"seq_{i}", i + 1, p)
            w.write_row(row)
            payloads.append((f"seq_{i}", i + 1, i))
    items = list(iter_aligned_shard(shard))
    assert len(items) == 12
    for i, (rk, count, back) in enumerate(items):
        assert rk == f"seq_{i}"
        assert count == i + 1
        assert back["best_match_score"] == float(i)


def _mock_aligner(args, seq, refs, ref_names, aln_matrix, pe_info):
    """Module-level mock aligner (picklable for multiprocessing tests)."""
    return dict(_make_realistic_payload())


class _MockArgs:
    crispresso_merge = False


def _mock_aligner_paired(args, s1, s2, q1, q2, refs, ref_names, aln_matrix, pe_info):
    """Module-level mock aligner for paired reads."""
    return dict(_make_realistic_payload())


def test_variant_parquet_generator_process_single_read(temp_dir):
    """Worker writes a parquet shard using a mock aligner (single-read path)."""
    read_items = [("ATCG", 3, None), ("GCTA", 1, None)]
    variant_parquet_generator_process(
        read_items, _mock_aligner, _MockArgs(), {}, ["ref1"], None, (0, None),
        0, temp_dir, is_paired=False,
    )
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    assert os.path.exists(shard)
    items = list(iter_aligned_shard(shard))
    assert len(items) == 2
    assert items[0][0] == "ATCG"
    assert items[0][1] == 3  # count carried from read_items
    assert items[1][0] == "GCTA"
    assert items[1][1] == 1


def test_variant_parquet_generator_process_paired(temp_dir):
    """Worker splits read_key on '+' for paired reads (matching the TSV worker)."""
    seen = []

    def mock_aligner(args, s1, s2, q1, q2, refs, ref_names, aln_matrix, pe_info):
        seen.append((s1, s2, q1, q2))
        return dict(_make_realistic_payload())

    class MockArgsPaired:
        crispresso_merge = True

    read_items = [("AAATTT+GGGCCC", 2, "IIIIII JJJJJJ")]
    variant_parquet_generator_process(
        read_items, mock_aligner, MockArgsPaired(), {}, ["ref1"], None, (0, None),
        0, temp_dir, is_paired=True,
    )
    assert seen == [("AAATTT", "GGGCCC", "IIIIII", "JJJJJJ")]
    items = list(iter_aligned_shard(os.path.join(temp_dir, "aligned_0.parquet")))
    assert items[0][1] == 2  # count carried


def test_variant_parquet_generator_process_empty_input(temp_dir):
    """Worker with 0 reads writes a valid empty shard (edge #6)."""
    variant_parquet_generator_process(
        [], lambda *a: {}, _MockArgs(), {}, ["ref1"], None, (0, None),
        0, temp_dir, is_paired=False,
    )
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    assert os.path.exists(shard)
    assert list(iter_aligned_shard(shard)) == []


def test_concurrency_multiple_workers(temp_dir):
    """N workers writing N shards concurrently (edge #5: no inter-worker locking)."""
    import multiprocessing as mp

    n_workers = 4
    reads_per_worker = 25
    procs = []
    for wid in range(n_workers):
        items = [(f"seq_{wid}_{i}", i + 1, None) for i in range(reads_per_worker)]
        p = mp.Process(
            target=variant_parquet_generator_process,
            args=(items, _mock_aligner, _MockArgs(), {}, ["ref1"], None, (0, None), wid, temp_dir),
            kwargs={"is_paired": False},
        )
        p.start()
        procs.append(p)
    for p in procs:
        p.join()
    total = 0
    for wid in range(n_workers):
        shard = os.path.join(temp_dir, f"aligned_{wid}.parquet")
        assert os.path.exists(shard), f"shard {wid} missing"
        count = sum(1 for _ in iter_aligned_shard(shard))
        assert count == reads_per_worker, f"shard {wid} has {count} rows, expected {reads_per_worker}"
        total += count
    assert total == n_workers * reads_per_worker


def test_shard_schema_uniformity(temp_dir):
    """All shards share the same schema (edge #7: scan_parquet concat works)."""
    for wid in range(3):
        items = [(f"seq_{i}", 1, None) for i in range(5)]
        variant_parquet_generator_process(
            items, _mock_aligner, _MockArgs(), {}, ["ref1"], None, (0, None),
            wid, temp_dir, is_paired=False,
        )
    # polars scan_parquet over the glob should concat cleanly
    import polars as pl
    df = pl.scan_parquet(os.path.join(temp_dir, "aligned_*.parquet")).select(
        pl.len()
    ).collect()
    assert df.item() == 15
