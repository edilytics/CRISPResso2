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


# -- Stage 3: streaming collapse (PR 5) --------------------------------------
#
# Parity strategy: each test builds the *reference* output by replicating the
# exact CRISPRessoCORE.main allele loop (3a re-key, 3b RC merge, 3c fan-out via
# get_allele_row) in pure Python on the same canned payloads, then asserts that
# VariantStore.collapse reproduces it. The reference implementation is small
# and mirrors the production code line-for-line, so a mismatch points directly
# at the collapse logic.

from CRISPResso2.storage import (
    collapse_aligned_shards,
)
import pyarrow as pa
import pyarrow.parquet as pq


def _shard_payload(ref_name="ref1", aln_seq="ATCGATCG--ATCGATCGAT",
                   aln_ref="ATCGATCGAAATCGATCGAT", classification="MODIFIED",
                   deletion_n=2, insertion_n=0, substitution_n=0):
    """A single-ref payload with the fields get_allele_row consumes."""
    sub = {
        "ref_name": ref_name,
        "aln_seq": aln_seq,
        "aln_ref": aln_ref,
        "aln_strand": "+",
        "classification": classification,
        "irregular_ends": False,
        "insertions_outside_window": 0,
        "deletions_outside_window": 0,
        "substitutions_outside_window": 0,
        "total_mods": deletion_n + substitution_n,
        "mods_in_window": deletion_n + substitution_n,
        "mods_outside_window": 0,
        "insertion_n": insertion_n,
        "deletion_n": deletion_n,
        "substitution_n": substitution_n,
        "ref_positions": [0, 1, 2, 3, 4, 5, 6, 7, -8, -8, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
        "all_insertion_positions": [],
        "all_insertion_left_positions": [],
        "insertion_positions": [],
        "insertion_coordinates": [],
        "insertion_sizes": [],
        "all_deletion_positions": [8, 9],
        "deletion_positions": [8, 9],
        "deletion_coordinates": [(8, 10)],
        "deletion_sizes": [2],
        "all_substitution_positions": [],
        "substitution_positions": [],
        "substitution_values": np.array([]),
    }
    return {
        "count": 1,
        "aln_ref_names": [ref_name],
        "aln_scores": [95.0],
        "best_match_score": 95.0,
        "class_name": ref_name + "_" + classification,
        "best_match_name": ref_name,
        "caching_is_ok": True,
        "variant_" + ref_name: sub,
    }


def _write_shard(path, rows):
    """rows: list of (read_key, count, payload_dict). Writes one aligned shard."""
    with AlignedShardWriter(path) as w:
        for read_key, count, payload in rows:
            w.write_row(payload_to_row(read_key, count, payload))


def _ref_allele_row(reference_name, variant_count, aln_ref_names_str,
                    aln_ref_scores_str, variant_payload, write_detailed):
    """Reference replica of CRISPRessoCORE.get_allele_row."""
    if write_detailed:
        return {
            "#Reads": variant_count,
            "Aligned_Sequence": variant_payload["aln_seq"],
            "Reference_Sequence": variant_payload["aln_ref"],
            "n_inserted": variant_payload["insertion_n"],
            "n_deleted": variant_payload["deletion_n"],
            "n_mutated": variant_payload["substitution_n"],
            "Reference_Name": reference_name,
            "Read_Status": variant_payload["classification"],
            "Aligned_Reference_Names": aln_ref_names_str,
            "Aligned_Reference_Scores": aln_ref_scores_str,
            "ref_positions": variant_payload["ref_positions"],
            "all_insertion_positions": variant_payload["all_insertion_positions"],
            "all_insertion_left_positions": variant_payload["all_insertion_left_positions"],
            "insertion_positions": variant_payload["insertion_positions"],
            "insertion_coordinates": variant_payload["insertion_coordinates"],
            "insertion_sizes": variant_payload["insertion_sizes"],
            "all_deletion_positions": variant_payload["all_deletion_positions"],
            "deletion_positions": variant_payload["deletion_positions"],
            "deletion_coordinates": variant_payload["deletion_coordinates"],
            "deletion_sizes": variant_payload["deletion_sizes"],
            "all_substitution_positions": variant_payload["all_substitution_positions"],
            "substitution_positions": variant_payload["substitution_positions"],
            "substitution_values": variant_payload["substitution_values"],
        }
    return {
        "#Reads": variant_count,
        "Aligned_Sequence": variant_payload["aln_seq"],
        "Reference_Sequence": variant_payload["aln_ref"],
        "n_inserted": variant_payload["insertion_n"],
        "n_deleted": variant_payload["deletion_n"],
        "n_mutated": variant_payload["substitution_n"],
        "Reference_Name": reference_name,
        "Read_Status": variant_payload["classification"],
        "Aligned_Reference_Names": aln_ref_names_str,
        "Aligned_Reference_Scores": aln_ref_scores_str,
        "ref_positions": variant_payload["ref_positions"],
    }


def _ref_collapse(variants, *, is_paired, discard_indel_reads, write_detailed):
    """Reference replica of the CRISPRessoCORE.main allele loop.

    ``variants``: list of (key, count, payload) in insertion order. For paired
    input the key is re-keyed by the primary aln_seq (3a); for single-read the
    key is used as-is. Then 3b RC merge + 3c fan-out, returning the same tuple
    shape as ``_collapse_fanout`` plus the sorted allele_rows.
    """
    rc = CRISPRessoShared.reverse_complement
    store = {}
    for key, count, payload in variants:
        if payload.get("best_match_score", 0) <= 0:
            continue
        if is_paired:
            new_key = payload["variant_" + payload["aln_ref_names"][0]]["aln_seq"]
        else:
            new_key = key
        if new_key in store:
            store[new_key]["count"] += count
        else:
            store[new_key] = {"count": count, "payload": payload}

    # 3b RC merge
    for key in list(store.keys()):
        rec = store[key]
        vc = rec["count"]
        if vc == 0:
            continue
        rck = rc(key)
        if rck in store and store[rck]["count"] > 0:
            vc += store[rck]["count"]
            store[rck]["count"] = 0
            rec["count"] = vc

    # 3c fan-out (always full dicts for the reference, like collapse does)
    rows = []
    n_total = 0
    class_counts = {}
    counts_total = {}
    counts_modified = {}
    counts_unmodified = {}
    counts_discarded = {}
    for key in store:
        rec = store[key]
        vc = rec["count"]
        if vc == 0:
            continue
        n_total += vc
        payload = rec["payload"]
        cn = payload.get("class_name")
        class_counts[cn] = class_counts.get(cn, 0) + vc
        arn = payload["aln_ref_names"]
        arn_str = "&".join(arn)
        ars = payload.get("aln_scores", [])
        ars_str = "&".join([str(x) for x in ars])
        if cn == "AMBIGUOUS":
            rows.append(_ref_allele_row("AMBIGUOUS_" + arn[0], vc, arn_str, ars_str,
                                        payload["variant_" + arn[0]], True))
            continue
        for ref_name in arn:
            vp = payload["variant_" + ref_name]
            if discard_indel_reads and (vp["deletion_n"] > 0 or vp["insertion_n"] > 0):
                counts_discarded[ref_name] = counts_discarded.get(ref_name, 0) + vc
                rows.append(_ref_allele_row("DISCARDED_" + arn[0], vc, arn_str, ars_str, vp, True))
                continue
            counts_total[ref_name] = counts_total.get(ref_name, 0) + vc
            if vp["classification"] == "MODIFIED":
                counts_modified[ref_name] = counts_modified.get(ref_name, 0) + vc
            else:
                counts_unmodified[ref_name] = counts_unmodified.get(ref_name, 0) + vc
            rows.append(_ref_allele_row(ref_name, vc, arn_str, ars_str, vp, True))
    rows.sort(key=lambda r: (-int(r["#Reads"]), r["Aligned_Sequence"], r["Reference_Sequence"]))
    if not write_detailed:
        non_det = ("#Reads", "Aligned_Sequence", "Reference_Sequence", "n_inserted",
                   "n_deleted", "n_mutated", "Reference_Name", "Read_Status",
                   "Aligned_Reference_Names", "Aligned_Reference_Scores", "ref_positions")
        rows = [{k: r[k] for k in non_det} for r in rows]
    return rows, n_total, class_counts, counts_total, counts_modified, counts_unmodified, counts_discarded


def _assert_rows_equal(a, b):
    """Type-tolerant allele-row equality (numpy arrays compared by value)."""
    assert len(a) == len(b), f"row count {len(a)} != {len(b)}"
    for i, (ra, rb) in enumerate(zip(a, b)):
        assert set(ra.keys()) == set(rb.keys()), f"row {i} keys {set(ra)} != {set(rb)}"
        for k in ra:
            va, vb = ra[k], rb[k]
            if isinstance(va, np.ndarray) or isinstance(vb, np.ndarray):
                va_l = va.tolist() if isinstance(va, np.ndarray) else list(va)
                vb_l = vb.tolist() if isinstance(vb, np.ndarray) else list(vb)
                assert va_l == vb_l, f"row {i} key {k}: {va_l!r} != {vb_l!r}"
            else:
                assert va == vb, f"row {i} key {k}: {va!r} != {vb!r}"


# -- 3a: aln_seq re-key (paired) --------------------------------------------


def test_3a_rekey_merges_same_aln_seq_paired(temp_dir):
    """Two different paired keys aligning to the same aln_seq → 1 row, counts summed."""
    p1 = _shard_payload()
    p2 = _shard_payload()  # identical payload -> same aln_seq
    rows = [("AAA+TTT", 3, p1), ("GGG+CCC", 2, p2)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=True)
    assert len(res.allele_rows) == 1
    assert res.allele_rows[0]["#Reads"] == 5
    assert res.n_total == 5


def test_3a_no_rekey_for_single_read(temp_dir):
    """Single-read: keys are raw seqs (no '+'); no re-key. Two distinct seqs → 2 rows."""
    p1 = _shard_payload()
    p2 = _shard_payload()
    # different aln_seq so they don't RC-merge; both read_keys are non-palindromic
    # and not reverse complements of each other (so 3b leaves them separate).
    p2["variant_ref1"]["aln_seq"] = "AAACCCGT--AAACCCGTAA"
    p2["variant_ref1"]["aln_ref"] = "AAACCCGTAAACCCGTAA"
    rows = [("AAAACCCC", 3, p1), ("AAACCCGT", 2, p2)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False)
    assert len(res.allele_rows) == 2
    assert {r["#Reads"] for r in res.allele_rows} == {3, 2}
    assert res.n_total == 5


def test_3a_rekey_parity_vs_reference(temp_dir):
    """3a re-key output matches the reference loop on canned payloads."""
    p1 = _shard_payload()
    p2 = _shard_payload()
    p3 = _shard_payload()
    p3["variant_ref1"]["aln_seq"] = "TTTTAAAA--TTTTAAAATT"
    rows = [("AAA+TTT", 3, p1), ("GGG+CCC", 4, p2), ("CCC+GGG", 1, p3)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=True, write_detailed_allele_table=True)
    ref = _ref_collapse(rows, is_paired=True, discard_indel_reads=False, write_detailed=True)
    _assert_rows_equal(res.allele_rows, ref[0])
    assert res.n_total == ref[1]
    assert res.class_counts == ref[2]


# -- 3b: reverse-complement merge -------------------------------------------


def test_3b_rc_merge_collapses_strand_pair(temp_dir):
    """A read and its reverse complement (as keys) merge → 1 row, counts summed."""
    p1 = _shard_payload(aln_seq="AAAATTTTGGGG", aln_ref="AAAATTTTGGGG", classification="UNMODIFIED",
                        deletion_n=0)
    p1["class_name"] = "ref1_UNMODIFIED"
    p1["variant_ref1"]["ref_positions"] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    p1["variant_ref1"]["all_deletion_positions"] = []
    p1["variant_ref1"]["deletion_positions"] = []
    p1["variant_ref1"]["deletion_coordinates"] = []
    p1["variant_ref1"]["deletion_sizes"] = []
    # RC of "AAAATTTTGGGG" is "CCCCAAAATTTT"
    p2 = _shard_payload(aln_seq="CCCCAAAATTTT", aln_ref="CCCCAAAATTTT", classification="UNMODIFIED",
                        deletion_n=0)
    p2["class_name"] = "ref1_UNMODIFIED"
    p2["variant_ref1"]["ref_positions"] = list(range(12))
    p2["variant_ref1"]["all_deletion_positions"] = []
    p2["variant_ref1"]["deletion_positions"] = []
    p2["variant_ref1"]["deletion_coordinates"] = []
    p2["variant_ref1"]["deletion_sizes"] = []
    # single-read: key = raw seq = aln_seq here
    rows = [("AAAATTTTGGGG", 3, p1), ("CCCCAAAATTTT", 2, p2)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False)
    assert len(res.allele_rows) == 1
    assert res.allele_rows[0]["#Reads"] == 5
    assert res.n_total == 5


def test_3b_rc_merge_representative_is_first_in_scan_order(temp_dir):
    """The first-iterated key of an RC pair is the representative (pandas tie-break)."""
    p1 = _shard_payload(aln_seq="AAAATTTT", aln_ref="AAAATTTT", classification="UNMODIFIED", deletion_n=0)
    p1["class_name"] = "ref1_UNMODIFIED"
    p1["variant_ref1"]["ref_positions"] = list(range(8))
    p1["variant_ref1"]["all_deletion_positions"] = []
    p1["variant_ref1"]["deletion_positions"] = []
    p1["variant_ref1"]["deletion_coordinates"] = []
    p1["variant_ref1"]["deletion_sizes"] = []
    # RC of AAAATTTT is AAATTTTT (AAAA->TTTT rc, TTTT->AAAA; reverse) = "AAAATTTT"? let's use a non-palindrome
    # "AAAAGCCC" rc = "GGGCTTTT"
    p1["variant_ref1"]["aln_seq"] = "AAAAGCCC"
    p1["variant_ref1"]["aln_ref"] = "AAAAGCCC"
    p2 = _shard_payload(aln_seq="GGGCTTTT", aln_ref="GGGCTTTT", classification="UNMODIFIED", deletion_n=0)
    p2["class_name"] = "ref1_UNMODIFIED"
    p2["variant_ref1"]["ref_positions"] = list(range(8))
    p2["variant_ref1"]["all_deletion_positions"] = []
    p2["variant_ref1"]["deletion_positions"] = []
    p2["variant_ref1"]["deletion_coordinates"] = []
    p2["variant_ref1"]["deletion_sizes"] = []
    rows = [("AAAAGCCC", 3, p1), ("GGGCTTTT", 2, p2)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False)
    assert len(res.allele_rows) == 1
    # representative = first in scan order = "AAAAGCCC" (p1's aln_seq)
    assert res.allele_rows[0]["Aligned_Sequence"] == "AAAAGCCC"
    assert res.allele_rows[0]["#Reads"] == 5


def test_3b_rc_merge_parity_vs_reference(temp_dir):
    """3b RC merge matches the reference loop (insertion-order tie-break)."""
    p1 = _shard_payload()
    p1["variant_ref1"]["aln_seq"] = "AAAAGCCC"
    p1["variant_ref1"]["aln_ref"] = "AAAAGCCC"
    p1["variant_ref1"]["ref_positions"] = list(range(8))
    p1["variant_ref1"]["all_deletion_positions"] = []
    p1["variant_ref1"]["deletion_positions"] = []
    p1["variant_ref1"]["deletion_coordinates"] = []
    p1["variant_ref1"]["deletion_sizes"] = []
    p1["variant_ref1"]["deletion_n"] = 0
    p1["class_name"] = "ref1_UNMODIFIED"
    p2 = _shard_payload()
    p2["variant_ref1"]["aln_seq"] = "GGGCTTTT"  # rc of AAAAGCCC
    p2["variant_ref1"]["aln_ref"] = "GGGCTTTT"
    p2["variant_ref1"]["ref_positions"] = list(range(8))
    p2["variant_ref1"]["all_deletion_positions"] = []
    p2["variant_ref1"]["deletion_positions"] = []
    p2["variant_ref1"]["deletion_coordinates"] = []
    p2["variant_ref1"]["deletion_sizes"] = []
    p2["variant_ref1"]["deletion_n"] = 0
    p2["class_name"] = "ref1_UNMODIFIED"
    p3 = _shard_payload()  # a third, unrelated allele
    p3["variant_ref1"]["aln_seq"] = "CCCCGGGG"
    p3["variant_ref1"]["aln_ref"] = "CCCCGGGG"
    p3["variant_ref1"]["ref_positions"] = list(range(8))
    p3["variant_ref1"]["all_deletion_positions"] = []
    p3["variant_ref1"]["deletion_positions"] = []
    p3["variant_ref1"]["deletion_coordinates"] = []
    p3["variant_ref1"]["deletion_sizes"] = []
    p3["variant_ref1"]["deletion_n"] = 0
    p3["class_name"] = "ref1_UNMODIFIED"
    rows = [("AAAAGCCC", 3, p1), ("GGGCTTTT", 2, p2), ("CCCCGGGG", 4, p3)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False, write_detailed_allele_table=True)
    ref = _ref_collapse(rows, is_paired=False, discard_indel_reads=False, write_detailed=True)
    _assert_rows_equal(res.allele_rows, ref[0])
    assert res.n_total == ref[1]


def test_3b_palindrome_doubles_count_parity(temp_dir):
    """A palindromic key (rc(key)==key) doubles its count — pandas behaviour, replicated."""
    # "ATAT" rc = "ATAT" (palindrome)
    p = _shard_payload(aln_seq="ATAT", aln_ref="ATAT", classification="UNMODIFIED", deletion_n=0)
    p["class_name"] = "ref1_UNMODIFIED"
    p["variant_ref1"]["ref_positions"] = [0, 1, 2, 3]
    p["variant_ref1"]["all_deletion_positions"] = []
    p["variant_ref1"]["deletion_positions"] = []
    p["variant_ref1"]["deletion_coordinates"] = []
    p["variant_ref1"]["deletion_sizes"] = []
    rows = [("ATAT", 5, p)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False, write_detailed_allele_table=True)
    ref = _ref_collapse(rows, is_paired=False, discard_indel_reads=False, write_detailed=True)
    # pandas doubles: 5 -> 10
    assert res.allele_rows[0]["#Reads"] == 10
    assert res.n_total == 10
    _assert_rows_equal(res.allele_rows, ref[0])


# -- 3c: multi-reference fan-out --------------------------------------------


def test_3c_multi_ref_fanout_two_rows(temp_dir):
    """A read aligning to 2 refs (expand_ambiguous) → 2 allele rows."""
    sub1 = _shard_payload("ref1")["variant_ref1"]
    sub2 = _shard_payload("ref2")["variant_ref2"]
    sub2["ref_name"] = "ref2"
    payload = {
        "count": 1, "aln_ref_names": ["ref1", "ref2"], "aln_scores": [95.0, 95.0],
        "best_match_score": 95.0, "class_name": "ref1_MODIFIED&ref2_MODIFIED",
        "best_match_name": "ref1", "caching_is_ok": True,
        "variant_ref1": sub1, "variant_ref2": sub2,
    }
    rows = [("KEY+YYY", 7, payload)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=True)
    assert len(res.allele_rows) == 2
    assert {r["Reference_Name"] for r in res.allele_rows} == {"ref1", "ref2"}
    assert all(r["#Reads"] == 7 for r in res.allele_rows)
    assert res.counts_total == {"ref1": 7, "ref2": 7}


def test_3c_ambiguous_one_row(temp_dir):
    """class_name == 'AMBIGUOUS' → one row 'AMBIGUOUS_<ref0>', no per-ref fan-out."""
    sub1 = _shard_payload("ref1")["variant_ref1"]
    sub2 = _shard_payload("ref2")["variant_ref2"]
    sub2["ref_name"] = "ref2"
    payload = {
        "count": 1, "aln_ref_names": ["ref1", "ref2"], "aln_scores": [95.0, 95.0],
        "best_match_score": 95.0, "class_name": "AMBIGUOUS",
        "best_match_name": "ref1", "caching_is_ok": True,
        "variant_ref1": sub1, "variant_ref2": sub2,
    }
    rows = [("KEY+YYY", 9, payload)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=True)
    assert len(res.allele_rows) == 1
    assert res.allele_rows[0]["Reference_Name"] == "AMBIGUOUS_ref1"
    assert res.allele_rows[0]["#Reads"] == 9
    # ambiguous reads do NOT contribute to counts_total
    assert res.counts_total == {}
    assert res.class_counts == {"AMBIGUOUS": 9}


def test_3c_discarded_indel_reads(temp_dir):
    """discard_indel_reads: indel-bearing refs → DISCARDED_<ref0> rows, not counts_total."""
    sub1 = _shard_payload("ref1", deletion_n=2)  # has deletion -> discarded
    p1 = sub1  # single ref
    payload = {
        "count": 1, "aln_ref_names": ["ref1"], "aln_scores": [95.0],
        "best_match_score": 95.0, "class_name": "ref1_MODIFIED",
        "best_match_name": "ref1", "caching_is_ok": True,
        "variant_ref1": sub1["variant_ref1"],
    }
    rows = [("ATCGATCG", 6, payload)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False, discard_indel_reads=True)
    assert len(res.allele_rows) == 1
    assert res.allele_rows[0]["Reference_Name"] == "DISCARDED_ref1"
    assert res.allele_rows[0]["#Reads"] == 6
    assert res.counts_total == {}
    assert res.counts_discarded == {"ref1": 6}


def test_3c_discarded_per_ref_when_multi_ref(temp_dir):
    """discard_indel_reads with multi-ref: a DISCARDED row per ref, counts_discarded keyed by actual ref."""
    sub1 = _shard_payload("ref1", deletion_n=2)["variant_ref1"]  # has deletion -> discarded
    sub2 = _shard_payload("ref2", deletion_n=0)["variant_ref2"]  # no deletion -> kept
    sub2["ref_name"] = "ref2"
    sub2["classification"] = "UNMODIFIED"  # match class_name suffix
    payload = {
        "count": 1, "aln_ref_names": ["ref1", "ref2"], "aln_scores": [95.0, 95.0],
        "best_match_score": 95.0, "class_name": "ref1_MODIFIED&ref2_UNMODIFIED",
        "best_match_name": "ref1", "caching_is_ok": True,
        "variant_ref1": sub1, "variant_ref2": sub2,
    }
    rows = [("KEY+YYY", 4, payload)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=True, discard_indel_reads=True)
    # ref1 discarded (DISCARDED_ref1), ref2 kept (ref2)
    names = sorted(r["Reference_Name"] for r in res.allele_rows)
    assert names == ["DISCARDED_ref1", "ref2"]
    assert res.counts_discarded == {"ref1": 4}
    assert res.counts_total == {"ref2": 4}
    assert res.counts_unmodified == {"ref2": 4}


def test_3c_fanout_parity_vs_reference(temp_dir):
    """3c fan-out matches the reference loop across normal/ambiguous/discard cases."""
    sub1 = _shard_payload("ref1", deletion_n=2)["variant_ref1"]
    sub2 = _shard_payload("ref2", deletion_n=0)["variant_ref2"]
    sub2["ref_name"] = "ref2"
    amb_payload = {
        "count": 1, "aln_ref_names": ["ref1", "ref2"], "aln_scores": [88.0, 88.0],
        "best_match_score": 88.0, "class_name": "AMBIGUOUS",
        "best_match_name": "ref1", "caching_is_ok": True,
        "variant_ref1": dict(_shard_payload("ref1")["variant_ref1"]),
        "variant_ref2": dict(_shard_payload("ref2")["variant_ref2"]),
    }
    amb_payload["variant_ref2"]["ref_name"] = "ref2"
    norm_payload = {
        "count": 1, "aln_ref_names": ["ref1", "ref2"], "aln_scores": [95.0, 95.0],
        "best_match_score": 95.0, "class_name": "ref1_MODIFIED&ref2_UNMODIFIED",
        "best_match_name": "ref1", "caching_is_ok": True,
        "variant_ref1": sub1, "variant_ref2": sub2,
    }
    rows = [("AMB+XXX", 9, amb_payload), ("NRM+YYY", 4, norm_payload)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    for discard in (False, True):
        for detailed in (False, True):
            res = collapse_aligned_shards(
                [shard], temp_dir, is_paired=True,
                discard_indel_reads=discard, write_detailed_allele_table=detailed,
            )
            ref = _ref_collapse(rows, is_paired=True, discard_indel_reads=discard,
                                write_detailed=detailed)
            _assert_rows_equal(res.allele_rows, ref[0])
            assert res.n_total == ref[1]
            assert res.class_counts == ref[2]
            assert res.counts_total == ref[3]
            assert res.counts_discarded == ref[6]


# -- output order + dtype + artifact ----------------------------------------


def test_output_sorted_by_reads_desc_then_aligned_then_ref(temp_dir):
    """Allele rows sorted (#Reads desc, Aligned_Sequence asc, Reference_Sequence asc)."""
    p1 = _shard_payload()
    p1["variant_ref1"]["aln_seq"] = "AAAC"
    p1["variant_ref1"]["aln_ref"] = "AAAC"
    p2 = _shard_payload()
    p2["variant_ref1"]["aln_seq"] = "AACA"
    p2["variant_ref1"]["aln_ref"] = "AACA"
    p3 = _shard_payload()
    p3["variant_ref1"]["aln_seq"] = "ACAA"
    p3["variant_ref1"]["aln_ref"] = "ACAA"
    for p in (p1, p2, p3):
        p["variant_ref1"]["ref_positions"] = list(range(4))
        p["variant_ref1"]["all_deletion_positions"] = []
        p["variant_ref1"]["deletion_positions"] = []
        p["variant_ref1"]["deletion_coordinates"] = []
        p["variant_ref1"]["deletion_sizes"] = []
        p["variant_ref1"]["deletion_n"] = 0
        p["class_name"] = "ref1_MODIFIED"
    rows = [("AAAC", 1, p1), ("AACA", 3, p2), ("ACAA", 2, p3)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False)
    seqs = [r["Aligned_Sequence"] for r in res.allele_rows]
    counts = [r["#Reads"] for r in res.allele_rows]
    assert counts == sorted(counts, reverse=True)  # desc by #Reads
    # counts 3,2,1 -> AACA, ACAA, AAAC
    assert seqs == ["AACA", "ACAA", "AAAC"]


def test_dtype_count_int64_positions_int64(temp_dir):
    """Collapsed parquet stores #Reads and position arrays as int64 (edge #17)."""
    p = _shard_payload()
    rows = [("ATCGATCG", 5, p)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False)
    assert res.parquet_path and os.path.exists(res.parquet_path)
    pf = pq.ParquetFile(res.parquet_path)
    schema = pf.schema_arrow
    assert schema.field("#Reads").type == pa.int64()
    assert schema.field("ref_positions").type == pa.list_(pa.int64())
    assert schema.field("all_deletion_positions").type == pa.list_(pa.int64())
    # coordinates are list of struct{start,end:int64} (edge #18)
    coord_t = schema.field("deletion_coordinates").type
    assert pa.types.is_list(coord_t)
    struct_t = coord_t.value_type
    assert pa.types.is_struct(struct_t)
    assert struct_t.field("start").type == pa.int64()
    assert struct_t.field("end").type == pa.int64()
    # substitution_values is list<string>
    assert schema.field("substitution_values").type == pa.list_(pa.string())
    # row count
    assert pf.read().num_rows == 1


def test_dtype_substitution_values_round_trip(temp_dir):
    """substitution_values round-trips as list<string> and reads back correctly."""
    p = _shard_payload()
    p["variant_ref1"]["substitution_values"] = np.array(["G", "T"])
    p["variant_ref1"]["substitution_n"] = 2
    p["variant_ref1"]["substitution_positions"] = [4, 7]
    p["variant_ref1"]["all_substitution_positions"] = [4, 7]
    rows = [("ATCGATCG", 5, p)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False, write_detailed_allele_table=True)
    # read back the parquet
    tbl = pq.read_table(res.parquet_path)
    row = tbl.to_pylist()[0]
    assert row["substitution_values"] == ["G", "T"]
    # and the returned allele row keeps the numpy array (parity with pandas get_allele_row)
    assert isinstance(res.allele_rows[0]["substitution_values"], np.ndarray)
    assert res.allele_rows[0]["substitution_values"].tolist() == ["G", "T"]


def test_non_detailed_allele_rows_subset(temp_dir):
    """Non-detailed allele_rows have only the 11-key subset; parquet still full."""
    p = _shard_payload()
    rows = [("ATCGATCG", 5, p)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False, write_detailed_allele_table=False)
    keys = set(res.allele_rows[0].keys())
    expected = {"#Reads", "Aligned_Sequence", "Reference_Sequence", "n_inserted",
                "n_deleted", "n_mutated", "Reference_Name", "Read_Status",
                "Aligned_Reference_Names", "Aligned_Reference_Scores", "ref_positions"}
    assert keys == expected
    # parquet is the full 23-column schema regardless
    schema = pq.ParquetFile(res.parquet_path).schema_arrow
    assert "all_insertion_positions" in schema.names
    assert "substitution_values" in schema.names


def test_detailed_allele_rows_full(temp_dir):
    """Detailed allele_rows have all 23 keys."""
    p = _shard_payload()
    rows = [("ATCGATCG", 5, p)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False, write_detailed_allele_table=True)
    assert len(res.allele_rows[0]) == 23


def test_unaligned_reads_filtered(temp_dir):
    """best_match_score <= 0 reads are excluded from collapse (removed pre-loop in pandas)."""
    p_aln = _shard_payload()
    p_unaln = {"count": 1, "aln_scores": [3.0], "best_match_score": 0}
    rows = [("ATCGATCG", 5, p_aln), ("NNNNNNNN", 2, p_unaln)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False)
    assert len(res.allele_rows) == 1
    assert res.allele_rows[0]["#Reads"] == 5
    assert res.n_total == 5


def test_empty_input(temp_dir):
    """No aligned reads → empty allele table, n_total=0, valid empty parquet."""
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    with AlignedShardWriter(shard) as w:
        pass
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False)
    assert res.allele_rows == []
    assert res.n_total == 0
    assert res.class_counts == {}
    assert os.path.exists(res.parquet_path)
    assert pq.ParquetFile(res.parquet_path).read().num_rows == 0


def test_glob_string_shard_paths(temp_dir):
    """Collapse accepts a glob string for shard_paths."""
    for i in range(3):
        p = _shard_payload()
        seq = ["AAAA", "AACC", "ACCA"][i]  # none palindromic, none RC pairs
        p["variant_ref1"]["aln_seq"] = seq
        p["variant_ref1"]["aln_ref"] = seq
        p["variant_ref1"]["ref_positions"] = list(range(4))
        p["variant_ref1"]["all_deletion_positions"] = []
        p["variant_ref1"]["deletion_positions"] = []
        p["variant_ref1"]["deletion_coordinates"] = []
        p["variant_ref1"]["deletion_sizes"] = []
        p["variant_ref1"]["deletion_n"] = 0
        p["class_name"] = "ref1_UNMODIFIED"
        p["variant_ref1"]["classification"] = "UNMODIFIED"
        rows = [(seq, i + 1, p)]
        _write_shard(os.path.join(temp_dir, f"aligned_{i}.parquet"), rows)
    res = collapse_aligned_shards(
        os.path.join(temp_dir, "aligned_*.parquet"), temp_dir, is_paired=False
    )
    assert len(res.allele_rows) == 3
    assert res.n_total == 1 + 2 + 3


def test_multi_shard_concat(temp_dir):
    """Multiple shards are concatenated by the collapse (one logical input)."""
    p1 = _shard_payload()
    p2 = _shard_payload()  # same aln_seq -> should merge (paired) to 1 row, count 7
    _write_shard(os.path.join(temp_dir, "aligned_0.parquet"), [("AAA+TTT", 3, p1)])
    _write_shard(os.path.join(temp_dir, "aligned_1.parquet"), [("GGG+CCC", 4, p2)])
    res = collapse_aligned_shards(
        [os.path.join(temp_dir, "aligned_0.parquet"), os.path.join(temp_dir, "aligned_1.parquet")],
        temp_dir, is_paired=True,
    )
    assert len(res.allele_rows) == 1
    assert res.allele_rows[0]["#Reads"] == 7


def test_allele_rows_dataframe_with_pct_reads(temp_dir):
    """allele_rows_dataframe adds %Reads and casts n_* to int (parity with df_alleles)."""
    p1 = _shard_payload()
    p1["variant_ref1"]["aln_seq"] = "AAAA"
    p2 = _shard_payload()
    p2["variant_ref1"]["aln_seq"] = "CCCC"
    for p in (p1, p2):
        p["variant_ref1"]["aln_ref"] = p["variant_ref1"]["aln_seq"]
        p["variant_ref1"]["ref_positions"] = list(range(4))
        p["variant_ref1"]["all_deletion_positions"] = []
        p["variant_ref1"]["deletion_positions"] = []
        p["variant_ref1"]["deletion_coordinates"] = []
        p["variant_ref1"]["deletion_sizes"] = []
        p["variant_ref1"]["deletion_n"] = 0
        p["class_name"] = "ref1_MODIFIED"
    rows = [("AAAA", 3, p1), ("CCCC", 1, p2)]
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    res = collapse_aligned_shards([shard], temp_dir, is_paired=False)
    df = res.allele_rows_dataframe()
    assert "%Reads" in df.columns
    # 3/4*100=75, 1/4*100=25
    assert sorted(df["%Reads"].tolist()) == [25.0, 75.0]
    assert df["n_deleted"].dtype.kind in "iu"


def test_full_parity_fuzz_vs_reference(temp_dir):
    """Fuzz: random canned payloads through collapse vs the reference loop."""
    import random
    rng = random.Random(123)
    bases = "ACGT"
    rows = []
    for i in range(40):
        ref = rng.choice(["ref1", "ref2"])
        seq = "".join(rng.choice(bases) for _ in range(8))
        # ~50% make a second ref (ambiguous or expanded)
        if rng.random() < 0.3:
            refs = ["ref1", "ref2"]
            cn = rng.choice(["AMBIGUOUS", "ref1_MODIFIED&ref2_UNMODIFIED"])
        else:
            refs = [ref]
            cn = ref + "_" + rng.choice(["MODIFIED", "UNMODIFIED"])
        payload = {
            "count": 1, "aln_ref_names": refs,
            "aln_scores": [float(rng.randint(80, 99)) for _ in refs],
            "best_match_score": float(rng.randint(80, 99)),
            "class_name": cn, "best_match_name": refs[0], "caching_is_ok": True,
        }
        for r in refs:
            sub = _shard_payload(r)["variant_" + r]
            sub["aln_seq"] = seq
            sub["aln_ref"] = seq
            sub["ref_positions"] = list(range(8))
            sub["all_deletion_positions"] = []
            sub["deletion_positions"] = []
            sub["deletion_coordinates"] = []
            sub["deletion_sizes"] = []
            sub["all_insertion_positions"] = []
            sub["insertion_positions"] = []
            sub["insertion_coordinates"] = []
            sub["insertion_sizes"] = []
            sub["all_substitution_positions"] = []
            sub["substitution_positions"] = []
            sub["substitution_values"] = np.array([])
            sub["deletion_n"] = 0
            sub["insertion_n"] = 0
            sub["substitution_n"] = 0
            sub["classification"] = "UNMODIFIED" if cn.endswith("UNMODIFIED") or cn == "AMBIGUOUS" else "MODIFIED"
            payload["variant_" + r] = sub
        rows.append((seq, rng.randint(1, 5), payload))
    shard = os.path.join(temp_dir, "aligned_0.parquet")
    _write_shard(shard, rows)
    for is_paired in (True, False):
        for discard in (False, True):
            for detailed in (False, True):
                res = collapse_aligned_shards(
                    [shard], temp_dir, is_paired=is_paired,
                    discard_indel_reads=discard, write_detailed_allele_table=detailed,
                )
                ref = _ref_collapse(rows, is_paired=is_paired, discard_indel_reads=discard,
                                    write_detailed=detailed)
                _assert_rows_equal(res.allele_rows, ref[0])
                assert res.n_total == ref[1]
                assert res.class_counts == ref[2]
                assert res.counts_total == ref[3]
                assert res.counts_modified == ref[4]
                assert res.counts_unmodified == ref[5]
                assert res.counts_discarded == ref[6]
