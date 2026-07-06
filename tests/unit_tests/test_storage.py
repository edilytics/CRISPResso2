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

import pytest

from CRISPResso2 import CRISPRessoShared
from CRISPResso2.storage import (
    DEFAULT_MEMORY_BUDGET_MB,
    ReadCounts,
    VariantStore,
    count_reads_from_fastq,
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
