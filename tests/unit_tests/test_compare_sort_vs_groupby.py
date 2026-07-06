"""Unit tests for the sort vs group_by comparison harness (S1b)."""
import importlib.util
from pathlib import Path

_HARNESS_PATH = Path(__file__).resolve().parents[2] / "scripts" / "compare_sort_vs_groupby.py"
_spec = importlib.util.spec_from_file_location("compare_sort_vs_groupby", _HARNESS_PATH)
cmp_harness = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cmp_harness)


def test_random_seqs_chunk_length_and_alphabet():
    seqs = cmp_harness._random_seqs_chunk(length=40, n=4, seed=1)
    assert len(seqs) == 4
    for s in seqs:
        assert len(s) == 40
        assert set(s) <= set("ACGT")


def test_random_seqs_chunk_deterministic():
    a = cmp_harness._random_seqs_chunk(length=30, n=3, seed=9)
    b = cmp_harness._random_seqs_chunk(length=30, n=3, seed=9)
    assert a == b


def test_generate_reads_parquet_row_count(temp_dir):
    import polars as pl

    shards = cmp_harness.generate_reads_parquet(
        Path(temp_dir), n_rows=250, read_length=15, chunk_size=100,
    )
    assert len(shards) == 3
    total = sum(pl.read_parquet(s).height for s in shards)
    assert total == 250


def test_peak_rss_self_positive():
    assert cmp_harness._peak_rss_self_mb() > 0


def test_peak_rss_children_nonnegative():
    # No children run yet, but the call must not error.
    assert cmp_harness._peak_rss_children_mb() >= 0


def test_format_report_rends_results():
    results = [
        {"engine": "polars_gb", "n_rows": 1000, "read_length": 100, "peak_rss_mb": 80.0,
         "self_rss_mb": 80.0, "child_rss_mb": 0.0, "wall_seconds": 0.1, "n_groups": 1000, "ok": True},
        {"engine": "sort_uniq", "n_rows": 1000, "read_length": 100, "peak_rss_mb": 60.0,
         "self_rss_mb": 60.0, "child_rss_mb": 5.0, "wall_seconds": 0.2, "n_groups": 1000, "ok": True},
    ]
    report = cmp_harness._format_report(results, ["polars_gb", "sort_uniq"], [1000], [100])
    assert "sort | uniq -c vs polars group_by" in report
    assert "polars_gb" in report
    assert "sort_uniq" in report


def test_format_report_handles_failed_row():
    results = [{"engine": "polars_gb", "n_rows": 1000, "read_length": 100, "ok": False, "error": "boom"}]
    report = cmp_harness._format_report(results, ["polars_gb"], [1000], [100])
    assert "FAIL" in report
    assert "boom" in report


def test_format_report_flat_scaling_verdict():
    results = [
        {"engine": "sort_uniq", "read_length": 1000, "n_rows": 10000, "peak_rss_mb": 90.0, "ok": True},
        {"engine": "sort_uniq", "read_length": 1000, "n_rows": 1000000, "peak_rss_mb": 190.0, "ok": True},
        {"engine": "polars_gb", "read_length": 1000, "n_rows": 10000, "peak_rss_mb": 137.0, "ok": True},
        {"engine": "polars_gb", "read_length": 1000, "n_rows": 1000000, "peak_rss_mb": 3017.0, "ok": True},
    ]
    report = cmp_harness._format_report(results, ["sort_uniq", "polars_gb"], [10000, 1000000], [1000])
    assert "[FLAT]" in report
    assert "[SCALES]" in report


def test_worker_rejects_unknown_engine(tmp_path):
    import argparse

    (tmp_path / "reads_00000.parquet").write_bytes(b"")
    ns = argparse.Namespace(
        _worker=True, _workdir=str(tmp_path), _engine="bogus",
        _rows=1, _read_length=10, _unique_fraction=1.0,
    )
    assert cmp_harness._worker(ns) == 1
