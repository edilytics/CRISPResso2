"""Unit tests for the streaming memory benchmark harness (PR 2 / Spike S1).

These tests cover the harness plumbing (generation, gate logic, formatting)
without running the expensive full matrix. The actual RSS measurements are
validated by running the harness end-to-end (see design_docs/STREAMING_GROUPBY_SPIKE.md).
"""
import importlib.util
from pathlib import Path


# The harness lives in scripts/, not the package. Import it by path so tests
# don't require it to be installed.
_HARNESS_PATH = Path(__file__).resolve().parents[2] / "scripts" / "bench_streaming.py"
_spec = importlib.util.spec_from_file_location("bench_streaming", _HARNESS_PATH)
bench = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bench)


def test_random_seqs_chunk_produces_correct_length_and_alphabet():
    """Generated reads must be the requested length and only contain ACGT."""
    seqs = bench._random_seqs_chunk(length=50, n=5, seed=1)
    assert len(seqs) == 5
    for s in seqs:
        assert len(s) == 50
        assert set(s) <= set("ACGT")


def test_random_seqs_chunk_is_deterministic_given_seed():
    """Same seed must produce the same sequences (reproducibility)."""
    a = bench._random_seqs_chunk(length=100, n=3, seed=42)
    b = bench._random_seqs_chunk(length=100, n=3, seed=42)
    assert a == b


def test_random_seqs_chunk_different_seeds_differ():
    a = bench._random_seqs_chunk(length=100, n=3, seed=1)
    b = bench._random_seqs_chunk(length=100, n=3, seed=2)
    assert a != b


def test_generate_reads_parquet_writes_shards_with_correct_row_count(temp_dir):
    """Generation must produce parquet shards whose total rows == n_rows."""
    import polars as pl

    shards = bench.generate_reads_parquet(
        Path(temp_dir), n_rows=250, read_length=20, unique_fraction=1.0, chunk_size=100,
    )
    # 250 rows / 100 per chunk => 3 shards
    assert len(shards) == 3
    total = 0
    for shard in shards:
        df = pl.read_parquet(shard)
        assert set(df.columns) == {"read_key", "count"}
        total += df.height
    assert total == 250


def test_generate_reads_parquet_with_duplicates_collapses_pool(temp_dir):
    """With unique_fraction < 1.0, the pool has fewer unique keys than rows."""
    import polars as pl

    shards = bench.generate_reads_parquet(
        Path(temp_dir), n_rows=1000, read_length=20, unique_fraction=0.1, chunk_size=500,
    )
    all_keys = []
    for shard in shards:
        all_keys.extend(pl.read_parquet(shard)["read_key"].to_list())
    assert len(all_keys) == 1000
    # 10% unique => at most ~100 distinct keys (sampling may repeat, so <= 100)
    assert len(set(all_keys)) <= 100


def test_gate_verdict_pass_when_rss_flat():
    """If streaming RSS stays flat, the gate passes."""
    results = [
        {"engine": "streaming", "read_length": 1000, "n_rows": 10000, "peak_rss_mb": 100.0, "ok": True},
        {"engine": "streaming", "read_length": 1000, "n_rows": 1000000, "peak_rss_mb": 250.0, "ok": True},
    ]
    overall, detail = bench._gate_verdict(results)
    assert overall == "PASS"
    assert "streaming" in detail


def test_gate_verdict_fail_when_rss_scales():
    """If every spilling engine's RSS scales linearly, the gate fails."""
    results = [
        {"engine": "streaming", "read_length": 1000, "n_rows": 10000, "peak_rss_mb": 100.0, "ok": True},
        {"engine": "streaming", "read_length": 1000, "n_rows": 1000000, "peak_rss_mb": 5000.0, "ok": True},
        {"engine": "duckdb", "read_length": 1000, "n_rows": 10000, "peak_rss_mb": 120.0, "ok": True},
        {"engine": "duckdb", "read_length": 1000, "n_rows": 1000000, "peak_rss_mb": 4000.0, "ok": True},
    ]
    overall, detail = bench._gate_verdict(results)
    assert overall == "FAIL"
    assert "streaming" in detail
    assert "duckdb" in detail


def test_gate_verdict_pass_if_duckdb_flat_even_when_streaming_scales():
    """The gate passes if ANY spilling engine stays flat — that engine is the chosen backend."""
    results = [
        {"engine": "streaming", "read_length": 1000, "n_rows": 10000, "peak_rss_mb": 100.0, "ok": True},
        {"engine": "streaming", "read_length": 1000, "n_rows": 1000000, "peak_rss_mb": 5000.0, "ok": True},
        {"engine": "duckdb", "read_length": 1000, "n_rows": 10000, "peak_rss_mb": 120.0, "ok": True},
        {"engine": "duckdb", "read_length": 1000, "n_rows": 1000000, "peak_rss_mb": 200.0, "ok": True},
    ]
    overall, _ = bench._gate_verdict(results)
    assert overall == "PASS"


def test_gate_verdict_indeterminate_with_no_successful_runs():
    results = [
        {"engine": "streaming", "read_length": 1000, "n_rows": 10000, "peak_rss_mb": 100.0, "ok": False},
    ]
    overall, detail = bench._gate_verdict(results)
    assert overall == "FAIL"
    assert "no successful runs" in detail or "insufficient" in detail


def test_format_table_handles_ok_and_failed_rows():
    results = [
        {"engine": "streaming", "n_rows": 1000, "read_length": 100, "peak_rss_mb": 50.0,
         "wall_seconds": 0.1, "n_groups": 1000, "input_mb": 0.3, "ok": True},
        {"engine": "eager", "n_rows": 1000, "read_length": 100, "input_mb": 0.3, "ok": False, "error": "boom"},
    ]
    table = bench._format_table(results)
    assert "streaming" in table
    assert "FAIL" in table
    assert "boom" in table


def test_peak_rss_mb_returns_positive_float():
    rss = bench._peak_rss_mb()
    assert isinstance(rss, float)
    assert rss > 0


def test_worker_rejects_unknown_engine(tmp_path):
    """The worker subprocess must reject an unknown engine name cleanly."""
    import argparse

    (tmp_path / "reads_00000.parquet").write_bytes(b"")  # placeholder so glob matches
    ns = argparse.Namespace(
        _worker=True, _workdir=str(tmp_path), _engine="bogus",
        _rows=1, _read_length=10, _unique_fraction=1.0,
    )
    rc = bench._worker(ns)
    assert rc == 1
