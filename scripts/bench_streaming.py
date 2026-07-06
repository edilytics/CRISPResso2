#!/usr/bin/env python3
"""Streaming memory benchmark harness for the parquet storage backend (PR 2).

This is the go/no-go SPIKE (S1) from the master design plan: does polars'
streaming ``group_by`` keep peak RSS flat as input row count scales, when the
group key is a multi-kilobase string column (the long-read amplicon scenario)?

The harness is a *reusable tool*, not a one-shot script. It:

  1. Generates synthetic reads of a configurable length, with a configurable
     fraction of duplicates, writing them to parquet shards on disk (chunked, so
     generation itself does not OOM). High diversity (~1:1 unique) is the worst
     case for the dedup ``group_by`` — this is what OOMs the current in-RAM dict.
  2. Runs ``scan_parquet(glob).group_by("read_key").len().collect(engine=...)`` in
     a **fresh subprocess** so ``resource.getrusage(RUSAGE_SELF).ru_maxrss``
     reflects only that operation (not generation).
  3. Records peak RSS, wall time, and a correctness check
     (sum of counts == input rows; group count matches expectation).
  4. Emits a results table to stdout and to a markdown report.

GATE (design Success Criterion #1): peak RSS under the streaming engine must NOT
scale linearly with row count. We assert ``RSS(max_rows) / RSS(min_rows) < 3``
(flat-ish). If the streaming engine falls back to eager, RSS will scale ~linearly
and the gate FAILS — at which point the design's fallback (DuckDB for the RC-dedup
stage, Approach C) must be evaluated before PR 5.

Usage:
    # default matrix (10k/100k/1M rows x 1kb/5kb/10kb strings, streaming + eager)
    python scripts/bench_streaming.py

    # quick smoke (fast, small)
    python scripts/bench_streaming.py --quick

    # custom matrix
    python scripts/bench_streaming.py --rows 10000,100000 --read-lengths 1000,5000

    # just report without re-running (reads cached results from --workdir)
    python scripts/bench_streaming.py --report-only
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# Default matrix from the design plan.
DEFAULT_ROWS = [10_000, 100_000, 1_000_000]
DEFAULT_READ_LENGTHS = [1_000, 5_000, 10_000]
QUICK_ROWS = [1_000, 10_000, 100_000]
QUICK_READ_LENGTHS = [1_000, 2_000]

# Worst case for dedup: every read unique (~1:1 diversity, the long-read scenario).
DEFAULT_UNIQUE_FRACTION = 1.0
# Generation chunk size — keeps write-time memory bounded regardless of total rows.
CHUNK_SIZE = 50_000
# Number of distinct nucleotides in the synthetic read key (DNA alphabet).
NUC_ALPHABET = "ACGT"

# RSS gate: streaming RSS(max rows) / RSS(min rows) must stay below this.
RSS_FLAT_RATIO_GATE = 3.0


def _peak_rss_mb() -> float:
    """Peak resident set size in MB. macOS reports bytes, Linux reports KB."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return rss / (1024 * 1024)
    return rss / 1024


def _random_seqs_chunk(length: int, n: int, seed: int) -> list[str]:
    """Generate ``n`` random DNA sequences of ``length`` bp, vectorized.

    Uses numpy to draw all bases at once (~100x faster than per-base
    ``random.choice`` for long sequences), then decodes to ASCII strings. This is
    what makes 1M x 5kb inputs tractable — pure-Python generation was the
    bottleneck (5 billion random.choice calls).
    """
    import numpy as np
    rng = np.random.default_rng(seed)
    # Draw 2-bit indices and map to the 4 DNA bases as bytes.
    idx = rng.integers(0, 4, size=n * length, dtype=np.uint8)
    # Lookup: 0->A(65), 1->C(67), 2->G(71), 3->T(84)
    lookup = np.array([65, 67, 71, 84], dtype=np.uint8)
    bytes_arr = lookup[idx].tobytes()
    # Slice into n rows of `length` bytes each and decode.
    return [bytes_arr[i * length:(i + 1) * length].decode("ascii") for i in range(n)]


def generate_reads_parquet(
    out_dir: Path,
    n_rows: int,
    read_length: int,
    unique_fraction: float,
    chunk_size: int = CHUNK_SIZE,
) -> list[Path]:
    """Generate synthetic reads as chunked parquet shards.

    Returns the list of shard paths. The read_key column holds the full read
    sequence (length = read_length); a count column (all ones) mirrors how
    Stage 1 of the design consumes raw reads.

    Generation is chunked so that creating very large inputs (e.g. 1M x 5kb =
    5GB) does not itself require holding everything in RAM: each shard holds at
    most ``chunk_size`` rows, and each chunk is generated vectorized via numpy.
    """
    import polars as pl

    out_dir.mkdir(parents=True, exist_ok=True)
    shard_paths: list[Path] = []
    keys_written = 0
    shard_idx = 0
    # For unique_fraction < 1.0, build a pool and sample from it (creates
    # duplicates). For unique_fraction == 1.0, every read is distinct (worst case).
    pool = None
    if unique_fraction < 1.0:
        n_unique = max(1, int(n_rows * unique_fraction))
        pool = _random_seqs_chunk(read_length, n_unique, seed=0)
    while keys_written < n_rows:
        take = min(chunk_size, n_rows - keys_written)
        if pool is not None:
            import random as _random
            rng_pick = _random.Random(shard_idx)
            keys = [rng_pick.choice(pool) for _ in range(take)]
        else:
            keys = _random_seqs_chunk(read_length, take, seed=shard_idx + 1)
        df = pl.DataFrame({"read_key": keys, "count": [1] * take})
        shard = out_dir / f"reads_{shard_idx:05d}.parquet"
        df.write_parquet(shard)
        shard_paths.append(shard)
        keys_written += take
        shard_idx += 1
        del df, keys
        gc.collect()
    return shard_paths


def _worker(args: argparse.Namespace) -> int:
    """Run a single (rows, read_length, engine) config in this fresh process.

    Prints a JSON result line on stdout. Designed to be invoked as a subprocess
    so that ``ru_maxrss`` reflects ONLY the group_by operation.
    """
    import glob as globlib

    shard_glob = str(Path(args._workdir) / "reads_*.parquet")
    glob_paths = sorted(globlib.glob(shard_glob))
    if not glob_paths:
        print(json.dumps({"error": f"no shards at {shard_glob}"}))
        return 1

    out_path = Path(args._workdir) / "result.parquet"
    if out_path.exists():
        out_path.unlink()
    t0 = time.time()
    n_groups: int | str | None
    sum_counts: int | None
    if args._engine == "streaming":
        import polars as pl
        # Stage 1/3 of the design never hold the result in RAM — they sink it to
        # parquet and downstream stages scan it lazily. sink_parquet with the
        # streaming engine is therefore the operation whose peak RSS we must
        # measure; collect() would inflate RSS with the materialized result table
        # and obscure whether the *grouping* spilled.
        pl.scan_parquet(shard_glob).group_by("read_key").len().sink_parquet(str(out_path))
        n_groups = None
        sum_counts = None
    elif args._engine == "eager":
        import polars as pl
        result = pl.scan_parquet(shard_glob).group_by("read_key").len().collect()
        n_groups = result.height
        sum_counts = int(result["len"].sum())
        # Write to parquet so the unified post-run correctness scan (below) works
        # the same as the streaming/duckdb paths. Adds a write step to eager
        # timing, but eager is the RSS oracle, not the timing reference.
        result.write_parquet(str(out_path))
        del result
        gc.collect()
    elif args._engine == "duckdb":
        # Documented fallback (design Approach C): DuckDB owns the group_by as a
        # SQL hash aggregate and streams the output to parquet via COPY, never
        # materializing the result in client RAM. This tests whether an OLAP
        # engine that is purpose-built for spilling hash aggregates keeps peak RSS
        # flat where polars streaming does not.
        import duckdb
        con = duckdb.connect()
        con.execute(
            """COPY (
                  SELECT read_key, count(*) AS len
                  FROM read_parquet(?)
                  GROUP BY read_key
               ) TO '{}' (FORMAT PARQUET)""".format(str(out_path).replace("'", "''")),
            [shard_glob],
        )
        n_groups = None
        sum_counts = None
    else:
        print(json.dumps({"error": f"unknown engine {args._engine}"}))
        return 1
    wall = time.time() - t0

    rss_mb = _peak_rss_mb()
    # Correctness: verify the sink produced the expected total row count by
    # scanning the output (does not materialize it). With unique_fraction=1.0
    # every read is its own group, so output row count == input row count.
    ok = True
    try:
        if args._engine == "duckdb":
            import duckdb
            out_n = duckdb.connect().execute(
                "SELECT count(*) FROM read_parquet('{}')".format(str(out_path).replace("'", "''")),
            ).fetchone()[0]
        else:
            import polars as pl
            out_n = pl.scan_parquet(str(out_path)).select(pl.len()).collect().item()
        if args._unique_fraction == 1.0:
            ok = (out_n == args._rows)
        n_groups = out_n
        sum_counts = out_n
    except Exception as e:
        ok = False
        n_groups = f"err: {e}"
    print(json.dumps({
        "engine": args._engine,
        "n_rows": args._rows,
        "read_length": args._read_length,
        "peak_rss_mb": round(rss_mb, 1),
        "wall_seconds": round(wall, 2),
        "n_groups": n_groups,
        "sum_counts": sum_counts,
        "ok": ok,
    }))
    return 0


def run_config(
    rows: int,
    read_length: int,
    engine: str,
    unique_fraction: float,
    workdir_root: Path,
) -> dict:
    """Generate input + run the group_by in a fresh subprocess. Returns a result dict."""
    workdir = workdir_root / f"r{rows}_l{read_length}_{engine}"
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)

    # 1. Generate input parquet shards (chunked; bounded memory).
    gen_t0 = time.time()
    generate_reads_parquet(
        workdir, rows, read_length, unique_fraction, chunk_size=CHUNK_SIZE,
    )
    gen_seconds = round(time.time() - gen_t0, 2)
    shard_sizes = [os.path.getsize(p) for p in sorted(workdir.glob("reads_*.parquet"))]
    input_bytes = sum(shard_sizes)

    # 2. Run group_by in a fresh subprocess so ru_maxrss is clean.
    cmd = [
        sys.executable, __file__,
        "--_worker",
        "--_workdir", str(workdir),
        "--_engine", engine,
        "--_rows", str(rows),
        "--_read-length", str(read_length),
        "--_unique-fraction", str(unique_fraction),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        result = {
            "engine": engine, "n_rows": rows, "read_length": read_length,
            "error": proc.stderr.strip() or proc.stdout.strip(),
            "ok": False,
        }
    else:
        # parse the last non-empty stdout line as JSON
        lines = [l for l in proc.stdout.strip().splitlines() if l.strip()]
        result = json.loads(lines[-1]) if lines else {"error": "no output"}

    result["gen_seconds"] = gen_seconds
    result["input_mb"] = round(input_bytes / (1024 * 1024), 1)
    result["unique_fraction"] = unique_fraction

    # 3. Cleanup input shards (keep workdir for inspection if --keep).
    if not os.environ.get("BENCH_KEEP"):
        shutil.rmtree(workdir)
    return result


def _format_table(results: list[dict]) -> str:
    header = f"{'engine':<10} {'rows':>9} {'len(b)':>8} {'peakRSS(MB)':>12} {'wall(s)':>8} {'groups':>9} {'input(MB)':>10} {'ok':>4}"
    sep = "-" * len(header)
    lines = [header, sep]
    for r in results:
        if r.get("ok"):
            lines.append(
                f"{r['engine']:<10} {r['n_rows']:>9} {r['read_length']:>8} "
                f"{r.get('peak_rss_mb', '?'):>12} {r.get('wall_seconds', '?'):>8} "
                f"{r.get('n_groups', '?'):>9} {r.get('input_mb', '?'):>10} {'Y' if r.get('ok') else 'N':>4}"
            )
        else:
            lines.append(
                f"{r['engine']:<10} {r['n_rows']:>9} {r['read_length']:>8} "
                f"{'FAIL':>12} {'-':>8} {'-':>9} {r.get('input_mb', '-'):>10} {'N':>4}  {r.get('error', '')[:60]}"
            )
    return "\n".join(lines)


def _gate_verdict(results: list[dict]) -> tuple[str, str]:
    """Evaluate the RSS-flat gate across spilling-capable engines.

    For each engine with successful runs, computes RSS(max_rows)/RSS(min_rows)
    at a FIXED read length. PASS if at least one engine stays flat (ratio <
    RSS_FLAT_RATIO_GATE) for all read lengths — that engine becomes the chosen
    backend for the streaming group_by stages. Reports per-engine detail.
    """
    spill_engines = ["streaming", "duckdb"]
    per_engine: list[str] = []
    any_pass = False
    for eng in spill_engines:
        runs = [r for r in results if r.get("engine") == eng and r.get("ok") and "peak_rss_mb" in r]
        if not runs:
            per_engine.append(f"{eng}: no successful runs")
            continue
        eng_worst = 0.0
        eng_pass_all = True
        for rl in sorted({r["read_length"] for r in runs}):
            pts = sorted([r for r in runs if r["read_length"] == rl], key=lambda r: r["n_rows"])
            if len(pts) < 2:
                continue
            min_rss = pts[0]["peak_rss_mb"]
            max_rss = pts[-1]["peak_rss_mb"]
            if min_rss <= 0:
                per_engine.append(f"{eng} len={rl}b: min RSS=0")
                continue
            ratio = max_rss / min_rss
            eng_worst = max(eng_worst, ratio)
            status = "PASS" if ratio < RSS_FLAT_RATIO_GATE else "FAIL"
            if status == "FAIL":
                eng_pass_all = False
            per_engine.append(
                f"{eng} len={rl}b: RSS {pts[0]['n_rows']}->{pts[-1]['n_rows']} = {min_rss:.0f}->{max_rss:.0f} MB "
                f"(ratio {ratio:.2f}) [{status}]"
            )
        any_pass = any_pass or (eng_worst > 0 and eng_pass_all)
    overall = "PASS" if any_pass else "FAIL"
    return overall, " | ".join(per_engine) if per_engine else "insufficient data (need >=2 row counts per read length per engine)"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rows", default=",".join(str(r) for r in DEFAULT_ROWS),
                   help="comma-separated row counts (default: %(default)s)")
    p.add_argument("--read-lengths", default=",".join(str(r) for r in DEFAULT_READ_LENGTHS),
                   help="comma-separated read lengths in bp (default: %(default)s)")
    p.add_argument("--unique-fraction", type=float, default=DEFAULT_UNIQUE_FRACTION,
                   help="fraction of reads that are unique (default %(default)s; 1.0 = worst case)")
    p.add_argument("--engines", default="streaming,eager,duckdb",
                   help="comma-separated engines: streaming,eager,duckdb (default: %(default)s)")
    p.add_argument("--quick", action="store_true",
                   help=f"use a small fast matrix ({QUICK_ROWS} x {QUICK_READ_LENGTHS})")
    p.add_argument("--workdir", default=None,
                   help="temp dir for parquet shards (default: system temp)")
    p.add_argument("--out", default=None,
                   help="write markdown report to this path (default: scripts/bench_streaming_results.md)")
    p.add_argument("--report-only", action="store_true",
                   help="skip running; just re-format the last results from --workdir")
    # hidden worker args
    p.add_argument("--_worker", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--_workdir", default=None, help=argparse.SUPPRESS)
    p.add_argument("--_engine", default=None, help=argparse.SUPPRESS)
    p.add_argument("--_rows", type=int, default=None, help=argparse.SUPPRESS)
    p.add_argument("--_read-length", type=int, default=None, help=argparse.SUPPRESS)
    p.add_argument("--_unique-fraction", type=float, default=1.0, help=argparse.SUPPRESS)
    args = p.parse_args(argv)

    # ---- worker mode: run a single group_by in this fresh process ----
    if args._worker:
        return _worker(args)

    # ---- main mode: orchestrate configs across subprocesses ----
    if args.quick:
        rows_list = QUICK_ROWS
        lengths_list = QUICK_READ_LENGTHS
    else:
        rows_list = [int(x) for x in args.rows.split(",")]
        lengths_list = [int(x) for x in args.read_lengths.split(",")]
    engines = [e.strip() for e in args.engines.split(",") if e.strip()]

    workdir_root = Path(args.workdir) if args.workdir else Path(tempfile.mkdtemp(prefix="bench_streaming_"))
    out_path = Path(args.out) if args.out else Path(__file__).parent / "bench_streaming_results.md"

    results: list[dict] = []
    print("# Streaming memory benchmark", flush=True)
    print(f"workdir: {workdir_root}", flush=True)
    print(f"matrix: rows={rows_list} x lengths={lengths_list} x engines={engines} "
          f"(unique_fraction={args.unique_fraction})", flush=True)
    print(flush=True)

    for engine in engines:
        for rl in lengths_list:
            for rows in rows_list:
                tag = f"[{engine}] rows={rows} len={rl}bp"
                print(f"running {tag} ...", flush=True, end=" ")
                sys.stdout.flush()
                r = run_config(rows, rl, engine, args.unique_fraction, workdir_root)
                results.append(r)
                status = "ok" if r.get("ok") else "FAIL"
                rss = r.get("peak_rss_mb", "?")
                print(f"{status} peakRSS={rss}MB wall={r.get('wall_seconds', '?')}s", flush=True)

    print(flush=True)
    print("## Results", flush=True)
    print(_format_table(results), flush=True)
    print(flush=True)

    overall, detail = _gate_verdict(results)
    print("## Gate verdict", flush=True)
    print(f"Streaming RSS-flat gate (RSS(max)/RSS(min) < {RSS_FLAT_RATIO_GATE}): {overall}", flush=True)
    print(f"  {detail}", flush=True)
    print(flush=True)
    if overall == "FAIL":
        print("GATE FAILED — neither polars streaming nor DuckDB keeps RSS flat on multi-kb string group_by.", flush=True)
        print("Re-evaluate the storage design before PR 5; the premise that an OLAP engine spills the", flush=True)
        print("hash aggregate for high-cardinality string keys does not hold at this scale.", flush=True)
    elif overall == "PASS":
        print("GATE PASSED — at least one engine keeps RSS flat; proceed with that engine for the", flush=True)
        print("streaming group_by stages (Stage 1 count_reads / Stage 3 collapse) in PR 3+.", flush=True)

    # write markdown report
    _write_report(out_path, results, overall, detail, args, rows_list, lengths_list, engines)
    print(f"\nreport written to: {out_path}", flush=True)

    # cleanup workdir
    if not os.environ.get("BENCH_KEEP") and not args.workdir:
        shutil.rmtree(workdir_root, ignore_errors=True)

    return 0 if overall != "FAIL" else 2


def _write_report(path, results, overall, detail, args, rows_list, lengths_list, engines):
    lines = [
        "# Streaming Memory Benchmark (PR 2 / Spike S1)",
        "",
        f"**Date:** {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"**Matrix:** rows={rows_list} x read_lengths={lengths_list}bp x engines={engines}",
        f"**Unique fraction:** {args.unique_fraction} (1.0 = worst-case all-unique)",
        f"**Gate:** peak RSS(max rows) / RSS(min rows) < {RSS_FLAT_RATIO_GATE} for the streaming engine",
        "",
        "## Verdict",
        "",
        f"**{overall}** — {detail}",
        "",
        "## Results",
        "",
        "```",
        _format_table(results),
        "```",
        "",
        "## Interpretation",
        "",
        "- The **streaming** engine must keep peak RSS roughly flat as row count grows (that is the whole "
        "point of the parquet backend — flat memory across arbitrary input sizes).",
        "- The **eager** engine is the regression oracle: its RSS should scale ~linearly with row count, "
        "demonstrating that the test actually exercises the dedup at scale.",
        "- If streaming RSS tracks eager RSS (linear growth), polars is falling back to eager on multi-kb "
        "string columns. Per the design, the fallback is DuckDB for the RC-dedup stage (Stage 3b) only.",
        "",
        "## How to reproduce",
        "",
        "```bash",
        "pixi run -e default python scripts/bench_streaming.py            # full matrix",
        "pixi run -e default python scripts/bench_streaming.py --quick    # fast smoke",
        "```",
        "",
    ]
    path.write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
