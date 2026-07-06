#!/usr/bin/env python3
"""Focused comparison: shell `sort | uniq -c` vs polars `group_by` for string dedup.

Answers S1b from the spike findings: does external merge-sort (which spills to
disk by design) keep peak RSS flat where hash `group_by` does not?

Two engines compared on the SAME parquet input:

  - polars_gb : scan_parquet → group_by(read_key).len() → sink_parquet
                (hash aggregate; Spike S1 showed RSS scales with unique count)
  - sort_uniq : polars streams read_key column to a text file, then
                `LC_ALL=C sort | uniq -c` (external merge-sort spills; uniq is
                O(1) streaming). Result written to text file.

Each config runs in a fresh subprocess so ru_maxrss is clean. For sort_uniq we
report BOTH RUSAGE_SELF (the python export process) and RUSAGE_CHILDREN (the
sort/uniq child processes — this is where the work happens).

Usage:
    python scripts/compare_sort_vs_groupby.py
    python scripts/compare_sort_vs_groupby.py --quick
    python scripts/compare_sort_vs_groupby.py --rows 10000,100000,1000000 --read-lengths 1000,5000
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

CHUNK_SIZE = 50_000


def _peak_rss_self_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return (rss / (1024 * 1024)) if sys.platform == "darwin" else rss / 1024


def _peak_rss_children_mb() -> float:
    """Max RSS of the largest child process (sort/uniq)."""
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return (rss / (1024 * 1024)) if sys.platform == "darwin" else rss / 1024


def _random_seqs_chunk(length: int, n: int, seed: int) -> list[str]:
    import numpy as np
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, 4, size=n * length, dtype=np.uint8)
    lookup = np.array([65, 67, 71, 84], dtype=np.uint8)  # ACGT
    b = lookup[idx].tobytes()
    return [b[i * length:(i + 1) * length].decode("ascii") for i in range(n)]


def generate_reads_parquet(out_dir: Path, n_rows: int, read_length: int,
                           chunk_size: int = CHUNK_SIZE) -> list[Path]:
    import polars as pl
    out_dir.mkdir(parents=True, exist_ok=True)
    shard_paths: list[Path] = []
    keys_written = 0
    shard_idx = 0
    while keys_written < n_rows:
        take = min(chunk_size, n_rows - keys_written)
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
    """Run one (rows, read_length, engine) config in this fresh process."""
    import glob as globlib
    shard_glob = str(Path(args._workdir) / "reads_*.parquet")
    glob_paths = sorted(globlib.glob(shard_glob))
    if not glob_paths:
        print(json.dumps({"error": f"no shards at {shard_glob}"}))
        return 1

    workdir = Path(args._workdir)
    out_path = workdir / "result.parquet"
    if out_path.exists():
        out_path.unlink()
    t0 = time.time()
    try:
        if args._engine == "polars_gb":
            import polars as pl
            pl.scan_parquet(shard_glob).group_by("read_key").len().sink_parquet(str(out_path))
            child_rss = 0.0
        elif args._engine == "sort_uniq":
            # 1. Stream read_key column to a text file in small batches via pyarrow.
            #    Batch size is ADAPTIVE: target ~50MB of string data per batch so
            #    the export process peak RSS is bounded by the byte budget, not by
            #    read_length. (50k rows x 5kb = 250MB would otherwise dominate RSS.)
            import pyarrow.parquet as pq
            BATCH_BYTES = 50_000_000  # 50 MB target per batch
            batch_size = max(1000, BATCH_BYTES // max(1, args._read_length))
            keys_file = str(workdir / "keys.txt")
            with open(keys_file, "w") as kf:
                for shard_path in glob_paths:
                    pf = pq.ParquetFile(shard_path)
                    for batch in pf.iter_batches(columns=["read_key"], batch_size=batch_size):
                        col = batch.column(0)
                        for s in col:
                            kf.write(s.as_py())
                            kf.write("\n")
                        del col, batch
                    gc.collect()
            self_rss_after_export = _peak_rss_self_mb()
            # 2. External merge-sort + streaming count. -S caps the in-memory sort
            #    buffer (forces spilling once exceeded); -T sets the spill dir.
            result_file = str(workdir / "counts.txt")
            tmpdir = str(workdir / "sort_tmp")
            os.makedirs(tmpdir, exist_ok=True)
            with open(result_file, "w") as fout:
                subprocess.run(
                    f"LC_ALL=C sort -S 64M -T {tmpdir!r} {keys_file!r} | uniq -c",
                    shell=True, stdout=fout, check=True,
                )
            child_rss = _peak_rss_children_mb()
        else:
            print(json.dumps({"error": f"unknown engine {args._engine}"}))
            return 1
        wall = time.time() - t0
        self_rss = _peak_rss_self_mb()

        # Correctness: count output groups. For unique_fraction=1.0, == input rows.
        if args._engine == "polars_gb":
            import polars as pl
            out_n = pl.scan_parquet(str(out_path)).select(pl.len()).collect().item()
        else:
            result_file = str(workdir / "counts.txt")
            out_n = int(subprocess.run(
                f"wc -l < {result_file!r}", shell=True, capture_output=True, text=True, check=True,
            ).stdout.strip() or "0")
        ok = (out_n == args._rows) if args._unique_fraction == 1.0 else True
    except Exception as e:
        wall = time.time() - t0
        self_rss = _peak_rss_self_mb()
        child_rss = 0.0
        out_n = None
        ok = False
        print(json.dumps({
            "engine": args._engine, "n_rows": args._rows, "read_length": args._read_length,
            "error": str(e)[:200], "ok": False,
        }))
        return 1

    # Peak RSS of the whole operation = max(self, children). For polars_gb there
    # are no children, so it's just self. For sort_uniq, the interesting number
    # is the children RSS (sort), but the operation peak is the max of both.
    peak_rss = max(self_rss, child_rss)
    print(json.dumps({
        "engine": args._engine, "n_rows": args._rows, "read_length": args._read_length,
        "peak_rss_mb": round(peak_rss, 1),
        "self_rss_mb": round(self_rss, 1),
        "child_rss_mb": round(child_rss, 1),
        "wall_seconds": round(wall, 2),
        "n_groups": out_n, "ok": ok,
    }))
    return 0


def run_config(rows, read_length, engine, workdir_root):
    workdir = workdir_root / f"r{rows}_l{read_length}_{engine}"
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    gen_t0 = time.time()
    generate_reads_parquet(workdir, rows, read_length)
    gen_seconds = round(time.time() - gen_t0, 2)
    input_mb = round(sum(os.path.getsize(p) for p in workdir.glob("reads_*.parquet")) / (1024 * 1024), 1)
    cmd = [sys.executable, __file__, "--_worker", "--_workdir", str(workdir),
           "--_engine", engine, "--_rows", str(rows),
           "--_read-length", str(read_length), "--_unique-fraction", "1.0"]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    lines = [l for l in proc.stdout.strip().splitlines() if l.strip()]
    result = json.loads(lines[-1]) if lines else {"error": proc.stderr.strip()}
    result["gen_seconds"] = gen_seconds
    result["input_mb"] = input_mb
    shutil.rmtree(workdir, ignore_errors=True)
    return result


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rows", default="10000,100000,1000000")
    p.add_argument("--read-lengths", default="1000,5000")
    p.add_argument("--engines", default="polars_gb,sort_uniq")
    p.add_argument("--quick", action="store_true")
    p.add_argument("--workdir", default=None)
    p.add_argument("--out", default=None, help="write markdown report to this path")
    p.add_argument("--_worker", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--_workdir", default=None, help=argparse.SUPPRESS)
    p.add_argument("--_engine", default=None, help=argparse.SUPPRESS)
    p.add_argument("--_rows", type=int, default=None, help=argparse.SUPPRESS)
    p.add_argument("--_read-length", type=int, default=None, help=argparse.SUPPRESS)
    p.add_argument("--_unique-fraction", type=float, default=1.0, help=argparse.SUPPRESS)
    args = p.parse_args(argv)
    if args._worker:
        return _worker(args)

    rows_list = [1000, 10000, 100000] if args.quick else [int(x) for x in args.rows.split(",")]
    lengths_list = [1000, 2000] if args.quick else [int(x) for x in args.read_lengths.split(",")]
    engines = [e.strip() for e in args.engines.split(",")]
    workdir_root = Path(args.workdir) if args.workdir else Path(tempfile.mkdtemp(prefix="cmp_sort_"))

    results = []
    print(f"# sort|uniq -c vs polars group_by  (workdir: {workdir_root})", flush=True)
    for rl in lengths_list:
        for rows in rows_list:
            for eng in engines:
                print(f"running {eng} rows={rows} len={rl}bp ... ", end="", flush=True)
                r = run_config(rows, rl, eng, workdir_root)
                results.append(r)
                print(f"{'ok' if r.get('ok') else 'FAIL'} peakRSS={r.get('peak_rss_mb', '?')}MB wall={r.get('wall_seconds', '?')}s", flush=True)

    print("\n## Results\n")
    print(f"{'engine':<12} {'rows':>9} {'len':>6} {'peakRSS(MB)':>11} {'self(MB)':>9} {'child(MB)':>10} {'wall(s)':>8} {'groups':>9} {'ok':>4}")
    print("-" * 90)
    for r in results:
        if r.get("ok"):
            print(f"{r['engine']:<12} {r['n_rows']:>9} {r['read_length']:>6} {r.get('peak_rss_mb', '?'):>11} "
                  f"{r.get('self_rss_mb', '?'):>9} {r.get('child_rss_mb', '?'):>10} "
                  f"{r.get('wall_seconds', '?'):>8} {r.get('n_groups', '?'):>9} {'Y':>4}")
        else:
            print(f"{r['engine']:<12} {r['n_rows']:>9} {r['read_length']:>6} {'FAIL':>11} - - {r.get('error', '')[:50]}")

    # RSS-flat analysis per engine
    print("\n## RSS scaling (max_rows / min_rows per read length)\n")
    for eng in engines:
        eng_results = [r for r in results if r.get("engine") == eng and r.get("ok") and "peak_rss_mb" in r]
        for rl in sorted({r["read_length"] for r in eng_results}):
            pts = sorted([r for r in eng_results if r["read_length"] == rl], key=lambda r: r["n_rows"])
            if len(pts) >= 2:
                ratio = pts[-1]["peak_rss_mb"] / pts[0]["peak_rss_mb"] if pts[0]["peak_rss_mb"] > 0 else float("inf")
                status = "FLAT" if ratio < 3.0 else "SCALES"
                print(f"  {eng:<12} len={rl}b: {pts[0]['n_rows']}->{pts[-1]['n_rows']} rows = "
                      f"{pts[0]['peak_rss_mb']:.0f}->{pts[-1]['peak_rss_mb']:.0f} MB (ratio {ratio:.2f}) [{status}]")

    if args.out:
        out_path = Path(args.out)
        out_path.write_text(_format_report(results, engines, rows_list, lengths_list))
        print(f"\nreport written to: {out_path}", flush=True)

    if not args.workdir:
        shutil.rmtree(workdir_root, ignore_errors=True)
    return 0


def _format_report(results, engines, rows_list, lengths_list) -> str:
    """Render a markdown report of the comparison results."""
    import time
    lines = [
        "# sort | uniq -c vs polars group_by — comparison",
        "",
        f"**Date:** {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"**Matrix:** rows={rows_list} x read_lengths={lengths_list}bp x engines={engines}",
        "**Unique fraction:** 1.0 (worst-case all-unique)",
        "**Gate:** peak RSS(max rows) / RSS(min rows) < 3.0 for a given read length",
        "",
        "## Results",
        "",
        "```",
        f"{'engine':<12} {'rows':>9} {'len':>6} {'peakRSS(MB)':>11} {'self(MB)':>9} {'child(MB)':>10} {'wall(s)':>8} {'groups':>9} {'ok':>4}",
        "-" * 90,
    ]
    for r in results:
        if r.get("ok"):
            lines.append(f"{r['engine']:<12} {r['n_rows']:>9} {r['read_length']:>6} {r.get('peak_rss_mb', '?'):>11} "
                         f"{r.get('self_rss_mb', '?'):>9} {r.get('child_rss_mb', '?'):>10} "
                         f"{r.get('wall_seconds', '?'):>8} {r.get('n_groups', '?'):>9} {'Y':>4}")
        else:
            lines.append(f"{r['engine']:<12} {r['n_rows']:>9} {r['read_length']:>6} FAIL {r.get('error', '')[:50]}")
    lines += ["", "## RSS scaling (max_rows / min_rows per read length)", ""]
    for eng in engines:
        eng_results = [r for r in results if r.get("engine") == eng and r.get("ok") and "peak_rss_mb" in r]
        for rl in sorted({r["read_length"] for r in eng_results}):
            pts = sorted([r for r in eng_results if r["read_length"] == rl], key=lambda r: r["n_rows"])
            if len(pts) >= 2:
                ratio = pts[-1]["peak_rss_mb"] / pts[0]["peak_rss_mb"] if pts[0]["peak_rss_mb"] > 0 else float("inf")
                status = "FLAT" if ratio < 3.0 else "SCALES"
                lines.append(f"  {eng:<12} len={rl}b: {pts[0]['n_rows']}->{pts[-1]['n_rows']} rows = "
                              f"{pts[0]['peak_rss_mb']:.0f}->{pts[-1]['peak_rss_mb']:.0f} MB (ratio {ratio:.2f}) [{status}]")
    lines += ["", "## Reproduce", "", "```bash",
              "pixi run -e default python scripts/compare_sort_vs_groupby.py",
              "pixi run -e default python scripts/compare_sort_vs_groupby.py --quick", "```", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    raise SystemExit(main())
