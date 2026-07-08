#!/usr/bin/env python3
"""Spike: does polars streaming ``sort`` on a multi-kb string key spill to disk?

Answers the one open implementation question in
``design_docs/STREAMING_SINGLE_READ_COLLAPSE.md``: for the streaming
single-read collapse we need to external-sort by ``canonical_key`` (a string
~amplicon-length, 1-10 kb). The proven fallback is system ``sort`` (external
merge-sort, ``-S`` buffer cap — spike S1b confirmed flat RSS). The preferred
mechanism is polars ``scan_parquet().sort("canonical_key").sink_parquet()`` —
*if* it spills (peak RSS flat in row count). This spike decides.

Two engines on the SAME generated parquet input (all-unique 2 kb
``canonical_key`` strings, worst case — no collapse, every row distinct):

  - polars_sort : scan_parquet → sort("canonical_key") → sink_parquet
                  (polars streaming sort; spills iff RSS stays flat)
  - system_sort : pyarrow streams canonical_key+seq_no to a text file, then
                  `LC_ALL=C sort -S 64M` (external merge-sort, proven flat)

Each (rows, engine) config runs in a FRESH subprocess so ``ru_maxrss`` reflects
only the sort (generation happens in a prior step, discarded from the
measurement). Input is parquet on disk, scanned lazily, so spilling is
testable.

Gate: ``RSS(max_rows) / RSS(min_rows) < 3.0`` for a spilling engine, fixed
key length. If polars_sort passes → use it (cleanest, one polars call, native
parquet round-trip). If it fails → system_sort is the proven fallback.

Usage:
    pixi run -e default python scripts/bench_polars_sort_spill.py
    pixi run -e default python scripts/bench_polars_sort_spill.py --quick
    pixi run -e default python scripts/bench_polars_sort_spill.py \\
        --rows 10000,100000,1000000 --key-length 2000
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
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return (rss / (1024 * 1024)) if sys.platform == "darwin" else rss / 1024


def _random_keys_chunk(length: int, n: int, seed: int) -> list[str]:
    import numpy as np
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, 4, size=n * length, dtype=np.uint8)
    lookup = np.array([65, 67, 71, 84], dtype=np.uint8)  # ACGT
    b = lookup[idx].tobytes()
    return [b[i * length:(i + 1) * length].decode("ascii") for i in range(n)]


def generate_keys_parquet(out_dir: Path, n_rows: int, key_length: int,
                          chunk_size: int = CHUNK_SIZE) -> list[Path]:
    """Write n_rows all-unique ACGT strings (key_length each) to parquet shards.

    Schema mirrors the streaming-collapse spill parquet: ``canonical_key``
    (string), ``seq_no`` (int64). One row per read.
    """
    import polars as pl
    out_dir.mkdir(parents=True, exist_ok=True)
    shard_paths: list[Path] = []
    keys_written = 0
    shard_idx = 0
    while keys_written < n_rows:
        take = min(chunk_size, n_rows - keys_written)
        keys = _random_keys_chunk(key_length, take, seed=shard_idx + 1)
        df = pl.DataFrame({
            "canonical_key": keys,
            "seq_no": list(range(keys_written, keys_written + take)),
        })
        shard = out_dir / f"keys_{shard_idx:05d}.parquet"
        df.write_parquet(shard)
        shard_paths.append(shard)
        keys_written += take
        shard_idx += 1
        del df, keys
        gc.collect()
    return shard_paths


def _worker(args: argparse.Namespace) -> int:
    """Run one (rows, key_length, engine) config in this fresh process."""
    import glob as globlib
    shard_glob = str(Path(args._workdir) / "keys_*.parquet")
    glob_paths = sorted(globlib.glob(shard_glob))
    if not glob_paths:
        print(json.dumps({"error": f"no shards at {shard_glob}"}))
        return 1

    workdir = Path(args._workdir)
    out_path = workdir / "sorted.parquet"
    if out_path.exists():
        out_path.unlink()
    t0 = time.time()
    out_n = None
    ok = True
    err = ""
    child_rss = 0.0
    try:
        if args._engine == "polars_sort":
            import polars as pl
            (pl.scan_parquet(shard_glob)
               .sort("canonical_key")
               .sink_parquet(str(out_path)))
            child_rss = 0.0
            out_n = pl.scan_parquet(str(out_path)).select(pl.len()).collect().item()
        elif args._engine == "system_sort":
            # 1. Stream canonical_key + seq_no to a text file in adaptive batches
            #    (target ~50MB of string data per batch so the export process peak
            #    RSS is bounded by the byte budget, not by key_length).
            import pyarrow.parquet as pq
            BATCH_BYTES = 50_000_000  # 50 MB target per batch
            batch_size = max(1000, BATCH_BYTES // max(1, args._key_length))
            keys_file = str(workdir / "keys.txt")
            with open(keys_file, "w") as kf:
                for shard_path in glob_paths:
                    pf = pq.ParquetFile(shard_path)
                    for batch in pf.iter_batches(columns=["canonical_key", "seq_no"],
                                                 batch_size=batch_size):
                        keys = batch.column(0)
                        seqs = batch.column(1)
                        for i in range(batch.num_rows):
                            kf.write(keys[i].as_py())
                            kf.write("\t")
                            kf.write(str(seqs[i].as_py()))
                            kf.write("\n")
                        del keys, seqs, batch
                    gc.collect()
            # 2. External merge-sort by canonical_key (field 1, TAB-separated).
            #    -S caps the in-memory buffer (forces disk spill); -T sets spill
            #    dir; -s stable so equal keys keep file order (parity: min seq_no
            #    wins within a group). LC_ALL=C forces byte-wise comparison.
            sorted_file = str(workdir / "keys.sorted")
            tmpdir = str(workdir / "sort_tmp")
            os.makedirs(tmpdir, exist_ok=True)
            subprocess.run(
                ["sort", "-S", "64M", "-T", tmpdir, "-s",
                 "-t", "\t", "-k1,1",
                 "-o", sorted_file, keys_file],
                check=True, env=dict(os.environ, LC_ALL="C"),
            )
            child_rss = _peak_rss_children_mb()
            out_n = int(subprocess.run(
                f"wc -l < {sorted_file!r}", shell=True, capture_output=True,
                text=True, check=True,
            ).stdout.strip() or "0")
        else:
            print(json.dumps({"error": f"unknown engine {args._engine}"}))
            return 1
    except Exception as e:
        err = str(e)[:300]
        ok = False
    wall = time.time() - t0
    self_rss = _peak_rss_self_mb()
    peak_rss = max(self_rss, child_rss)
    print(json.dumps({
        "engine": args._engine, "n_rows": args._rows, "key_length": args._key_length,
        "peak_rss_mb": round(peak_rss, 1),
        "self_rss_mb": round(self_rss, 1),
        "child_rss_mb": round(child_rss, 1),
        "wall_seconds": round(wall, 2),
        "n_out": out_n, "ok": ok, "error": err,
    }))
    return 0


def run_config(rows, key_length, engine, workdir_root):
    workdir = workdir_root / f"r{rows}_l{key_length}_{engine}"
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    gen_t0 = time.time()
    generate_keys_parquet(workdir, rows, key_length)
    gen_seconds = round(time.time() - gen_t0, 2)
    input_mb = round(sum(os.path.getsize(p) for p in workdir.glob("keys_*.parquet")) / (1024 * 1024), 1)
    cmd = [sys.executable, __file__, "--_worker", "--_workdir", str(workdir),
           "--_engine", engine, "--_rows", str(rows),
           "--_key-length", str(key_length)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    lines = [l for l in proc.stdout.strip().splitlines() if l.strip()]
    result = json.loads(lines[-1]) if lines else {"error": proc.stderr.strip(), "ok": False}
    result["gen_seconds"] = gen_seconds
    result["input_mb"] = input_mb
    shutil.rmtree(workdir, ignore_errors=True)
    return result


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rows", default="10000,100000,1000000")
    p.add_argument("--key-lengths", default="2000")
    p.add_argument("--engines", default="polars_sort,system_sort")
    p.add_argument("--quick", action="store_true")
    p.add_argument("--workdir", default=None)
    p.add_argument("--out", default=None, help="write markdown report to this path")
    p.add_argument("--_worker", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--_workdir", default=None, help=argparse.SUPPRESS)
    p.add_argument("--_engine", default=None, help=argparse.SUPPRESS)
    p.add_argument("--_rows", type=int, default=None, help=argparse.SUPPRESS)
    p.add_argument("--_key-length", type=int, default=None, help=argparse.SUPPRESS)
    args = p.parse_args(argv)
    if args._worker:
        return _worker(args)

    rows_list = [1000, 10000, 100000] if args.quick else [int(x) for x in args.rows.split(",")]
    lengths_list = [1000, 2000] if args.quick else [int(x) for x in args.key_lengths.split(",")]
    engines = [e.strip() for e in args.engines.split(",")]
    workdir_root = Path(args.workdir) if args.workdir else Path(tempfile.mkdtemp(prefix="spike_sort_"))

    results = []
    print(f"# polars sort spill spike  (workdir: {workdir_root})", flush=True)
    for kl in lengths_list:
        for rows in rows_list:
            for eng in engines:
                print(f"running {eng} rows={rows} len={kl}bp ... ", end="", flush=True)
                r = run_config(rows, kl, eng, workdir_root)
                results.append(r)
                print(f"{'ok' if r.get('ok') else 'FAIL'} peakRSS={r.get('peak_rss_mb', '?')}MB wall={r.get('wall_seconds', '?')}s", flush=True)

    print("\n## Results\n")
    print(f"{'engine':<14} {'rows':>9} {'len':>6} {'peakRSS(MB)':>11} {'self(MB)':>9} {'child(MB)':>10} {'wall(s)':>8} {'n_out':>9} {'ok':>4}")
    print("-" * 95)
    for r in results:
        if r.get("ok"):
            print(f"{r['engine']:<14} {r['n_rows']:>9} {r['key_length']:>6} {r.get('peak_rss_mb', '?'):>11} "
                  f"{r.get('self_rss_mb', '?'):>9} {r.get('child_rss_mb', '?'):>10} "
                  f"{r.get('wall_seconds', '?'):>8} {r.get('n_out', '?'):>9} {'Y':>4}")
        else:
            print(f"{r['engine']:<14} {r['n_rows']:>9} {r['key_length']:>6} {'FAIL':>11} - - {r.get('error', '')[:50]}")

    print("\n## RSS scaling (max_rows / min_rows per key length)\n")
    for eng in engines:
        eng_results = [r for r in results if r.get("engine") == eng and r.get("ok") and "peak_rss_mb" in r]
        for kl in sorted({r["key_length"] for r in eng_results}):
            pts = sorted([r for r in eng_results if r["key_length"] == kl], key=lambda r: r["n_rows"])
            if len(pts) >= 2:
                ratio = pts[-1]["peak_rss_mb"] / pts[0]["peak_rss_mb"] if pts[0]["peak_rss_mb"] > 0 else float("inf")
                status = "FLAT" if ratio < 3.0 else "SCALES"
                print(f"  {eng:<14} len={kl}b: {pts[0]['n_rows']}->{pts[-1]['n_rows']} rows = "
                      f"{pts[0]['peak_rss_mb']:.0f}->{pts[-1]['peak_rss_mb']:.0f} MB (ratio {ratio:.2f}) [{status}]")

    if args.out:
        out_path = Path(args.out)
        out_path.write_text(_format_report(results, engines, rows_list, lengths_list))
        print(f"\nreport written to: {out_path}", flush=True)

    if not args.workdir:
        shutil.rmtree(workdir_root, ignore_errors=True)
    return 0


def _format_report(results, engines, rows_list, lengths_list) -> str:
    import time
    lines = [
        "# polars sort spill spike",
        "",
        f"**Date:** {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"**Matrix:** rows={rows_list} x key_lengths={lengths_list}bp x engines={engines}",
        "**Key cardinality:** all-unique (worst case — no sort collapse)",
        "**Gate:** peak RSS(max rows) / RSS(min rows) < 3.0 for a given key length",
        "",
        "Context: `design_docs/STREAMING_SINGLE_READ_COLLAPSE.md` step 3 — the one",
        "open implementation question. polars `sort`+`sink_parquet` is preferred",
        "(cleanest); system `sort` (external merge-sort, `-S 64M`) is the proven",
        "fallback (spike S1b).",
        "",
        "## Results",
        "",
        "```",
        f"{'engine':<14} {'rows':>9} {'len':>6} {'peakRSS(MB)':>11} {'self(MB)':>9} {'child(MB)':>10} {'wall(s)':>8} {'n_out':>9} {'ok':>4}",
        "-" * 95,
    ]
    for r in results:
        if r.get("ok"):
            lines.append(f"{r['engine']:<14} {r['n_rows']:>9} {r['key_length']:>6} {r.get('peak_rss_mb', '?'):>11} "
                         f"{r.get('self_rss_mb', '?'):>9} {r.get('child_rss_mb', '?'):>10} "
                         f"{r.get('wall_seconds', '?'):>8} {r.get('n_out', '?'):>9} {'Y':>4}")
        else:
            lines.append(f"{r['engine']:<14} {r['n_rows']:>9} {r['key_length']:>6} FAIL {r.get('error', '')[:50]}")
    lines += ["", "## RSS scaling (max_rows / min_rows per key length)", ""]
    for eng in engines:
        eng_results = [r for r in results if r.get("engine") == eng and r.get("ok") and "peak_rss_mb" in r]
        for kl in sorted({r["key_length"] for r in eng_results}):
            pts = sorted([r for r in eng_results if r["key_length"] == kl], key=lambda r: r["n_rows"])
            if len(pts) >= 2:
                ratio = pts[-1]["peak_rss_mb"] / pts[0]["peak_rss_mb"] if pts[0]["peak_rss_mb"] > 0 else float("inf")
                status = "FLAT" if ratio < 3.0 else "SCALES"
                lines.append(f"  {eng:<14} len={kl}b: {pts[0]['n_rows']}->{pts[-1]['n_rows']} rows = "
                              f"{pts[0]['peak_rss_mb']:.0f}->{pts[-1]['peak_rss_mb']:.0f} MB (ratio {ratio:.2f}) [{status}]")
    lines += ["", "## Reproduce", "", "```bash",
              "pixi run -e default python scripts/bench_polars_sort_spill.py",
              "pixi run -e default python scripts/bench_polars_sort_spill.py --quick", "```", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    raise SystemExit(main())
