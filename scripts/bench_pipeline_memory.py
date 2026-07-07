#!/usr/bin/env python3
"""End-to-end peak-RSS benchmark for the CRISPResso2 storage backends.

Measures Success Criterion #1 from the flat-memory design: *peak RSS does not
increase as read count scales*. Generates synthetic high-diversity FASTQs on a
fixed amplicon, runs ``CRISPResso`` end-to-end under each backend via
``/usr/bin/time -l``, and reports the peak resident set size.

Why high diversity: the memory ceiling is the per-unique-read payload held in
RAM. With ~1:1 read diversity (every read unique — the long-read regime) the
pandas ``variantCache`` grows linearly in the read count, while the parquet
backend streams and should stay flat. (Heavy dedup would hide the ceiling:
few unique reads → small cache → no scaling.)

Usage:
    python scripts/bench_pipeline_memory.py --reads 10000,100000,1000000 \\
        --backends pandas,parquet --workdir /tmp/crispresso_bench

Outputs a table and a markdown summary at ``<workdir>/results.md``.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

# A ~2 kb synthetic amplicon with a 20 bp guide embedded near the center.
# Fixed seed so runs are reproducible.
_GUID = "GGAATCCCTTCTGCAGCACC"


def _make_amplicon(length: int = 2000, seed: int = 42) -> str:
    """Random 2kb amplicon (ACGT) with the guide embedded at position length//2."""
    import random

    rng = random.Random(seed)
    bases = "ACGT"
    seq = "".join(rng.choice(bases) for _ in range(length))
    mid = length // 2
    return seq[:mid] + _GUID + seq[mid + len(_GUID):]


def _generate_fastq(path: Path, n_reads: int, amplicon: str, seed: int = 7) -> int:
    """Write a FASTQ of n_reads high-diversity reads.

    Each read = the amplicon with 1-3 random single-base substitutions, so
    ~every read is unique (worst case for dedup / the long-read regime).
    Quality is all 'I' (phred 40). Returns the file size in bytes.
    """
    import random

    rng = random.Random(seed)
    bases = "ACGT"
    alen = len(amplicon)
    amp_bytes = list(amplicon)
    with open(path, "w") as fh:
        for i in range(n_reads):
            n_subs = rng.randint(1, 3)
            arr = amp_bytes[:]
            for _ in range(n_subs):
                pos = rng.randrange(alen)
                arr[pos] = rng.choice(bases)
            seq = "".join(arr)
            fh.write(f"@read{i}\n{seq}\n+\n{'I' * alen}\n")
    return os.path.getsize(path)


def _peak_rss_bytes(child_stdout: str) -> int | None:
    """Parse 'maximum resident set size' (bytes) from /usr/bin/time -l stderr.

    macOS prints it as a line ``  <bytes>  maximum resident set size``. The
    time wrapper interleaves it with the child's own stderr; we scan all lines.
    """
    for line in child_stdout.splitlines():
        if "maximum resident set size" in line:
            try:
                return int(line.split("maximum resident set size")[0].strip())
            except ValueError:
                pass
    return None


def _run_one(backend: str, fastq: Path, amplicon: str, out_dir: Path,
             n_proc: int = 2, timeout_s: int = 1800) -> dict:
    """Run CRISPResso under one backend via /usr/bin/time -l; return RSS + timing.

    Uses the installed ``CRISPResso`` console script (not ``python -m``) because
    ``CRISPRessoCORE.py`` has no ``__main__`` guard — ``main()`` is only
    invoked by the console-script entry point.
    """
    crispresso = str(Path(sys.executable).parent / "CRISPResso")
    cmd = [
        "/usr/bin/time", "-l",
        crispresso,
        "-r1", str(fastq),
        "-a", amplicon,
        "-g", _GUID,
        "--n_processes", str(n_proc),
        "--storage_backend", backend,
        "-n", f"bench_{backend}",
        "--place_report_in_output_folder",
        "--halt_on_plot_fail",
        "--suppress_plots",          # skip plot rendering — we only measure the pipeline
        "--suppress_report",
        "--output_folder", str(out_dir),
    ]
    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout_s, check=False,
        )
    except subprocess.TimeoutExpired:
        return {"backend": backend, "rss_mb": None, "wall_s": timeout_s,
                "status": "TIMEOUT"}
    wall = time.time() - t0
    rss = _peak_rss_bytes(proc.stderr)
    # /usr/bin/time writes its own lines to stderr (where the child's stderr also goes).
    if rss is None:
        # also check stdout just in case
        rss = _peak_rss_bytes(proc.stdout)
    status = "OK" if proc.returncode == 0 else f"FAIL({proc.returncode})"
    return {
        "backend": backend,
        "rss_mb": (rss / 1024 / 1024) if rss else None,
        "wall_s": wall,
        "status": status,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads", default="10000,100000,1000000",
                    help="comma-separated read counts")
    ap.add_argument("--backends", default="pandas,parquet",
                    help="comma-separated backends")
    ap.add_argument("--workdir", default="/tmp/crispresso_bench",
                    help="working directory for fastqs + outputs")
    ap.add_argument("--amplicon-length", type=int, default=2000)
    ap.add_argument("--n-processes", type=int, default=2)
    ap.add_argument("--timeout", type=int, default=1800,
                    help="per-run timeout (seconds)")
    ap.add_argument("--skip-1m-pandas", action="store_true",
                    help="skip the pandas backend for read counts >= 1M (it is "
                         "expected to be slow / OOM; the parquet number is the point)")
    args = ap.parse_args()

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    amplicon = _make_amplicon(args.amplicon_length)
    read_counts = [int(x) for x in args.reads.split(",")]
    backends = [b.strip() for b in args.backends.split(",")]

    rows = []
    for n in read_counts:
        fq = workdir / f"reads_{n}.fastq"
        if not fq.exists():
            print(f"[gen] {n} reads → {fq} ...", flush=True)
            sz = _generate_fastq(fq, n, amplicon)
            print(f"       {sz / 1e6:.1f} MB", flush=True)
        for backend in backends:
            if args.skip_1m_pandas and n >= 1_000_000 and backend == "pandas":
                print(f"[skip] pandas @ {n} (expected slow/OOM; --skip-1m-pandas)")
                rows.append({"n_reads": n, "backend": backend, "rss_mb": None,
                             "wall_s": None, "status": "SKIPPED"})
                continue
            print(f"[run] {backend} @ {n} reads ...", flush=True)
            out_dir = workdir / f"out_{n}_{backend}"
            if out_dir.exists():
                shutil.rmtree(out_dir)
            r = _run_one(backend, fq, amplicon, out_dir,
                         n_proc=args.n_processes, timeout_s=args.timeout)
            r["n_reads"] = n
            rss = r["rss_mb"]
            print(f"      → RSS={rss if rss is not None else 'N/A'} MB  "
                  f"wall={r['wall_s']:.1f}s  status={r['status']}", flush=True)
            rows.append(r)
            # clean up the (large) output dir to keep disk bounded
            if out_dir.exists():
                shutil.rmtree(out_dir, ignore_errors=True)

    # Report
    print("\n" + "=" * 70)
    print(f"{'reads':>10} {'backend':>9} {'peak RSS (MB)':>16} {'wall (s)':>10} {'status':>10}")
    print("-" * 70)
    for r in rows:
        rss = f"{r['rss_mb']:.1f}" if r.get("rss_mb") is not None else "—"
        wall = f"{r['wall_s']:.1f}" if r.get("wall_s") is not None else "—"
        print(f"{r['n_reads']:>10} {r['backend']:>9} {rss:>16} {wall:>10} {r['status']:>10}")
    print("=" * 70)

    # Flat-memory check for the parquet backend
    parquet = [r for r in rows if r["backend"] == "parquet" and r.get("rss_mb")]
    if len(parquet) >= 2:
        smallest = min(r["n_reads"] for r in parquet)
        largest = max(r["n_reads"] for r in parquet)
        rss_small = next(r["rss_mb"] for r in parquet if r["n_reads"] == smallest)
        rss_large = next(r["rss_mb"] for r in parquet if r["n_reads"] == largest)
        ratio = rss_large / rss_small if rss_small else float("inf")
        print(f"\nSC #1 (parquet): RSS({largest})/RSS({smallest}) = "
              f"{rss_large:.1f}/{rss_small:.1f} = {ratio:.2f}  "
              f"({'PASS' if ratio < 3 else 'FAIL — scales linearly, not flat'})")

    # Markdown summary
    md = workdir / "results.md"
    with open(md, "w") as fh:
        fh.write("# CRISPResso2 storage-backend peak-RSS benchmark\n\n")
        fh.write(f"Amplicon: synthetic {args.amplicon_length} bp, ~1:1 read diversity "
                 f"(1-3 random subs per read). Guide: {_GUID}. "
                 f"n_processes={args.n_processes}.\n\n")
        fh.write("| reads | backend | peak RSS (MB) | wall (s) | status |\n")
        fh.write("|------:|:-------|--------------:|---------:|:-------|\n")
        for r in rows:
            rss = f"{r['rss_mb']:.1f}" if r.get("rss_mb") is not None else "—"
            wall = f"{r['wall_s']:.1f}" if r.get("wall_s") is not None else "—"
            fh.write(f"| {r['n_reads']} | {r['backend']} | {rss} | {wall} | {r['status']} |\n")
        if len(parquet) >= 2:
            fh.write(f"\n**SC #1 (parquet): RSS({largest})/RSS({smallest}) = "
                     f"{ratio:.2f} — {'PASS (<3)' if ratio < 3 else 'FAIL'}**\n")
    print(f"\nmarkdown summary → {md}")


if __name__ == "__main__":
    main()
