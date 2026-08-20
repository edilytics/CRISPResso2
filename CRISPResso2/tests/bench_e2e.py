#!/usr/bin/env python3
"""End-to-end CRISPResso benchmark: wall time + peak RAM, A/B current vs master.

Unlike ``CRISPResso2/tests/bench_align.py`` (which times only the ``global_align`` kernel),
this runs the FULL ``CRISPResso`` Core pipeline so the optimization's effect --
and its RAM footprint -- can be seen in a realistic run, not just the microbench.

It generates a deterministic synthetic FASTQ derived from the FANC amplicon
(CRISPR-like indels near the cut + per-read sequencing noise -> many unique
alleles, so per-read alignment is exercised at realistic volume), then runs
``CRISPResso`` on it N times per version (current branch and genuine master),
interleaved round-robin to smooth system drift. Each run is wrapped in
``/usr/bin/time -l`` to capture wall/user/sys time and peak RSS.

Two configs are reported:
  - compute : --suppress_plots --suppress_report  (isolates the alignment +
             processing core; the optimization's signal is cleanest here)
  - full     : normal pipeline incl. plots/report  (realistic per-run cost)

The A/B delta is compared against the within-version spread (noise floor), like
``CRISPResso2/tests/compare_bench.py``.

PREREQ: a built+installed master worktree (see CRISPResso2/tests/ab_vs_master.py docstring):
  git worktree add /tmp/c2-master master
  cd /tmp/c2-master && pixi install -e test && pixi run -e test pip install -e . --no-build-isolation

USAGE
-----
  pixi run -e test python CRISPResso2/tests/bench_e2e.py                 # default: 12k reads, 3 reps, both configs
  pixi run -e test python CRISPResso2/tests/bench_e2e.py --reads 20000 --reps 5 --config compute
"""
import argparse
import os
import re
import shlex
import statistics
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))
CURRENT_CRISPRESSO = os.path.join(REPO_ROOT, ".pixi", "envs", "test", "bin", "CRISPResso")
DEFAULT_MASTER = "/tmp/c2-master"
MASTER_CRISPRESSO = os.path.join(DEFAULT_MASTER, ".pixi", "envs", "test", "bin", "CRISPResso")
FASTQ = os.path.join(REPO_ROOT, "tests", "FANC.Cas9.fastq")

AMP = ("CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCA"
       "CCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCT"
       "ACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATC"
       "CATCGGCGCTTTGGTCGG")
GUIDE = "GGAATCCCTTCTGCAGCACC"


# --------------------------------------------------------------------------- #
# Deterministic synthetic FASTQ (alignment-heavy, realistic CRISPR edit model)
# --------------------------------------------------------------------------- #
import random


def find_cut(amp, guide):
    """Cas9 cut ~3bp upstream of the 3' end of the guide match."""
    idx = amp.find(guide)
    if idx < 0:
        return len(amp) // 2
    return idx + len(guide) - 6  # ~3bp before PAM-proximal end


def make_amplicon(length, seed):
    """Synthesize a random amplicon of the given length with a guide embedded.

    Used for long-read benchmarks (e.g. 5 kb) where a real amplicon of that
    size isn't shipped. Biology is irrelevant here -- we only need a reference
    sequence + guide so CRISPResso exercises global_align at large N*M.
    """
    rng = random.Random(seed ^ 0x5EED)
    amp = "".join(rng.choice("ACGT") for _ in range(length))
    guide = "GGAATCCCTTCTGCAGCACC"
    pos = length // 2
    amp = amp[:pos] + guide + amp[pos + len(guide):]  # overwrite mid with guide
    return amp, guide


def generate_fastq(path, n_reads, seed, amp, guide):
    """Deterministic synthetic FASTQ maximising UNIQUE variants.

    CRISPResso deduplicates reads (variantCache), so alignment work scales with
    the number of UNIQUE alleles, not raw read count. To make per-read
    alignment a dominant, measurable phase, this emits near-unique reads: a
    realistic CRISPR edit model (wildtype / indel near cut / substitution) PLUS
    enough per-read sequencing noise that ~all reads are distinct.
    """
    rng = random.Random(seed)
    cut = find_cut(amp, guide)
    n = len(amp)
    with open(path, "w") as fh:
        for i in range(n_reads):
            seq = list(amp)
            r = rng.random()
            if r < 0.20:                      # wildtype (unedited)
                pass
            elif r < 0.85:                    # indel near the cut (NHEJ-like)
                pos = max(1, min(n - 2, cut + rng.randint(-3, 3)))
                if rng.random() < 0.7:        # deletion
                    dlen = rng.randint(1, 12)
                    del seq[pos:pos + dlen]
                else:                         # insertion
                    ilen = rng.randint(1, 6)
                    seq[pos:pos] = list("".join(rng.choice("ACGT") for _ in range(ilen)))
            else:                             # substitution near cut
                pos = max(0, min(len(seq) - 1, cut + rng.randint(-3, 3)))
                seq[pos] = rng.choice([b for b in "ACGT" if b != seq[pos]])
            # always add 2-3 random subs -> ~all reads unique -> alignment
            # work scales linearly with read count (no dedup collapse)
            for _ in range(rng.randint(2, 3)):
                p = rng.randrange(len(seq))
                seq[p] = rng.choice("ACGT")
            s = "".join(seq)
            fh.write(f"@bench_read_{i}\n{s}\n+\n{'I' * len(s)}\n")
    return n_reads


# --------------------------------------------------------------------------- #
# Measurement
# --------------------------------------------------------------------------- #
def run_once(crispresso, fastq, out_dir, name, suppress, amp, guide):
    """Run CRISPResso once; return dict with wall/user/sys (s) and peak_rss_mb.

    Wraps in /usr/bin/time -l; CRISPResso's own output is sent to /dev/null so
    /usr/bin/time's stderr summary is clean and parseable.
    """
    flags = []
    if suppress:
        flags += ["--suppress_plots", "--suppress_report"]
    inner = " ".join(shlex.quote(x) for x in [
        crispresso, "-r1", fastq, "-a", amp, "-g", guide,
        "-n", name, "-o", out_dir] + flags) + " >/dev/null 2>&1"
    t0 = time.monotonic()
    proc = subprocess.run(["/usr/bin/time", "-l", "sh", "-c", inner],
                          capture_output=True, text=True, check=False)
    wall_mono = time.monotonic() - t0
    if proc.returncode != 0:
        return {"failed": True, "returncode": proc.returncode,
                "stderr": proc.stderr[-400:]}
    m = re.search(r"([\d.]+)\s*real\s+([\d.]+)\s*user\s+([\d.]+)\s*sys",
                  proc.stderr)
    rss = re.search(r"(\d+)\s+maximum resident set size", proc.stderr)
    if not m or not rss:
        return {"failed": True, "returncode": proc.returncode,
                "stderr": "unparseable /usr/bin/time output: " + proc.stderr[-400:]}
    return {
        "wall": float(m.group(1)), "user": float(m.group(2)),
        "sys": float(m.group(3)), "wall_mono": wall_mono,
        "peak_rss_mb": int(rss.group(1)) / (1024 * 1024),
    }


def stats(values):
    if not values:
        return {}
    s = sorted(values)

    def pct(p):
        i = min(len(s) - 1, max(0, round((p / 100.0) * (len(s) - 1))))
        return s[i]
    return {"min": s[0], "median": statistics.median(values),
            "mean": statistics.mean(values), "p95": pct(95), "max": s[-1]}


def fmt_stats(st, unit):
    return (f"min={st['min']:.2f} med={st['median']:.2f} "
            f"mean={st['mean']:.2f} p95={st['p95']:.2f} {unit}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads", type=int, default=12000)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--config", choices=["compute", "full", "both"], default="both")
    ap.add_argument("--amplicon-len", type=int, default=0,
                    help="synthesize a random amplicon of this length (e.g. 5000 "
                         "for long-read); 0 = use the real ~220bp FANC amplicon")
    ap.add_argument("--master", default=DEFAULT_MASTER)
    ap.add_argument("--seed", type=int, default=20260730)
    args = ap.parse_args(argv)

    versions = [("current", CURRENT_CRISPRESSO), ("master", MASTER_CRISPRESSO)]
    for name, exe in versions:
        if not os.path.exists(exe):
            sys.exit(f"ERROR: {name} CRISPResso not found at {exe}")

    print(f"current HEAD: {subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=REPO_ROOT, text=True).strip()}")
    print(f"master  HEAD: {subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=args.master, text=True).strip()}")

    with tempfile.TemporaryDirectory(prefix="c2bench_") as tmp:
        if args.amplicon_len > 0:
            amp, guide = make_amplicon(args.amplicon_len, args.seed)
            print(f"amplicon: synthesized {len(amp)} bp (long-read mode)")
        else:
            amp, guide = AMP, GUIDE
            print(f"amplicon: real FANC {len(amp)} bp")
        fastq = os.path.join(tmp, "synth.fastq")
        print(f"generating synthetic FASTQ: {args.reads} reads -> {fastq}")
        generate_fastq(fastq, args.reads, args.seed, amp, guide)
        print(f"  size: {os.path.getsize(fastq) / 1024:.0f} KB")

        configs = (["compute", "full"] if args.config == "both"
                   else [args.config])
        for cfg in configs:
            suppress = (cfg == "compute")
            print("\n" + "=" * 78)
            print(f"CONFIG: {cfg} "
                  f"({'--suppress_plots --suppress_report' if suppress else 'full pipeline incl. plots'})")
            print("=" * 78)
            # round-robin across versions to smooth drift
            results = {name: [] for name, _ in versions}
            failures = {name: [] for name, _ in versions}
            for rep in range(args.reps):
                for name, exe in versions:
                    run_out = os.path.join(tmp, f"{cfg}_{name}_{rep}")
                    os.makedirs(run_out, exist_ok=True)
                    res = run_once(exe, fastq, run_out, f"r{rep}", suppress, amp, guide)
                    if res.get("failed"):
                        failures[name].append(res)
                        print(f"  rep {rep} {name:<8}: FAILED rc={res['returncode']} "
                              f"{res.get('stderr', '')[:120]}")
                    else:
                        results[name].append(res)
                        print(f"  rep {rep} {name:<8}: wall={res['wall']:.2f}s "
                              f"(cpu={res['user'] + res['sys']:.2f}s) "
                              f"peakRSS={res['peak_rss_mb']:.0f} MB")

            print("\n-- aggregate (per version) --")
            agg = {}
            for name, _ in versions:
                runs = results[name]
                if not runs:
                    print(f"  {name}: NO SUCCESSFUL RUNS")
                    continue
                w = stats([r["wall"] for r in runs])
                c = stats([r["user"] + r["sys"] for r in runs])
                rss = max(r["peak_rss_mb"] for r in runs)
                agg[name] = {"wall_med": w["median"], "wall_min": w["min"],
                             "wall_spread": w["max"] - w["min"],
                             "cpu_med": c["median"], "rss_peak": rss}
                print(f"  {name:<8} WALL : {fmt_stats(w, 's')}")
                print(f"  {'':8} CPU  : {fmt_stats(c, 's')}")
                print(f"  {'':8} peakRSS (max of reps): {rss:.0f} MB")

            if "current" in agg and "master" in agg:
                print("\n-- A/B (current vs master) --")
                cur, mas = agg["current"], agg["master"]
                d_wall = mas["wall_med"] - cur["wall_med"]
                pct = 100.0 * d_wall / mas["wall_med"] if mas["wall_med"] else 0
                noise = max(cur["wall_spread"], mas["wall_spread"])
                d_rss = mas["rss_peak"] - cur["rss_peak"]
                verdict = ("real" if abs(d_wall) > 2 * noise and noise > 0
                           else ("noise-level" if noise > 0 else "n/a"))
                print(f"  wall : master {mas['wall_med']:.2f}s -> current {cur['wall_med']:.2f}s "
                      f"= {d_wall:+.2f}s ({pct:+.1f}%)   [noise floor ~{noise:.2f}s -> {verdict}]")
                print(f"  cpu  : master {mas['cpu_med']:.2f}s -> current {cur['cpu_med']:.2f}s "
                      f"= {cur['cpu_med'] - mas['cpu_med']:+.2f}s")
                print(f"  RAM  : master {mas['rss_peak']:.0f} MB -> current {cur['rss_peak']:.0f} MB "
                      f"= {d_rss:+.0f} MB")

            if any(failures.values()):
                print("\n!! some runs failed; see above. A/B may be incomplete.")


if __name__ == "__main__":
    main()
