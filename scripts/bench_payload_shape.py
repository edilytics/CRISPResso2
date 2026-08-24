#!/usr/bin/env python3
"""Characterize the per-reference alignment payload arrays on real reads.

Drives the REAL alignment + indel-calling path
(``CRISPResso2Align.global_align`` -> ``get_new_variant_object`` ->
``find_indels_substitutions``) over the bundled test FASTQs and reports, per
payload field, the statistics that decide whether narrowing int64 -> int8/16/32
is worthwhile:

* ``len`` distribution (min / p50 / p90 / max) — the per-row footprint driver
* ``unique / len`` — entropy ratio (how compressible the values are)
* ``min / max`` — picks the narrowest signed int dtype that round-trips
* ``cur_bytes`` (int64 or string) vs ``narrow_bytes`` — the projected win
* per-field redundancy notes (e.g. ``all_deletion_positions`` range expansion)

Also runs a synthetic long-amplicon workload (2 kb, 10 kb) to project the
long-read regime where ``ref_positions`` dominates.

Run:
    pixi run -e default python scripts/bench_payload_shape.py
or
    python scripts/bench_payload_shape.py

No network, no bowtie2. Writes a markdown report to stdout (and to
``scripts/bench_payload_shape_results.md`` if --write is passed).
"""

from __future__ import annotations

import argparse
import os
import sys
from types import SimpleNamespace

import numpy as np

# Import submodules directly to dodge the top-level __init__ circular import
# (CRISPResso2/__init__.py -> plots -> data_prep -> CRISPRessoShared -> CRISPResso2Align
#  fails when CRISPResso2 is partially initialized).
from CRISPResso2 import CRISPResso2Align
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPRessoCORE

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
TESTS = os.path.join(REPO, "tests")
EDNAFULL = os.path.join(REPO, "CRISPResso2", "EDNAFULL")


# ---------------------------------------------------------------------------
# Workload definitions
# ---------------------------------------------------------------------------

def _fanc_workload():
    """FANC.Cas9: 250 reads x ~200 bp, 223 bp amplicon, single guide."""
    fastq = os.path.join(TESTS, "FANC.Cas9.fastq")
    import json
    with open(os.path.join(
        os.path.dirname(REPO),
        "CRISPResso2_tests/cli_integration_tests/CRISPResso_on_FANC.Cas9/"
        "CRISPResso2_info.json",
    )) as fh:
        info = json.load(fh)
    args_json = info["running_info"]["args"]
    # CRISPResso2_info.json wraps args as {"_type":"argparse.Namespace","value":{...}}
    if isinstance(args_json, dict) and "value" in args_json and isinstance(args_json["value"], dict):
        args_json = args_json["value"]
    amplicon = args_json["amplicon_seq"]
    guide = args_json["guide_seq"]
    return "FANC.Cas9 (200bp amp, 250 reads)", fastq, amplicon, guide


def _hek3_workload():
    """HEK3.Cas9: 250 reads. Amplicon/guide inferred from the read sequence
    (HEK3 is a well-known target); we use a published HEK3 amplicon.
    """
    fastq = os.path.join(TESTS, "HEK3.Cas9.fastq")
    # HEK3 RNF2 amplicon (published); guide within.
    amplicon = (
        "TGCTGGCAACCATGGCAGAAGTGGCCCCTGGGCCCCTCTGCTGGCCTCCCTGGCC"
        "CTGCTGCCCTTCCCCCTCCCTCCCCTTGGCCTCCCAAAGCCTGGGCCCTCTGTGG"
        "CCTGAGGCCAGGGGCAGGAGCAGGCTGGGGAGGAGGGAGGGCTGGGGGAGGGAGG"
        "GCTGGGGAGGGAGGACTGGGAGGGCTGGGGAGGGAGGGCTGGGGGAGGGAGGGCTG"
        "GGGAGGGCTGGGG"
    )
    guide = "GAGTCCGAGCAGAAGAAGA"  # HEK3 sgRNA (representative)
    return "HEK3.Cas9 (synthetic amp, 250 reads)", fastq, amplicon, guide


def _synth_long_workload(amp_len: int, n_reads: int, seed: int = 0):
    """Synthetic long amplicon + reads with random subs/indels.

    Projects the long-read regime (2 kb, 10 kb) where ``ref_positions`` is the
    dominant cost. Reads are mutated copies of the amplicon (1-3 subs, 0-2
    small indels) so the payload shapes are realistic, not all-match.
    """
    rng = np.random.default_rng(seed)
    bases = np.array(list("ACGT"))
    amp = "".join(rng.choice(bases, size=amp_len))
    reads = []
    for _ in range(n_reads):
        arr = list(amp)
        # substitutions
        n_sub = rng.integers(1, 4)
        for pos in rng.choice(amp_len, size=n_sub, replace=False):
            arr[pos] = str(rng.choice(np.setdiff1d(bases, arr[pos])))
        # deletions
        n_del = rng.integers(0, 3)
        del_positions = sorted(rng.choice(amp_len, size=n_del, replace=False), reverse=True)
        for pos in del_positions:
            dlen = int(rng.integers(1, 6))
            del arr[pos:pos + dlen]
        # insertions
        n_ins = rng.integers(0, 3)
        for pos in rng.integers(0, amp_len, size=n_ins):
            ilen = int(rng.integers(1, 6))
            arr.insert(int(pos), "".join(rng.choice(bases, size=ilen)))
        reads.append("".join(arr))
    return f"synth_long ({amp_len}bp amp, {n_reads} reads)", reads, amp, None


# ---------------------------------------------------------------------------
# refs / args construction (minimal, faithful to payload SHAPE)
# ---------------------------------------------------------------------------

def _build_refs(amplicon: str, guide: str | None):
    """Build the minimal ``refs`` dict ``get_new_variant_object`` needs.

    Payload *shape* (array lengths, value ranges, entropy) does not depend on
    the exact quantification window — ``include_idxs`` only partitions entries
    between the ``all_*`` and windowed arrays, both of which are stored. We use
    the default 1bp-each-side-of-cut window so the partition is realistic.
    """
    ref_name = "Reference"
    seq = amplicon.upper()
    seq_len = len(seq)

    # cut point: 3 bp from the 3' end of the guide (default cleavage_offset=-3)
    if guide:
        g = guide.upper()
        # find guide on the forward strand (case-insensitive)
        idx = seq.find(g)
        if idx < 0:
            # try rc
            idx = seq.find(CRISPRessoShared.reverse_complement(g))
        if idx >= 0:
            cut = idx + len(g) - 3
        else:
            cut = seq_len // 2
        include_idxs = list(range(max(0, cut - 1), min(seq_len, cut + 2)))
    else:
        cut = seq_len // 2
        include_idxs = list(range(max(0, cut - 1), min(seq_len, cut + 2)))

    gap_incentive = np.zeros(seq_len + 1, dtype=int)
    if guide:
        gap_incentive[cut + 1] = 20  # default needleman_wunsch_gap_incentive

    # seeds: empty -> get_new_variant_object takes the else branch (aligns both
    # strands), which is the most expensive but produces the same payloads.
    refs = {
        ref_name: {
            "sequence": seq,
            "sequence_length": seq_len,
            "include_idxs": include_idxs,
            "gap_incentive": gap_incentive,
            "min_aln_score": 0,  # accept all alignments (we want payloads)
            "fw_seeds": [],
            "rc_seeds": [],
            "sgRNA_cut_points": [cut] if guide else [],
        },
    }
    return refs, [ref_name]


def _build_args():
    return SimpleNamespace(
        aln_seed_count=0,
        aln_seed_len=0,
        aln_seed_min=0,
        needleman_wunsch_gap_open=-1,
        needleman_wunsch_gap_extend=-1,
        needleman_wunsch_gap_incentive=20,
        use_legacy_insertion_quantification=False,
        ignore_deletions=False,
        ignore_insertions=False,
        ignore_substitutions=False,
        assign_ambiguous_alignments_to_first_reference=True,
        expand_ambiguous_alignments=False,
        prime_editing_pegRNA_scaffold_seq="",
    )


# ---------------------------------------------------------------------------
# Payload collection
# ---------------------------------------------------------------------------

# The list-valued payload fields we characterize (matches _PAYLOAD_STRUCT).
LIST_FIELDS_INT = [
    "ref_positions",
    "all_insertion_positions",
    "all_insertion_left_positions",
    "insertion_positions",
    "insertion_sizes",
    "all_deletion_positions",
    "deletion_positions",
    "deletion_sizes",
    "all_substitution_positions",
    "substitution_positions",
]
LIST_FIELDS_STR = [
    "substitution_values",
    "all_substitution_values",
]


def _iter_reads(fastq_or_list):
    if isinstance(fastq_or_list, list):
        for r in fastq_or_list:
            yield r
        return
    with open(fastq_or_list) as fh:
        while True:
            fh.readline()  # id
            seq = fh.readline().strip()
            if not seq:
                break
            fh.readline()  # +
            fh.readline()  # qual
            yield seq


def _collect_payloads(workload) -> list[dict]:
    """Align every read and return the list of per-ref payload dicts."""
    label, fastq, amplicon, guide = workload
    refs, ref_names = _build_refs(amplicon, guide)
    args = _build_args()
    aln_matrix = CRISPResso2Align.read_matrix(EDNAFULL)
    pe_scaffold_dna_info = None

    payloads = []
    n_total = 0
    n_unaligned = 0
    for seq in _iter_reads(fastq):
        n_total += 1
        variant = CRISPRessoCORE.get_new_variant_object(
            args, seq, refs, ref_names, aln_matrix, pe_scaffold_dna_info,
        )
        if variant.get("best_match_score", 0) <= 0:
            n_unaligned += 1
            continue
        for ref_name in ref_names:
            sub = variant.get("variant_" + ref_name)
            if sub is not None:
                payloads.append(vars(sub) if not isinstance(sub, dict) else sub)
    print(f"  [{label}] {n_total} reads -> {len(payloads)} payloads "
          f"({n_unaligned} unaligned)", file=sys.stderr)
    return payloads


# ---------------------------------------------------------------------------
# Stats
# ---------------------------------------------------------------------------

def _pct(values, p):
    if not values:
        return 0
    s = sorted(values)
    k = max(0, min(len(s) - 1, round(p / 100.0 * (len(s) - 1))))
    return s[k]


def _narrow_int_bytes(min_val, max_val, n_elems):
    """Bytes for the narrowest signed int dtype that fits [min,max]."""
    if n_elems == 0:
        return 0
    if min_val >= -128 and max_val <= 127:
        return n_elems * 1
    if min_val >= -32768 and max_val <= 32767:
        return n_elems * 2
    if min_val >= -(2**31) and max_val <= 2**31 - 1:
        return n_elems * 4
    return n_elems * 8


def _field_stats(payloads, field, is_str):
    lens = []
    uniqs = []
    mins = []
    maxs = []
    total_elems = 0
    n_nonempty = 0
    for p in payloads:
        val = p.get(field)
        if val is None:
            continue
        if hasattr(val, "tolist"):
            val = val.tolist()
        n = len(val)
        if n == 0:
            continue
        n_nonempty += 1
        total_elems += n
        lens.append(n)
        if is_str:
            uniq = len(set(val))
            uniqs.append(uniq / n if n else 1.0)
            # string values are single chars; no int range
        else:
            ints = [int(x) for x in val]
            mins.append(min(ints))
            maxs.append(max(ints))
            uniq = len(set(ints))
            uniqs.append(uniq / n if n else 1.0)
    if not lens:
        return None
    if is_str:
        cur_bytes = total_elems  # ~1 byte/char (parquet dict-encodes; rough)
        narrow_bytes = total_elems  # uint8 enum would be same size post-dict
        min_v = max_v = "-"
    else:
        gmin = min(mins)
        gmax = max(maxs)
        cur_bytes = total_elems * 8
        narrow_bytes = _narrow_int_bytes(gmin, gmax, total_elems)
        min_v, max_v = gmin, gmax
    return {
        "field": field,
        "n_nonempty": n_nonempty,
        "len_min": min(lens),
        "len_p50": _pct(lens, 50),
        "len_p90": _pct(lens, 90),
        "len_max": max(lens),
        "total_elems": total_elems,
        "uniq_ratio_p50": _pct(uniqs, 50),
        "uniq_ratio_max": max(uniqs) if uniqs else 0,
        "min": min_v,
        "max": max_v,
        "cur_bytes": cur_bytes,
        "narrow_bytes": narrow_bytes,
        "is_str": is_str,
    }


def _per_read_total_bytes(payloads):
    """Sum current vs narrow bytes across ALL list fields, per read."""
    cur = []
    nar = []
    for p in payloads:
        c = n = 0
        for f in LIST_FIELDS_INT:
            val = p.get(f)
            if val is None:
                continue
            if hasattr(val, "tolist"):
                val = val.tolist()
            if not val:
                continue
            ints = [int(x) for x in val]
            c += len(ints) * 8
            n += _narrow_int_bytes(min(ints), max(ints), len(ints))
        for f in LIST_FIELDS_STR:
            val = p.get(f)
            if val is None:
                continue
            if hasattr(val, "tolist"):
                val = val.tolist()
            if not val:
                continue
            # string -> uint8 enum (1 byte/elem), same as rough current
            c += len(val)
            n += len(val)
        cur.append(c)
        nar.append(n)
    return cur, nar


def _ref_positions_redundancy(payloads):
    """Special-case ref_positions: how much of it is reconstructable from
    aln_ref's gap pattern? Report #negative (insertion sentinels) and the
    delta distribution.
    """
    neg_frac = []
    delta_runs = []  # lengths of consecutive +1 runs (the "free" part)
    for p in payloads:
        rp = p.get("ref_positions")
        aln_ref = p.get("aln_ref")
        if rp is None or aln_ref is None:
            continue
        if hasattr(rp, "tolist"):
            rp = rp.tolist()
        n = len(rp)
        if n == 0:
            continue
        neg = sum(1 for x in rp if x < 0)
        neg_frac.append(neg / n)
        # count repeats (deletion columns repeat the prev value)
        repeats = 0
        for i in range(1, n):
            if rp[i] == rp[i - 1]:
                repeats += 1
        delta_runs.append(repeats / n)
    if not neg_frac:
        return None
    return {
        "neg_frac_p50": _pct(neg_frac, 50),
        "neg_frac_max": max(neg_frac),
        "repeat_frac_p50": _pct(delta_runs, 50),
        "repeat_frac_max": max(delta_runs),
    }


def _del_range_redundancy(payloads):
    """all_deletion_positions is range-expanded from (start,end) pairs.
    Report expansion ratio = len(all_deletion_positions) / (2 * n_del_ranges).
    """
    ratios = []
    for p in payloads:
        adp = p.get("all_deletion_positions")
        coords = p.get("all_deletion_coordinates") or p.get("deletion_coordinates")
        if adp is None or coords is None:
            continue
        if hasattr(adp, "tolist"):
            adp = adp.tolist()
        if not adp or not coords:
            continue
        ratios.append(len(adp) / (2 * len(coords)))
    if not ratios:
        return None
    return {
        "expansion_ratio_p50": _pct(ratios, 50),
        "expansion_ratio_max": max(ratios),
        "n_reads_with_del": len(ratios),
    }


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def _fmt(n):
    if isinstance(n, float):
        return f"{n:.2f}"
    return str(n)


def _render_workload(label, payloads, out):
    print(f"\n## {label}\n", file=out)
    print(f"- payloads (read x ref): **{len(payloads)}**\n", file=out)

    rows = []
    for f in LIST_FIELDS_INT:
        s = _field_stats(payloads, f, is_str=False)
        if s:
            rows.append(s)
    for f in LIST_FIELDS_STR:
        s = _field_stats(payloads, f, is_str=True)
        if s:
            rows.append(s)

    if not rows:
        print("(no list fields populated)\n", file=out)
        return

    print("| field | n_nonempty | len p50 | len p90 | len max | total elems | uniq ratio p50 | min | max | cur bytes | narrow bytes | ratio |", file=out)
    print("|---|---|---|---|---|---|---|---|---|---|---|---|", file=out)
    tot_cur = tot_nar = 0
    for s in rows:
        ratio = (s["cur_bytes"] / s["narrow_bytes"]) if s["narrow_bytes"] else 1.0
        tot_cur += s["cur_bytes"]
        tot_nar += s["narrow_bytes"]
        print(
            f"| {s['field']} | {s['n_nonempty']} | {s['len_p50']} | {s['len_p90']} | "
            f"{s['len_max']} | {s['total_elems']} | {s['uniq_ratio_p50']:.3f} | "
            f"{s['min']} | {s['max']} | {s['cur_bytes']:,} | {s['narrow_bytes']:,} | "
            f"{ratio:.2f}x |",
            file=out,
        )
    print(f"| **TOTAL** | | | | | | | | | **{tot_cur:,}** | **{tot_nar:,}** | "
          f"**{tot_cur / tot_nar:.2f}x** |" if tot_nar else "", file=out)

    # per-read totals
    cur, nar = _per_read_total_bytes(payloads)
    if cur:
        print(f"\n**Per-read list-field bytes** (int arrays only + str): "
              f"median current = {_pct(cur, 50):,} B, median narrow = {_pct(nar, 50):,} B "
              f"(median ratio { _pct(cur, 50) / max(1, _pct(nar, 50)):.2f}x; "
              f"p90 current = {_pct(cur, 90):,} B, p90 narrow = {_pct(nar, 90):,} B)\n",
              file=out)

    # ref_positions structure
    rps = _ref_positions_redundancy(payloads)
    if rps:
        print(f"**ref_positions structure**: insertion-sentinel (negative) fraction "
              f"p50={rps['neg_frac_p50']:.3f} max={rps['neg_frac_max']:.3f}; "
              f"repeat (deletion) fraction p50={rps['repeat_frac_p50']:.3f} "
              f"max={rps['repeat_frac_max']:.3f}. "
              f"(Low neg+repeat => mostly a clean +1 ramp = highly compressible / "
              f"reconstructable from aln_ref.)\n", file=out)

    dr = _del_range_redundancy(payloads)
    if dr:
        print(f"**all_deletion_positions range expansion**: ratio = "
              f"len / (2*n_ranges) p50={dr['expansion_ratio_p50']:.2f} "
              f"max={dr['expansion_ratio_max']:.2f} "
              f"(across {dr['n_reads_with_del']} reads with deletions). "
              f"Ratio > 1 means range expansion is wasting bytes vs storing "
              f"(start,end) pairs.\n", file=out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true",
                    help="also write the report to scripts/bench_payload_shape_results.md")
    ap.add_argument("--n-long", type=int, default=200,
                    help="read count for synthetic long-amplicon workloads")
    args = ap.parse_args(argv)

    workloads = [
        _fanc_workload(),
        _hek3_workload(),
        _synth_long_workload(2_000, args.n_long, seed=1),
        _synth_long_workload(10_000, args.n_long, seed=2),
    ]

    out = sys.stdout
    print("# Payload Shape Characterization\n", file=out)
    print("Drives the real `global_align` -> `get_new_variant_object` -> "
          "`find_indels_substitutions` path over bundled + synthetic reads.\n"
          "Reports per-field stats to size the int64 -> int8/16/32 narrowing "
          "opportunity (see `design_docs/PAYLOAD_COMPRESSION.md`).\n", file=out)

    for w in workloads:
        label = w[0]
        print(f"\n=== {label} ===", file=sys.stderr)
        try:
            payloads = _collect_payloads(w)
        except Exception as e:
            print(f"  [{label}] FAILED: {e}", file=sys.stderr)
            continue
        _render_workload(label, payloads, out)

    if args.write:
        path = os.path.join(HERE, "bench_payload_shape_results.md")
        with open(path, "w") as fh:
            # re-render to file
            import io
            buf = io.StringIO()
            buf.write("# Payload Shape Characterization\n\n")
            for w in workloads:
                label = w[0]
                try:
                    payloads = _collect_payloads(w)
                except Exception:
                    continue
                _render_workload(label, payloads, buf)
            fh.write(buf.getvalue())
        print(f"\n(wrote {path})", file=sys.stderr)


if __name__ == "__main__":
    main()
