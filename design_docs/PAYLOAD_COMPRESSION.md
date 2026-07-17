# Payload Array Compression — Investigation & Plan

**Status:** Steps 3 (int16 narrowing) + 4 (all_deletion_coordinates compact storage) IMPLEMENTED + parity-tested. Step 5 deferred.
**Scope:** The per-reference alignment payload stored in parquet shards
(`aligned_{i}.parquet` → `COLLAPSED_SCHEMA`) and streamed between Stage 2
(workers) → Stage 3 (collapse) → Stage 4 (count-vector stream-out). This is
the "payload for the alignment" the user is asking about.

---

## 1. What the payload actually contains

The per-reference struct (`_PAYLOAD_STRUCT` in `storage.py`) carries ~30
fields. The list-valued ones (the ones worth compressing) are:

| field | length | dtype | built by | consumed by |
|---|---|---|---|---|
| `ref_positions` | **alignment length L** | int64 | `find_indels_substitutions` | per-base count vectors (`all_base_count_vectors`) — `CRISPRessoCORE.py:4334` |
| `all_substitution_values` | **#subs** (NOT L) | string | appended only at substitution sites | per-base sub vectors — `CRISPRessoCORE.py:4311` |
| `all_substitution_positions` | #subs | int64 | same | sub count vectors — `:4300` |
| `all_deletion_positions` | **sum of del sizes** (range-expanded) | int64 | `extend(range(s,e))` | del count vectors — `:4283` |
| `all_insertion_positions` | 2 × #ins (start,end pairs) | int64 | appended as pairs | ins count vectors |
| `deletion_coordinates` / `insertion_coordinates` | #del / #ins | struct(s,e) | — | length/exon vectors |
| `substitution_values` | #window-subs | string | — | detailed TSV |

The dominant arrays by bytes are **`ref_positions`** (length = alignment
length, 8 bytes/elem) and, for reads with large deletions, **`all_deletion_positions`** (range-expanded).

---

## 2. Key findings (some contradict the existing comments)

### Finding 1 — `all_substitution_values` is NOT amplicon-length.
The `storage.py` docstring (`_PAYLOAD_STRUCT` comment, lines 607–614) claims it
is "per-position substituted base across the WHOLE amplicon (length = ref_len;
'-' where not substituted)." **This is wrong.** Both `find_indels_substitutions`
and the `_legacy` variant append to it *only at substitution sites*
(`CRISPRessoCOREResources.pyx:113` / `:224`), so `len(all_substitution_values)
== len(all_substitution_positions) == #subs`. It is already sparse/compact.
Values are single chars drawn from `{A,C,G,T}` (the code excludes `-` and `N`).
=> Generic compression here is low-value. The actionable fix is (a) correct the
misleading comment, and (b) optionally store as `uint8` enum (4 values) instead
of `string` — cheap, but small absolute win because the array is already short.

### Finding 2 — `ref_positions` is the real cost, and it is highly structured.
It is a monotonic-ish int stream over the alignment columns:
- `+1` step at every match/mismatch (non-gap ref char),
- repeated value at deletion columns (read char `-`),
- a negative sentinel (`-1` or `-idx`) at insertion columns (ref char `-`).

So its *information* content is tiny (essentially "where are the gaps"), but it
is stored as a full int64 array of length L. For a 200 bp amplicon that is
1.6 KB/read; for a 10 kb long-read amplicon it is 80 KB/read — and at a 50 k-row
adaptive batch that is the ~4 GB in-flight footprint the `_adaptive_batch_size`
machinery is already bending over backwards to avoid.

It is consumed in exactly one place downstream: the per-base count-vector loop
(`CRISPRessoCORE.py:4334`):
```python
for i in range(len(aln_seq)):
    if ref_pos[i] < 0: continue
    all_base_count_vectors[ref_name + "_" + aln_seq[i]][ref_pos[i]] += variant_count
```
That loop only needs (a) which columns are insertions (ref_pos < 0), and (b) the
running ref index per column — both derivable from `aln_ref`'s gap pattern. So
`ref_positions` is *reconstructable* from `aln_ref` alone; it does not need to be
persisted at all.

### Finding 3 — `all_deletion_positions` is a redundant range expansion.
It is built by `all_deletion_positions.extend(range(start, end))` for each
deletion, then consumed by numpy fancy-indexing
(`all_deletion_count_vectors[ref][all_deletion_positions] += count`). It is
exactly the flattened expansion of `(start,end)` pairs — i.e. it is fully
described by a small list of `(start, size)` pairs (which is what
`deletion_coordinates` + `deletion_sizes` already are, but window-scoped). The
"all_" (whole-amplicon) variant currently has no coordinate-pair equivalent
stored (`all_deletion_coordinates` is deliberately dropped, per the
`_PAYLOAD_STRUCT` comment). For a 100 bp deletion that is 100 int64s (800 bytes)
replaceable by 2 int16s (4 bytes).

### Finding 4 — Generic Huffman on top of parquet is the wrong tool.
Parquet already applies dictionary encoding + Snappy/Zstd per page. Stacking a
custom Huffman/entropy layer on top would (a) double-compress (CPU waste, no
size win — Snappy already exploits the byte redundancy), and (b) break columnar
pushdown / `iter_batches` projection, which the streaming design depends on.
The user's *instinct* (lots of duplicates / sparsity) is correct; the right
expression of it is **schema-level redundancy elimination + narrower physical
types**, not a generic entropy codec.

---

## 3. Opportunities, ranked by leverage

| # | change | lever | expected win | cost/risk |
|---|---|---|---|---|
| **A** | Drop `ref_positions` from the persisted payload; reconstruct on read from `aln_ref` (or keep an int8/16 gap-mask). | eliminates the single largest array | ~8–10× reduction in per-row payload bytes for long reads; biggest batch-size / RSS win | parity: reconstruction must match `find_indels_substitutions` exactly (incl. `-1` vs `-idx` sentinel + the `idx==0` edge). Needs a parity test vs the Cython output. |
| **B** | Narrow int64 → int8/int16 for the position arrays that remain (`all_*_positions`, `insertion_*`). Values are bounded by ref_len; delta-encode `ref_positions` if kept. | 4–8× physical shrink of remaining int arrays | moderate; compounds with A | arrow `list<int16>` round-trips fine; numpy fancy-indexing needs `.astype(np.intp)` on read (cheap). |
| **C** | Store `all_deletion_positions` as `(start,size)` pairs (RLE/range form) and expand on read at the single fancy-indexing call site. | removes range expansion — big for large deletions | large for large-del regimes; ~0 for small dels | one-line expand at `:4283`; parity test. |
| **D** | `all_substitution_values`/`substitution_values`: string → `uint8` enum (A/C/G/T/N → 0–4). | removes per-elem string overhead | small (array is already short) | trivial; do for tidiness, not perf. |
| **E** | Correct the misleading `storage.py` docstring (Finding 1). | docs accuracy | — | none. |

A+B together attack the dominant term; C attacks the second; D/E are cleanup.

---

## 4. Validated results (Step 3 — int16 narrowing)

Characterization script: `scripts/bench_payload_shape.py` (results in
`scripts/bench_payload_shape_results.md`). Drove the real
`global_align` → `get_new_variant_object` → `find_indels_substitutions`
path over FANC.Cas9, HEK3.Cas9, and synthetic 2 kb / 10 kb amplicons.

**Narrowing is a uniform 4× win on the in-flight arrow footprint** (the RSS
driver that `_adaptive_batch_size` targets):

| workload | per-row int-list bytes (int64) | (int16) | ratio |
|---|---|---|---|
| FANC.Cas9 (223 bp amp) | 6,332 B | 1,583 B | 4.0× |
| HEK3.Cas9 (~230 bp amp) | 6,594 B | 1,525 B | 4.3× |
| synth 2 kb amp | 16,144 B | 4,036 B | 4.0× |
| synth 10 kb amp | 80,144 B | 20,036 B | 4.0× |

Projected for a 50 k-row adaptive batch (FANC shape): **317 MB → 79 MB** freed
in the int list columns alone.

**On-disk reduction is negligible (~1.01×)** — parquet's Snappy + dictionary
encoding already crushes the repetitive int64 arrays (monotonic ramps,
repeated deletion values). The win is purely in-RSS during streaming, which is
exactly what the flat-memory design exists to bound. This is why a generic
Huffman layer would have been wasted effort (Finding 4): parquet's existing
entropy stage already captures the on-disk redundancy; the untapped lever was
the *physical element width*, not the entropy.

**`.index()` safety (your concern, confirmed):** pyarrow deserializes
`list<int16>` cells to plain Python `list[int]` — values are regular ints, so
`ref_positions.index(cut_point)` (called in ~15 places: plot code,
`writers/vcf.py`, prime-editing scaffold lookup, insertion-coordinate
remapping) works identically. Verified empirically across int8/int16/int32/
int64 — all produce the same Python list on read-back. The narrowing is purely
an on-disk + in-arrow physical encoding; invisible to every consumer.

**Guardrail:** `_to_int_list` validates the value range before casting and
raises a clear `ValueError` (naming the offending value + the amplicon-
bounded constraint) if any value exceeds int16 range (±32 767). pyarrow would
also raise `ArrowInvalid` on overflow, but the guardrail produces a better
error. The parquet backend is amplicon-bounded by design (amplicons up to
~32 kb); larger workloads use the default pandas backend.

**Latent bug noted (out of scope, pre-existing):** two non-`get_slice`
read-back paths (`_row_to_variant_payload` ~line 3168, and the
`_allele_row_to_*` path ~line 3537) reconstruct `ref_positions` via
`_to_np_int_array` → numpy array, which would break `.index()`. The plot
path uses `get_slice` / `_reconstruct_slice_cell` (correctly returns a list)
so it's unaffected today. Worth fixing when those paths are wired into `main`,
but orthogonal to narrowing (they'd break at int64 too).

## 5. Proposed approach (incremental, parity-gated)

Each step lands behind `--storage_backend parquet` (default pandas path is the
parity oracle, unchanged) and is gated by a unit test that diffs the
reconstructed payload against the current `find_indels_substitutions` output on
canned + real alignments.

1. **Step 1 — Characterize (DONE).** `scripts/bench_payload_shape.py` ran
   over FANC/HEK3/synth workloads; results above validate the 4× narrowing
   leverage and the `all_deletion_positions` range-expansion redundancy on
   real data.

2. **Step 2 — Fix docs + cheap wins (D, E).** Correct the
   `all_substitution_values` comment (Finding 1: it is NOT amplicon-length);
   switch the two string position-value arrays to `uint8` enum. Low risk.

3. **Step 3 — Narrow int types (DONE).** `_to_int_list` now produces int16
   numpy arrays; both schemas (`_PAYLOAD_STRUCT`, `COLLAPSED_SCHEMA`) use
   `pa.list_(pa.int16())` for all position/size columns. Read-back unchanged.
   4× in-flight reduction confirmed. Tests: `test_narrow_int16_schema_and_round_trip`,
   `test_narrow_int16_guardrail_rejects_overflow`, updated
   `test_dtype_count_int64_positions_narrow`; full 891-test unit suite green.

4. **Step 4 — Range-pair encoding for `all_deletion_positions` (DONE, aligned-shard scoped).**
   The aligned shard (`_PAYLOAD_STRUCT`) now stores `all_deletion_coordinates`
   (compact, 2 int16s/deletion) instead of the range-expanded
   `all_deletion_positions` (N int16s/deletion). On FANC.Cas9 this is a **110×
   reduction** in the deletion column's in-flight arrow footprint (51,051
   expanded int16s → 462 coord tuples). `_struct_to_payload` reconstructs the
   expanded form via `_expand_deletion_coords` at read-back so the payload dict
   matches the pandas-path shape — zero downstream changes.

   **Scope decision:** only the aligned shard (per-read, highest cardinality +
   the in-flight arrow RSS driver) was compacted. The collapsed parquet
   (per-allele, low cardinality) keeps `all_deletion_positions` because
   compacting it would require an expand/compress dance across the
   parity-tested allele-row boundary (`_assert_rows_equal` enforces exact
   key-set equality with the pandas `get_allele_row` replica) — not worth the
   complexity for the small per-allele RSS term. The aligned shard captures
   essentially all the win.

   The producer's `all_deletion_coordinates` is used directly when present
   (always, for real payloads); a `_compress_deletion_positions` RLE fallback
   derives it from `all_deletion_positions` for canned test fixtures. Tests:
   `test_all_deletion_coordinates_compact_storage_round_trip` (expand/compress
   inverse, multi-deletion, trailing-deletion semantics, derive-fallback);
   full 892-test unit suite green.

5. **Step 5 — Drop / reconstruct `ref_positions` (A) — DEFERRED.**
   `ref_positions` is used via `.index()` in ~15 places (plot code, VCF
   writer, prime-editing scaffold lookup, insertion-coordinate remapping) —
   it's a ref-coord → alignment-column lookup table, not just a per-base
   vector input. The bench shows for the short-amplicon heavy-editing regime
   (FANC/HEK3) 50% of `ref_positions` is insertion sentinels + 50% deletion
   repeats, so reconstruction from `aln_ref` is non-trivial there (though
   trivial for the long-read regime where it's a ~99% clean +1 ramp).
   Narrowing (Step 3) already captures the 4× win with zero risk, so dropping
   is not worth the `.index()` reimplementation cost now. Revisit only if
   the long-read RSS budget demands more than 4×.

---

## 6. Success criteria (Step 3)

- **Parity:** all 101 `test_storage.py` tests + full 891-test unit suite
  pass. Reconstructed payloads value-equal the int64 path (round-trip test
  covers negative sentinels, repeats, and int16 boundary values 32000/-32000).
- **Memory:** 4× reduction in in-flight arrow int-list footprint confirmed on
  FANC/HEK3/2 kb/10 kb workloads (317 MB → 79 MB projected for a 50 k-row
  FANC-shape batch).
- **No regression** on the default `pandas` storage path (untouched).

## 7. Open questions / follow-ups

1. **Step 4 leverage check:** the 55× `all_deletion_positions` expansion on
   FANC is the single biggest remaining redundancy. Worth doing next?
2. Should the `_COORD_STRUCT` (start/end) fields also narrow to int16? Left
   int64 for now (small count per read; not the dominant cost).
3. Fix the pre-existing `_to_np_int_array` → numpy paths (latent `.index()`
   bug) when those read-back paths are wired into `main`.
