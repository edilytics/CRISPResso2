# Multi-Mode Plot Plugin & Data Prep Extraction

## Overview

Bring the 5 remaining CRISPResso2 modes (Batch, Aggregate, Pooled, WGS, Compare) up to the same level of refactoring as CRISPRessoCORE: PlotContext construction, data_prep extraction, and full CRISPRessoPro plugin hook ownership of plots and reports. Also includes a small cleanup branch for CRISPRessoCORE itself.

**Depends on:** Completed CRISPRessoCORE plot plugin work (see `DATA_PREP_EXTRACTION.md`, `PLOT_CONTEXT_HIERARCHY.md`, and CRISPRessoPro's `PLOT_PLUGIN_ARCHITECTURE.md` / `FULL_DATA_FUNC_COVERAGE.md`).

## Branch Structure

6 branches, each independent (no cross-branch dependencies). Ordered by complexity:

| # | Branch | Scope | Complexity |
|---|--------|-------|------------|
| 0 | `core-plot-cleanup` | 3 polish items in CRISPRessoCORE | Trivial |
| 1 | `pooled-plot-plugin` | 2 plots | Trivial |
| 2 | `wgs-plot-plugin` | 2 plots (reuses Pooled's prep functions) | Trivial |
| 3 | `compare-plot-plugin` | ~10 plot calls across 3 plot types | Medium |
| 4 | `batch-plot-plugin` | ~6 plot types, branching, Pro-only plots | High |
| 5 | `aggregate-plot-plugin` | Same as Batch + pagination + reads summary | High |

Batch → Aggregate is the only natural ordering (code similarity makes Aggregate easier after Batch). All other branches are independent.

## Architecture

### Per-Branch Work Pattern

For each mode, in order:

1. **Write data_prep functions** — Extract plot input computation from the CORE module into `data_prep.py` functions that take the mode's context
2. **Refactor the `elif` branch** — CORE module calls `prep_*()` functions instead of inline dict construction
3. **Construct context in CORE** — Build the `XxxPlotContext` before the plotting section (classes already defined in `plot_context.py`)
4. **Add the Pro `if` branch** — `if C2PRO_INSTALLED: hook(ctx)`
5. **Pro side** — hooks, PlotDef registrations, data_funcs (which call the data_prep functions), report hook
6. **Tests** — data_prep unit tests, context construction tests, Pro data_func tests, integration

### Hook Pattern

Each mode gets two hook functions in Pro's `hooks.py`:

```python
on_pooled_plots_complete(ctx: PooledPlotContext, logger)
on_wgs_plots_complete(ctx: WGSPlotContext, logger)
on_compare_plots_complete(ctx: ComparePlotContext, logger)
on_batch_plots_complete(ctx: BatchPlotContext, logger)
on_aggregate_plots_complete(ctx: AggregatePlotContext, logger)
```

Plus a corresponding `make_xxx_report()` for each. These follow the exact same lifecycle as Core's existing hooks:

1. Build a mode-specific registry (e.g., `register_batch_plots(registry, plot_module)`)
2. Iterate scopes, call data_funcs, generate plots
3. Write metadata to `crispresso2_info` using `_record_plugin_plot`
4. For report: assemble figures from registry + metadata, render the shared template

### CORE Module Pattern

Each `XxxCORE.py` gets this at the plotting boundary:

```python
ctx = XxxPlotContext(args=args, run_data=crispresso2_info, ...)

if C2PRO_INSTALLED:
    try:
        from CRISPRessoPro import hooks as pro_hooks
        pro_hooks.on_xxx_plots_complete(ctx, logger)
    except Exception as e:
        if args.halt_on_plot_fail:
            raise
        logger.warning(f"CRISPRessoPro plugin hook failed: {e}")
elif not args.suppress_plots:
    # existing matplotlib plotting code, refactored to use prep functions
    ...
```

Report generation follows the same pattern:

```python
if not args.suppress_report:
    if C2PRO_INSTALLED:
        from CRISPRessoPro import hooks as pro_hooks
        pro_hooks.make_xxx_report(crispresso2_info, report_name, OUTPUT_DIRECTORY, _ROOT, logger)
    else:
        CRISPRessoReport.make_xxx_report_from_folder(report_name, crispresso2_info, OUTPUT_DIRECTORY, _ROOT, logger)
```

### Registry & Scopes Per Mode

Each mode gets its own PlotDef registration function. Scope names (`global`, `amplicon`, `sgRNA`) are reused but have mode-specific semantics:

- **Core:** `amplicon` = per reference amplicon in one sample; `sgRNA` = per guide within a reference
- **Batch/Aggregate:** `amplicon` = per consensus amplicon across all batches/samples; `sgRNA` = per consensus guide within a consensus amplicon
- **Pooled/WGS:** `global` only (2 summary plots)
- **Compare:** `amplicon` = per amplicon present in both samples

The iteration logic is defined per mode in each hook function — the hook knows what `amplicon` and `sgRNA` mean for its mode.

### Report Template

The existing Pro `report.html` is the single entry point for all modes (option A). Mode-specific differences are handled by:

- Passing a `mode` variable to the template (e.g., `'batch'`, `'pooled'`)
- Mode-specific partials (e.g., `partials/batch_runs_table.html`) rendered conditionally
- The figure iteration loop is identical across all modes

### Metadata Convention

Each Pro hook writes metadata to `crispresso2_info` using the same key conventions that the mode's existing `CRISPRessoReport.make_xxx_report_from_folder` reads, preserving backward compatibility.

### Shared Prep Functions

Several prep functions are reused across modes:

- `prep_reads_total(ctx, prefix)` / `prep_unmod_mod_pcts(ctx, prefix)`: Pooled, WGS, Aggregate
- `prep_batch_nuc_quilt_around_sgRNA` / `prep_batch_nuc_quilt`: Batch, Aggregate
- `prep_batch_allele_modification_heatmap` / `prep_batch_allele_modification_line`: Batch, Aggregate

Functions shared between Batch and Aggregate use explicit union types: `BatchPlotContext | AggregatePlotContext`.

---

## Branch 0: `core-plot-cleanup`

### Current Gaps

Three polish items from the CRISPRessoCORE data_prep extraction:

1. **Plot 10i** (upset plot, ~line 5215) constructs its input dict inline rather than using a prep function return
2. **Redundant local variable unpacking** at the top of the per-ref loop (~line 4770): `ref_seq`, `include_idxs_list`, `tot_aln_reads`, etc. used for CSV writes
3. **Plot 10i calls `CRISPRessoPlot.plot_combination_upset` directly** rather than going through the `plot()` multiprocessing wrapper

### Changes

**`data_prep.py`:** No new functions needed — `prep_base_edit_upset` already exists. The issue is that CORE doesn't use its return value to build the plot input.

**`CRISPRessoCORE.py`:**
- Refactor plot 10i to use `prep_base_edit_upset`'s return value for the plot input dict (or add a `prep_base_edit_upset_plot` that wraps the result)
- Remove local variable unpacking block; CSV write code accesses `refs[ref_name]` directly
- Route plot 10i through the `plot()` wrapper

---

## Branch 1: `pooled-plot-plugin`

### Current State

Two plots in CRISPRessoPooledCORE.py (~lines 1613–1635), called synchronously:

```python
if not args.suppress_plots:
    CRISPRessoPlot.plot_reads_total(df_summary_quantification=..., fig_filename_root=..., save_png=..., cutoff=...)
    CRISPRessoPlot.plot_unmod_mod_pcts(df_summary_quantification=..., fig_filename_root=..., save_png=..., cutoff=...)
```

Metadata: `crispresso2_info['results']['general_plots']['summary_plot_root']`, `summary_plot_titles`, `summary_plot_labels`, `summary_plot_datas`.

### Data Prep Functions

Two new functions in `data_prep.py`, typed to accept `PlotContext` subclasses with `df_summary_quantification`:

**`prep_reads_total(ctx, prefix: str)`** — Returns kwargs dict for `plot_reads_total`:
- `df_summary_quantification`, `fig_filename_root` (using prefix), `save_png`, `cutoff` (from `args.min_reads_to_use_region`)

**`prep_unmod_mod_pcts(ctx, prefix: str)`** — Same shape, different filename root.

Both are ~5 lines. The `prefix` parameter enables reuse across Pooled (`'CRISPRessoPooled'`), WGS (`'CRISPRessoWGS'`), and Aggregate (`'CRISPRessoAggregate'`).

### CORE Module Flow

```python
ctx = PooledPlotContext(
    args=args, run_data=crispresso2_info, output_directory=OUTPUT_DIRECTORY,
    save_png=save_png, _jp=_jp, custom_config=custom_config,
    df_summary_quantification=df_summary_quantification,
)

if C2PRO_INSTALLED:
    try:
        from CRISPRessoPro import hooks as pro_hooks
        pro_hooks.on_pooled_plots_complete(ctx, logger)
    except Exception as e:
        if args.halt_on_plot_fail:
            raise
        logger.warning(f"CRISPRessoPro plugin hook failed: {e}")
elif not args.suppress_plots:
    reads_input = prep_reads_total(ctx, prefix='CRISPRessoPooled')
    CRISPRessoPlot.plot_reads_total(**reads_input)
    # ... metadata registration (same keys as today)

    pcts_input = prep_unmod_mod_pcts(ctx, prefix='CRISPRessoPooled')
    CRISPRessoPlot.plot_unmod_mod_pcts(**pcts_input)
    # ... metadata registration
```

Report hook added at the report generation site.

### Pro Side

- **`hooks.py`:** `on_pooled_plots_complete(ctx, logger)` — builds registry, generates 2 plots via data_funcs, records metadata
- **`hooks.py`:** `make_pooled_report(run_data, report_file, output_dir, _root, logger)` — renders data-driven template with `mode='pooled'`
- **PlotDefs:** 2 registrations (`pooled_reads_summary`, `pooled_modification_summary`), both scope `global`
- **data_funcs:** 2 functions calling `prep_reads_total(ctx, 'CRISPRessoPooled')` and `prep_unmod_mod_pcts(ctx, 'CRISPRessoPooled')`

---

## Branch 2: `wgs-plot-plugin`

Identical pattern to Pooled. Same 2 plots, same prep functions (reused via `prefix='CRISPRessoWGS'`), `WGSPlotContext` instead of `PooledPlotContext`.

Only difference: the metadata key for reads summary is `reads_summary_plot` (WGS) vs `summary_plot_root` (Pooled) — existing inconsistency preserved for backward compat.

- **`hooks.py`:** `on_wgs_plots_complete`, `make_wgs_report`
- **PlotDefs:** 2 registrations
- **data_funcs:** 2 functions

---

## Branch 3: `compare-plot-plugin`

### Current State

Three plot types in CRISPRessoCompareCORE.py (~lines 200–470), all inside `for amplicon_name in amplicon_names_in_both:`:

| Plot | Per amplicon | Details |
|------|-------------|---------|
| Editing comparison barchart | 1 | Computes N_TOTAL/UNMODIFIED/MODIFIED per sample |
| Quantification positions | 4 (one per mod type) | Fisher exact test + Bonferroni per position; writes CSV |
| Allele table comparison | 2 per allele file pair (top/bottom enriched) | Merge DataFrames, compute LFC, sort both directions |

### ComparePlotContext Additions

New mutable scope-like fields, all `Optional` with `None` defaults. CORE populates these per-amplicon iteration:

```python
# Per-amplicon loaded data (set by CORE each iteration)
profile_1: Optional[np.ndarray] = None
profile_2: Optional[np.ndarray] = None
mod_freqs_1: Optional[dict] = None
mod_freqs_2: Optional[dict] = None
consensus_sequence: Optional[str] = None
quant_windows_1: Optional[np.ndarray] = None
quant_windows_2: Optional[np.ndarray] = None
cut_points: Optional[list] = None
sgRNA_intervals: Optional[list] = None

# Per-mod-type scope
mod_type: Optional[str] = None

# Allele comparison data
allele_pairs: Optional[list] = None  # list of (df1, df2, file_root, is_base_edit)
```

### Data Prep Functions

**`prep_compare_editing_barchart(ctx: ComparePlotContext)`**
- Reads `ctx.amplicon_info_1[ctx.amplicon_name]` for N_TOTAL/UNMODIFIED/MODIFIED from each sample
- Returns kwargs dict for `plot_quantification_comparison_barchart`

**`prep_compare_modification_positions(ctx: ComparePlotContext)`**
- Uses `ctx.mod_freqs_1`, `ctx.mod_freqs_2`, `ctx.mod_type`
- Computes Fisher exact test per position + Bonferroni correction
- Returns dict with:
  - `plot_kwargs` for `plot_quantification_positions`
  - `mod_df` for CSV write
  - `sig_count` and `sig_count_quant_window` for summary table

**`prep_compare_allele_table(ctx: ComparePlotContext, pair_index: int)`**
- Uses `ctx.allele_pairs[pair_index]` to get (df1, df2, file_root, is_base_edit)
- Merges DataFrames, computes LFC
- Returns dict with:
  - `merged_df` for CSV write
  - `ref_seq_around_cut`
  - `plot_top_kwargs` (sorted for enrichment in sample_1)
  - `plot_bottom_kwargs` (sorted for enrichment in sample_2)
  - `is_base_edit`, `file_root` for metadata

### CORE Module Flow

```python
ctx = ComparePlotContext(
    args=args, run_data=crispresso2_info, ...,
    amplicon_names=amplicon_names_in_both,
    sample_1_name=sample_1_name, sample_2_name=sample_2_name,
    run_info_1=run_info_1, run_info_2=run_info_2,
    amplicon_info_1=amplicon_info_1, amplicon_info_2=amplicon_info_2,
)

if C2PRO_INSTALLED:
    pro_hooks.on_compare_plots_complete(ctx, logger)
elif:
    for amplicon_name in amplicon_names_in_both:
        ctx.amplicon_name = amplicon_name
        ctx.profile_1 = parse_profile(...)
        ctx.profile_2 = parse_profile(...)
        ctx.mod_freqs_1 = parse_count_file(...)  # etc.

        # Plot 1
        barchart_data = prep_compare_editing_barchart(ctx)
        CRISPRessoPlot.plot_quantification_comparison_barchart(**barchart_data)

        # Plot 2 (×4)
        for mod_type in ['Insertions', 'Deletions', 'Substitutions', 'All_modifications']:
            ctx.mod_type = mod_type
            positions_data = prep_compare_modification_positions(ctx)
            positions_data['mod_df'].to_csv(...)
            CRISPRessoPlot.plot_quantification_positions(**positions_data['plot_kwargs'])
            sig_counts[amplicon_name][mod_type] = positions_data['sig_count']

        # Plot 3 (per allele file pair)
        ctx.allele_pairs = load_matching_allele_pairs(...)
        for pair_idx in range(len(ctx.allele_pairs)):
            allele_data = prep_compare_allele_table(ctx, pair_idx)
            allele_data['merged_df'].to_csv(...)
            CRISPRessoPlot.plot_alleles_table_compare(**allele_data['plot_top_kwargs'])
            CRISPRessoPlot.plot_alleles_table_compare(**allele_data['plot_bottom_kwargs'])
```

### Pro Side

- **`hooks.py`:** `on_compare_plots_complete` — iterates amplicons, loads data onto ctx, generates all plots via data_funcs
- **PlotDefs:** 3 registrations, all scope `amplicon`
- **data_funcs:** 3 functions. `modification_positions` returns a list (one per mod_type). `allele_table` returns a list (one per pair × 2 sort orders).

---

## Branch 4: `batch-plot-plugin`

### Current State

~350 lines of plotting with significant branching. Data aggregation loop (~lines 430–570) builds per-amplicon summary DataFrames. Plotting section (~lines 570–860) generates plots.

| Plot | When | Iteration |
|------|------|-----------|
| Nucleotide quilt (per-sgRNA) | `guides_all_same` | per amplicon × per sgRNA |
| Nucleotide quilt (whole-region) | always | per amplicon |
| Conversion map (per-sgRNA) | `guides_all_same` and `base_editor_output` | per amplicon × per sgRNA |
| Conversion map (whole-region) | `base_editor_output` | per amplicon |
| Allele modification heatmap | **Pro-only** (`C2PRO_INSTALLED and not use_matplotlib`) | per amplicon × per mod_type (3) |
| Allele modification line | **Pro-only** (same gate) | per amplicon × per mod_type (3) |

### BatchPlotContext Additions

```python
mod_type: Optional[str] = None  # scope field for heatmap/line iteration
```

All other fields already defined (amplicon_names, summary DataFrames, consensus guide data, guides_all_same — all as dicts keyed by amplicon name).

### CORE Restructuring

The current code interleaves data aggregation and plotting in a single per-amplicon loop. This must be separated into two passes:

1. **Aggregation pass:** Iterates amplicons × batch folders, builds all summary DataFrames, stores in dicts keyed by amplicon name
2. **Context construction:** One `BatchPlotContext` with all amplicons' data
3. **Plotting pass:** `if C2PRO_INSTALLED` branch or `elif` with prep function calls

This is the same transformation Core did when it moved all analysis before `CorePlotContext` construction.

### Data Prep Functions

Six functions, four of which are typed `BatchPlotContext | AggregatePlotContext` for reuse:

**`prep_batch_nuc_quilt_around_sgRNA(ctx: BatchPlotContext | AggregatePlotContext)`**
- Uses `ctx.amplicon_name`, `ctx.sgRNA_ind`
- Computes sub-sgRNA coordinate mapping: `sub_sgRNA_intervals`, `sub_include_idxs` by mapping full amplicon coordinates to the sgRNA plot window
- Slices `nucleotide_percentage_summary_df` and `modification_percentage_summary_df` to sgRNA columns
- Returns kwargs dict for `plot_nucleotide_quilt`

**`prep_batch_nuc_quilt(ctx: BatchPlotContext | AggregatePlotContext)`**
- Uses `ctx.amplicon_name`
- Packages full DataFrames with sgRNA_intervals and include_idxs
- Returns kwargs dict for `plot_nucleotide_quilt`

**`prep_batch_conversion_map_around_sgRNA(ctx: BatchPlotContext)`**
- Same sub-indexing as sgRNA quilt, returns kwargs for `plot_conversion_map`
- Includes `conversion_nuc_from`, `conversion_nuc_to` from args

**`prep_batch_conversion_map(ctx: BatchPlotContext)`**
- Whole-region conversion map kwargs

**`prep_batch_allele_modification_heatmap(ctx: BatchPlotContext | AggregatePlotContext)`**
- Uses `ctx.mod_type` to filter `modification_frequency_summary_df`
- Formats DataFrame: renames columns to `'{base} ({position})'`, renames rows to `'{batch} ({index})'`
- Computes sgRNA_intervals list for all rows
- Returns kwargs for `plot_allele_modification_heatmap`

**`prep_batch_allele_modification_line(ctx: BatchPlotContext | AggregatePlotContext)`**
- Same data shaping as heatmap, returns kwargs for `plot_allele_modification_line`

### CORE Module Flow

```python
# Phase 1: Aggregation (builds dicts of DataFrames per amplicon)
nuc_freq_dfs, nuc_pct_dfs = {}, {}
mod_freq_dfs, mod_pct_dfs = {}, {}
consensus_guides_by_amp, consensus_include_idxs_by_amp = {}, {}
# ... iterate amplicons × batch folders, populate dicts ...

# Phase 2: Context construction
ctx = BatchPlotContext(
    args=args, run_data=crispresso2_info, ...,
    amplicon_names=amplicon_name_list,
    nucleotide_frequency_summary_dfs=nuc_freq_dfs,
    nucleotide_percentage_summary_dfs=nuc_pct_dfs,
    # ... etc
)

# Phase 3: Plotting
if C2PRO_INSTALLED:
    pro_hooks.on_batch_plots_complete(ctx, logger)
elif not args.suppress_plots:
    for amplicon_name in ctx.amplicon_names:
        ctx.amplicon_name = amplicon_name
        if ctx.guides_all_same[amplicon_name]:
            for sgRNA_ind, sgRNA in enumerate(ctx.consensus_guides[amplicon_name]):
                ctx.sgRNA_ind = sgRNA_ind
                quilt_input = prep_batch_nuc_quilt_around_sgRNA(ctx)
                plot(CRISPRessoPlot.plot_nucleotide_quilt, quilt_input)
                # ... metadata, CSV writes
                if args.base_editor_output:
                    conv_input = prep_batch_conversion_map_around_sgRNA(ctx)
                    plot(CRISPRessoPlot.plot_conversion_map, conv_input)

        quilt_input = prep_batch_nuc_quilt(ctx)
        plot(CRISPRessoPlot.plot_nucleotide_quilt, quilt_input)
        if args.base_editor_output:
            conv_input = prep_batch_conversion_map(ctx)
            plot(CRISPRessoPlot.plot_conversion_map, conv_input)
```

Note: Pro-only heatmap/line plots **disappear from the `elif` branch**. They only existed when `C2PRO_INSTALLED` was true, so they move exclusively into Pro's hook.

### Pro Side

- **`hooks.py`:** `on_batch_plots_complete(ctx, logger)` — iterates amplicons × sgRNA × mod_types, generates all 6 plot types
- **PlotDefs:** 6 registrations. Quilts and conversion maps use scope `amplicon` / `sgRNA` with appropriate `applicable` lambdas. Heatmap/line use scope `amplicon`.
- **data_funcs:** 6 functions calling the prep functions

---

## Branch 5: `aggregate-plot-plugin`

### Current State

Structurally similar to Batch with differences:

- **Same as Batch:** Nucleotide quilt (per-sgRNA + whole-region), allele modification heatmap/line (Pro-only)
- **No conversion map plots** — Aggregate doesn't have them
- **Adds `plot_reads_total` and `plot_unmod_mod_pcts`** — summary plots Batch lacks
- **Pagination** — when sample count exceeds `max_samples_per_summary_plot`, whole-region quilt splits across pages with suffix `_1`, `_2`, etc.

### AggregatePlotContext Additions

```python
mod_type: Optional[str] = None  # scope field for heatmap/line iteration
```

All other fields already defined (same per-amplicon dict structure as Batch, plus `df_summary_quantification` and `sample_count`).

### Shared Prep Functions

Reuses from other branches:

- `prep_batch_nuc_quilt_around_sgRNA(ctx: BatchPlotContext | AggregatePlotContext)` — from Batch
- `prep_batch_nuc_quilt(ctx: BatchPlotContext | AggregatePlotContext)` — from Batch
- `prep_batch_allele_modification_heatmap(ctx: BatchPlotContext | AggregatePlotContext)` — from Batch
- `prep_batch_allele_modification_line(ctx: BatchPlotContext | AggregatePlotContext)` — from Batch
- `prep_reads_total(ctx, prefix='CRISPRessoAggregate')` — from Pooled
- `prep_unmod_mod_pcts(ctx, prefix='CRISPRessoAggregate')` — from Pooled

No new Aggregate-specific prep functions needed.

### Pagination

Pagination is a display concern, not data preparation. The prep function (`prep_batch_nuc_quilt`) returns the full DataFrame in its kwargs. CORE's pagination loop slices the returned DataFrame before dispatching each page:

```python
quilt_input = prep_batch_nuc_quilt(ctx)
full_df = quilt_input['nuc_pct_df']
for start in range(0, len(full_df), max_rows):
    page_input = dict(quilt_input)
    page_input['nuc_pct_df'] = full_df.iloc[start:start + max_rows]
    page_input['fig_filename_root'] = base_root + suffix
    plot(CRISPRessoPlot.plot_nucleotide_quilt, page_input)
```

Pro's hook can paginate the same way or handle large datasets natively (plotly).

### CORE Restructuring

Same two-pass separation as Batch: aggregation → context construction → plotting.

### Pro Side

- **`hooks.py`:** `on_aggregate_plots_complete` — generates quilts, heatmap/line, reads_total, unmod_mod_pcts
- **`hooks.py`:** `make_aggregate_report` — renders data-driven template with `mode='aggregate'`
- **PlotDefs:** ~6 registrations (quilts + heatmap/line + 2 summary plots)
- **data_funcs:** ~6 functions calling shared prep functions

---

## Summary

### Files Changed Per Branch

**Branch 0 (`core-plot-cleanup`):**

| File | Change |
|------|--------|
| `CRISPResso2/CRISPRessoCORE.py` | Refactor 10i, remove variable unpacking, route 10i through plot() |

**Branch 1 (`pooled-plot-plugin`):**

| File | Change |
|------|--------|
| `CRISPResso2/plots/data_prep.py` | Add `prep_reads_total`, `prep_unmod_mod_pcts` |
| `CRISPResso2/CRISPRessoPooledCORE.py` | Context construction, `if C2PRO_INSTALLED` branch, prep calls, report hook |
| `CRISPRessoPro/hooks.py` | Add `on_pooled_plots_complete`, `make_pooled_report` |
| `CRISPRessoPro/plots/builtin_plots.py` (or new file) | Register 2 PlotDefs |
| `CRISPRessoPro/plots/data_funcs.py` (or new file) | Add 2 data_funcs |
| `CRISPRessoPro/templates/partials/` | Add pooled sub-run table partial |
| `tests/unit_tests/test_data_prep.py` | Tests for 2 prep functions |
| `CRISPRessoPro/tests/` | Tests for 2 data_funcs, hook |

**Branch 2 (`wgs-plot-plugin`):**

| File | Change |
|------|--------|
| `CRISPResso2/CRISPRessoWGSCORE.py` | Context construction, `if C2PRO_INSTALLED` branch, prep calls, report hook |
| `CRISPRessoPro/hooks.py` | Add `on_wgs_plots_complete`, `make_wgs_report` |
| `CRISPRessoPro/plots/builtin_plots.py` | Register 2 PlotDefs |
| `CRISPRessoPro/plots/data_funcs.py` | Add 2 data_funcs |
| `CRISPRessoPro/templates/partials/` | Add WGS sub-run table partial |
| `CRISPRessoPro/tests/` | Tests for 2 data_funcs, hook |

**Branch 3 (`compare-plot-plugin`):**

| File | Change |
|------|--------|
| `CRISPResso2/plots/plot_context.py` | Add ~11 scope fields to `ComparePlotContext` |
| `CRISPResso2/plots/data_prep.py` | Add 3 prep functions |
| `CRISPResso2/CRISPRessoCompareCORE.py` | Context construction, `if C2PRO_INSTALLED` branch, prep calls, report hook |
| `CRISPRessoPro/hooks.py` | Add `on_compare_plots_complete`, `make_compare_report` |
| `CRISPRessoPro/plots/builtin_plots.py` | Register 3 PlotDefs |
| `CRISPRessoPro/plots/data_funcs.py` | Add 3 data_funcs |
| `CRISPRessoPro/templates/partials/` | Add compare sub-run table partial |
| `tests/unit_tests/test_data_prep.py` | Tests for 3 prep functions |
| `tests/unit_tests/test_plot_context.py` | Tests for new scope fields |
| `CRISPRessoPro/tests/` | Tests for 3 data_funcs, hook |

**Branch 4 (`batch-plot-plugin`):**

| File | Change |
|------|--------|
| `CRISPResso2/plots/plot_context.py` | Add `mod_type` scope field to `BatchPlotContext` |
| `CRISPResso2/plots/data_prep.py` | Add 6 prep functions |
| `CRISPResso2/CRISPRessoBatchCORE.py` | Two-pass restructure, context construction, `if C2PRO_INSTALLED` branch, prep calls, report hook |
| `CRISPRessoPro/hooks.py` | Add `on_batch_plots_complete`, `make_batch_report` |
| `CRISPRessoPro/plots/builtin_plots.py` | Register 6 PlotDefs |
| `CRISPRessoPro/plots/data_funcs.py` | Add 6 data_funcs |
| `CRISPRessoPro/templates/partials/` | Add batch sub-run table partial |
| `tests/unit_tests/test_data_prep.py` | Tests for 6 prep functions |
| `tests/unit_tests/test_plot_context.py` | Tests for `mod_type` field |
| `CRISPRessoPro/tests/` | Tests for 6 data_funcs, hook |

**Branch 5 (`aggregate-plot-plugin`):**

| File | Change |
|------|--------|
| `CRISPResso2/plots/plot_context.py` | Add `mod_type` scope field to `AggregatePlotContext` |
| `CRISPResso2/CRISPRessoAggregateCORE.py` | Two-pass restructure, context construction, `if C2PRO_INSTALLED` branch, prep calls, pagination, report hook |
| `CRISPRessoPro/hooks.py` | Add `on_aggregate_plots_complete`, `make_aggregate_report` |
| `CRISPRessoPro/plots/builtin_plots.py` | Register ~6 PlotDefs |
| `CRISPRessoPro/plots/data_funcs.py` | Add ~6 data_funcs |
| `CRISPRessoPro/templates/partials/` | Add aggregate sub-run table partial |
| `tests/unit_tests/test_plot_context.py` | Tests for `mod_type` field |
| `CRISPRessoPro/tests/` | Tests for ~6 data_funcs, hook |

### Totals

| Metric | Count |
|--------|-------|
| New unique prep functions | ~11 |
| New PlotDefs | ~19 |
| New data_funcs | ~19 |
| New hook functions | 10 (5 plot hooks + 5 report hooks) |
| New report partials | 5 |

---

## Acceptance Criteria

### Shared Infrastructure

1. The system shall define mode-specific scope values for Batch and Aggregate PlotDef registrations, where `amplicon` means "per consensus amplicon across samples" and `sgRNA` means "per consensus guide within an amplicon."

2. When a Pro hook function for any mode generates plots, the system shall write metadata to `crispresso2_info` using the same key conventions that the mode's existing `CRISPRessoReport.make_xxx_report_from_folder` function reads, preserving backward compatibility.

3. When Pro's data-driven report template renders a non-Core mode, the system shall accept a `mode` variable and load mode-specific partials (e.g., sub-run tables) conditionally, reusing the shared figure iteration logic.

### Branch 0: `core-plot-cleanup`

4. When plot 10i (base edit upset) is generated in the built-in plot branch, the system shall produce its input data by calling a prep function rather than constructing the input dict inline.

5. The system shall not contain redundant local variable unpacking from `refs[ref_name]` at the top of the per-reference plotting loop — the remaining CSV write code shall access reference data through `refs[ref_name]` directly or through the context.

6. When plot 10i is generated, the system shall dispatch it through the `plot()` multiprocessing wrapper, consistent with all other plots.

7. When `make all test diff-plots` is run in `../CRISPResso2_tests`, the system shall produce identical output to the expected results.

### Branch 1: `pooled-plot-plugin`

8. The system shall provide `prep_reads_total(ctx, prefix)` in `data_prep.py` that takes a `PlotContext` subclass with a `df_summary_quantification` field and a filename prefix string, and returns the kwargs dict for `plot_reads_total`.

9. The system shall provide `prep_unmod_mod_pcts(ctx, prefix)` in `data_prep.py` following the same pattern as `prep_reads_total`.

10. When `CRISPRessoPooledCORE.main()` runs, the system shall construct a `PooledPlotContext` before the plotting section, after `df_summary_quantification` is built.

11. When CRISPRessoPro is installed and `CRISPRessoPooledCORE.main()` runs, the system shall call `pro_hooks.on_pooled_plots_complete(ctx, logger)` and skip the built-in matplotlib plotting section.

12. When CRISPRessoPro is not installed and plots are not suppressed, the system shall call `prep_reads_total` and `prep_unmod_mod_pcts` to build plot inputs, dispatch the plots, and register metadata — replacing the current inline calls.

13. When `on_pooled_plots_complete` is called, the system shall build a mode-specific registry, generate both plots via data_funcs that call the data_prep functions, and record metadata to `crispresso2_info`.

14. When `make_pooled_report` is called, the system shall render the data-driven report template with `mode='pooled'` and a sub-run table partial showing per-region CRISPResso runs.

15. If `on_pooled_plots_complete` raises an exception and `halt_on_plot_fail` is False, then the system shall log a warning and continue.

16. When the report is generated and CRISPRessoPro is installed, `CRISPRessoPooledCORE` shall call `pro_hooks.make_pooled_report` instead of `CRISPRessoReport.make_pooled_report_from_folder`.

17. When `pytest tests/unit_tests/` is run, the system shall include tests for `prep_reads_total` and `prep_unmod_mod_pcts` verifying correct kwargs structure.

18. When `make pooled test` is run, the system shall produce identical output to the expected results.

### Branch 2: `wgs-plot-plugin`

19. When `CRISPRessoWGSCORE.main()` runs, the system shall construct a `WGSPlotContext` before the plotting section.

20. When CRISPRessoPro is installed and `CRISPRessoWGSCORE.main()` runs, the system shall call `pro_hooks.on_wgs_plots_complete(ctx, logger)` and skip the built-in plotting section.

21. When CRISPRessoPro is not installed and plots are not suppressed, the system shall call `prep_reads_total` and `prep_unmod_mod_pcts` (reused from Pooled) with prefix `'CRISPRessoWGS'`.

22. When `make_wgs_report` is called by Pro, the system shall render the data-driven report template with `mode='wgs'`.

23. When `make wgs test` is run, the system shall produce identical output to the expected results.

### Branch 3: `compare-plot-plugin`

24. The `ComparePlotContext` dataclass shall include mutable scope fields for per-amplicon loaded data: `profile_1`, `profile_2`, `mod_freqs_1`, `mod_freqs_2`, `consensus_sequence`, `quant_windows_1`, `quant_windows_2`, `cut_points`, `sgRNA_intervals`, `mod_type`, and `allele_pairs`, all defaulting to `None`.

25. The system shall provide `prep_compare_editing_barchart(ctx: ComparePlotContext)` in `data_prep.py` that returns kwargs for `plot_quantification_comparison_barchart`.

26. The system shall provide `prep_compare_modification_positions(ctx: ComparePlotContext)` that uses `ctx.mod_type` to select the modification, computes Fisher exact tests and Bonferroni correction, and returns a dict with `plot_kwargs` for `plot_quantification_positions`, `mod_df` for CSV write, `sig_count`, and `sig_count_quant_window`.

27. The system shall provide `prep_compare_allele_table(ctx: ComparePlotContext, pair_index: int)` that merges allele DataFrames from both samples, computes LFC, and returns `merged_df`, `ref_seq_around_cut`, `plot_top_kwargs`, and `plot_bottom_kwargs`.

28. When `CRISPRessoCompareCORE.main()` runs, the system shall construct a `ComparePlotContext` after loading both run infos and before the per-amplicon loop.

29. When CRISPRessoPro is installed, the system shall call `pro_hooks.on_compare_plots_complete(ctx, logger)` and skip the built-in plotting section.

30. When CRISPRessoPro is not installed, the system shall load per-amplicon data onto the context's scope fields and call prep functions to build plot inputs, replacing inline dict construction and computation.

31. When `on_compare_plots_complete` is called, the system shall iterate amplicons, load per-amplicon data onto context scope fields, and generate all three plot types via data_funcs.

32. When `make_compare_report` is called, the system shall render the data-driven report template with `mode='compare'`.

33. When `pytest tests/unit_tests/` is run, the system shall include tests for all three Compare prep functions, including Fisher test correctness and LFC computation.

34. When `make compare test` is run, the system shall produce identical output to the expected results.

### Branch 4: `batch-plot-plugin`

35. The `BatchPlotContext` dataclass shall include a `mod_type: Optional[str]` scope field defaulting to `None`.

36. The system shall provide `prep_batch_nuc_quilt_around_sgRNA(ctx: BatchPlotContext | AggregatePlotContext)` that computes sub-sgRNA coordinate mapping and returns kwargs for `plot_nucleotide_quilt` with the sliced summary DataFrames.

37. The system shall provide `prep_batch_nuc_quilt(ctx: BatchPlotContext | AggregatePlotContext)` that returns kwargs for `plot_nucleotide_quilt` with the full-region summary DataFrames.

38. The system shall provide `prep_batch_conversion_map_around_sgRNA(ctx: BatchPlotContext)` and `prep_batch_conversion_map(ctx: BatchPlotContext)` that return kwargs for `plot_conversion_map`.

39. The system shall provide `prep_batch_allele_modification_heatmap(ctx: BatchPlotContext | AggregatePlotContext)` and `prep_batch_allele_modification_line(ctx: BatchPlotContext | AggregatePlotContext)` that use `ctx.mod_type`, format the modification DataFrame, and return kwargs for their respective plot functions.

40. When `CRISPRessoBatchCORE.main()` runs, the system shall separate the per-amplicon data aggregation phase from the plotting phase, storing all per-amplicon summary DataFrames in dicts keyed by amplicon name before constructing the context.

41. When `CRISPRessoBatchCORE.main()` runs, the system shall construct a `BatchPlotContext` after data aggregation completes and before the plotting section begins.

42. When CRISPRessoPro is installed, the system shall call `pro_hooks.on_batch_plots_complete(ctx, logger)` and skip the entire built-in plotting section, including the formerly Pro-only heatmap and line plots.

43. When CRISPRessoPro is not installed and plots are not suppressed, the system shall call prep functions to build plot inputs for all nucleotide quilt and conversion map plots, replacing inline dict construction.

44. When CRISPRessoPro is not installed, the built-in plotting section shall not contain allele modification heatmap or line plot code — those are Pro-only and live exclusively in Pro's hook.

45. When `on_batch_plots_complete` is called, the system shall iterate amplicons, then per-sgRNA within each amplicon (when guides_all_same), then per-mod_type for heatmap/line plots, generating all 6 plot types via data_funcs.

46. When `make_batch_report` is called, the system shall render the data-driven report template with `mode='batch'` and a sub-run table partial showing per-batch CRISPResso runs.

47. When `pytest tests/unit_tests/` is run, the system shall include tests for all 6 Batch prep functions, with the sgRNA coordinate mapping test verifying correct sub-interval computation.

48. When `make batch test` is run, the system shall produce identical output to the expected results.

### Branch 5: `aggregate-plot-plugin`

49. The `AggregatePlotContext` dataclass shall include a `mod_type: Optional[str]` scope field defaulting to `None`.

50. The system shall reuse `prep_batch_nuc_quilt_around_sgRNA`, `prep_batch_nuc_quilt`, `prep_batch_allele_modification_heatmap`, and `prep_batch_allele_modification_line` from the Batch branch by accepting `BatchPlotContext | AggregatePlotContext` as the type hint.

51. The system shall reuse `prep_reads_total` and `prep_unmod_mod_pcts` from the Pooled branch with prefix `'CRISPRessoAggregate'`.

52. When `CRISPRessoAggregateCORE.main()` runs, the system shall separate data aggregation from plotting, construct an `AggregatePlotContext`, and follow the same `if C2PRO_INSTALLED` branch pattern.

53. When CRISPRessoPro is not installed and the number of samples exceeds `max_samples_per_summary_plot`, the system shall paginate the nucleotide quilt by slicing the returned kwargs' DataFrame before dispatching each page to the plot function.

54. When `on_aggregate_plots_complete` is called, the system shall generate all plot types — quilts, heatmap/line (if applicable), reads_total, and unmod_mod_pcts — via data_funcs.

55. When `make_aggregate_report` is called, the system shall render the data-driven report template with `mode='aggregate'`.

56. When `pytest tests/unit_tests/` is run, the system shall include tests verifying that Batch prep functions work when called with an `AggregatePlotContext`.

57. When `make all test` is run, the system shall produce identical output to the expected results.

### Cross-Branch Correctness

58. The system shall not change any plot output, CSV file content, or `crispresso2_info` metadata structure as a result of any branch's refactoring.

59. When CRISPRessoPro is not installed, the system shall produce identical output to the current behavior for all modes.
