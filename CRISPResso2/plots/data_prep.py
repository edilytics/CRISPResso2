"""Extracted data preparation functions for CRISPResso2 plots.

Each function takes a :class:`PlotContext` as its only argument and returns
the kwargs dict that the corresponding CRISPRessoPlot function expects.
These are pure functions — no file I/O, no side effects — with the exception
of ``prep_alleles_around_cut`` and ``prep_base_edit_quilt`` which
intentionally mutate their input DataFrames (see their docstrings).

Only functions with non-trivial computation are extracted here.
Trivial pass-through plots (where the dict is just packaging existing
variables) stay inline in CORE — there's nothing to reuse.

CORE calls these to build plot inputs. CRISPRessoPro can call them
from PlotContext to generate plots independently.

Scope fields on PlotContext
---------------------------
Most functions require ``ctx.ref_name`` to be set (which reference amplicon
is being processed). Functions that operate per-sgRNA additionally require
``ctx.sgRNA_ind``. ``prep_amino_acid_table`` also requires
``ctx.coding_seq_ind``.
"""

from __future__ import annotations

import os
from collections import Counter
import numpy as np
import pandas as pd

from CRISPResso2.plots.plot_context import PlotContext


# =============================================================================
# Private helpers (unchanged)
# =============================================================================


def _to_numeric_ignore_columns(df, ignore_columns):
    """Convert DataFrame columns to numeric, ignoring specified columns."""
    for col in df.columns:
        if col not in ignore_columns:
            df[col] = df[col].apply(pd.to_numeric, errors='raise')
    return df


def plot_title_with_ref_name(title, ref_name, num_refs):
    """Add ref_name suffix to plot title when there are multiple references."""
    if num_refs > 1:
        return title + ": " + ref_name
    return title


def _clip_to_percentile(values, counts, percentile=0.99, total=None):
    """Find the value at which cumulative counts reach the given percentile.

    Scans ``values`` and ``counts`` in order, returning the value where
    the running sum of counts first exceeds ``percentile * total``.
    If the cutoff is never reached, returns the last value.

    Parameters
    ----------
    total : float, optional
        Override the total used for cutoff calculation. If None, uses
        ``counts.sum()``.

    """
    if total is None:
        total = counts.sum()
    if total == 0:
        return values[-1] if len(values) > 0 else 0
    cutoff = percentile * total
    running = 0
    for val, cnt in zip(values, counts):
        running += cnt
        if running > cutoff:
            return val
    return values[-1]


# =============================================================================
# New private helpers for PlotContext extraction
# =============================================================================


def _make_plot_root(ctx: PlotContext, filename: str) -> str:
    """Construct plot root path, or return *filename* if _jp is unavailable."""
    if ctx._jp is not None:
        return ctx._jp(filename)
    return filename


def _ref(ctx: PlotContext) -> dict:
    """Shortcut for the current reference's data dict."""
    return ctx.refs[ctx.ref_name]


def _ref_plot_name(ctx: PlotContext) -> str:
    """Return the current reference's plot-name prefix."""
    return _ref(ctx)['ref_plot_name']


def _sgRNA_label(ctx: PlotContext) -> str:
    """Compute the file-name label for the current sgRNA."""
    from CRISPResso2 import CRISPRessoShared

    ref = _ref(ctx)
    sgRNA = ref['sgRNA_orig_sequences'][ctx.sgRNA_ind]
    sgRNA_name = ref['sgRNA_names'][ctx.sgRNA_ind]
    label = "sgRNA_" + sgRNA
    if sgRNA_name != "":
        label = sgRNA_name
    return CRISPRessoShared.slugify(label)


def _sgRNA_legend(ctx: PlotContext) -> str:
    """Compute the legend string for the current sgRNA."""
    ref = _ref(ctx)
    sgRNA = ref['sgRNA_orig_sequences'][ctx.sgRNA_ind]
    sgRNA_name = ref['sgRNA_names'][ctx.sgRNA_ind]
    legend = "sgRNA " + sgRNA
    if sgRNA_name != "":
        legend = sgRNA_name + " (" + sgRNA + ")"
    return legend


_CLIPPING_NOTE = (
    " Note that histograms are clipped to show 99% of the data."
    " To show all data, run using the parameter '--plot_histogram_outliers'. "
)


def _build_indel_clipped_string(xmin, xmax, raw_xmin, raw_xmax):
    """Build the clipping annotation for plot_3a (indel size distribution).

    Compares the clipped xmin/xmax to the raw data range and returns a
    string describing which extremes were hidden, or ``""`` if nothing
    was clipped.
    """
    parts = []
    if xmax < raw_xmax:
        parts.append(f" (Maximum {int(raw_xmax)} not shown)")
    if xmin > raw_xmin:
        parts.append(f" (Minimum {int(raw_xmin)} not shown)")
    if parts:
        return _CLIPPING_NOTE + "".join(parts)
    return ""


def _build_freq_del_ins_clipped_string(xmax_ins, xmax_del, xmax_mut,
                                       raw_xmax_ins, raw_xmax_del,
                                       raw_xmax_mut):
    """Build the clipping annotation for plot_3b (ins/del/sub histograms).

    Compares the clipped xmax values to the raw maxima and returns a
    string describing which categories were truncated, or ``""`` if
    nothing was clipped.
    """
    parts = []
    if xmax_ins < raw_xmax_ins:
        parts.append(f" (Insertion maximum {int(raw_xmax_ins)} not shown)")
    if xmax_del < raw_xmax_del:
        parts.append(f" (Deletion minimum -{int(raw_xmax_del)} not shown)")
    if xmax_mut < raw_xmax_mut:
        parts.append(f" (Mutation maximum {int(raw_xmax_mut)} not shown)")
    if parts:
        return _CLIPPING_NOTE + "".join(parts)
    return ""


def _build_nuc_freq_df(ctx: PlotContext, ref_name: str | None = None):
    """Build nucleotide frequency and percentage DataFrames for a reference.

    Returns ``(df_nuc_freq, df_nuc_pct)`` where rows are nucleotides
    ``['A', 'C', 'G', 'T', 'N', '-']`` and columns are the reference
    sequence positions.
    """
    if ref_name is None:
        ref_name = ctx.ref_name
    ref_seq = ctx.refs[ref_name]['sequence']
    tot = float(ctx.counts_total[ref_name])

    df_nuc_freq = pd.DataFrame([
        ctx.all_base_count_vectors[ref_name + '_A'],
        ctx.all_base_count_vectors[ref_name + '_C'],
        ctx.all_base_count_vectors[ref_name + '_G'],
        ctx.all_base_count_vectors[ref_name + '_T'],
        ctx.all_base_count_vectors[ref_name + '_N'],
        ctx.all_base_count_vectors[ref_name + '_-'],
    ])
    df_nuc_freq.index = ['A', 'C', 'G', 'T', 'N', '-']
    df_nuc_freq.columns = list(ref_seq)
    df_nuc_pct = df_nuc_freq.divide(tot) if tot > 0 else df_nuc_freq * 0
    return df_nuc_freq, df_nuc_pct


def _compute_half_windows(cut_point, plot_window_size, ref_len):
    """Compute asymmetric left/right window sizes, clamping to amplicon bounds.

    Returns ``(left, right, window_truncated)``.
    """
    window_truncated = False
    if cut_point - plot_window_size + 1 >= 0:
        left = plot_window_size
    else:
        left = cut_point + 1
        window_truncated = True
    if cut_point + plot_window_size < ref_len:
        right = plot_window_size
    else:
        right = ref_len - cut_point - 1
        window_truncated = True
    return left, right, window_truncated


# =============================================================================
# Public prep functions — each takes a single PlotContext
# =============================================================================


def prep_indel_size_distribution(ctx: PlotContext):
    """Prepare kwargs for plot_indel_size_distribution (plot_3a).

    Requires ``ctx.ref_name``.

    Computes xmin/xmax by clipping to the 99th percentile of the
    density distribution (unless ``args.plot_histogram_outliers`` is True),
    then clamping to at least ±15.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)

    hdensity = ref['hdensity']
    hlengths = ref['hlengths']
    center_index = ref['center_index']
    n_this_category = ctx.counts_total[ref_name]
    plot_histogram_outliers = ctx.args.plot_histogram_outliers

    raw_xmin = min(hlengths)
    raw_xmax = max(hlengths)
    xmin = raw_xmin
    xmax = raw_xmax

    if not plot_histogram_outliers:
        xmax = _clip_to_percentile(hlengths, hdensity, 0.99)
        xmin = _clip_to_percentile(hlengths[::-1], hdensity[::-1], 0.99)

    xmin = min(xmin, -15)
    xmax = max(xmax, 15)

    clipped_string = _build_indel_clipped_string(xmin, xmax, raw_xmin, raw_xmax)

    num_refs = len(ctx.ref_names)
    plot_root = _make_plot_root(
        ctx, '3a.' + _ref_plot_name(ctx) + 'Indel_size_distribution',
    )

    return {
        'hdensity': hdensity,
        'hlengths': hlengths,
        'center_index': center_index,
        'n_this_category': n_this_category,
        'xmin': xmin,
        'xmax': xmax,
        'title': plot_title_with_ref_name(
            'Indel size distribution', ref_name, num_refs,
        ),
        'plot_root': plot_root,
        'save_also_png': ctx.save_png,
        'ref_name': ref_name,
        'clipped_string': clipped_string,
    }


def prep_frequency_deletions_insertions(ctx: PlotContext):
    """Prepare kwargs for plot_frequency_deletions_insertions (plot_3b).

    Requires ``ctx.ref_name``.

    Clips each histogram's xmax to the 99th percentile of the overall
    density (unless ``args.plot_histogram_outliers`` is True), clamped to 15.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)

    x_bins_ins = ref['x_bins_ins']
    y_values_ins = ref['y_values_ins']
    x_bins_del = ref['x_bins_del']
    y_values_del = ref['y_values_del']
    x_bins_mut = ref['x_bins_mut']
    y_values_mut = ref['y_values_mut']
    hdensity = ref['hdensity']
    counts_total = ctx.counts_total[ref_name]
    plot_histogram_outliers = ctx.args.plot_histogram_outliers

    raw_xmax_ins = max(x_bins_ins)
    raw_xmax_del = max(x_bins_del)
    raw_xmax_mut = max(x_bins_mut)
    xmax_ins = raw_xmax_ins
    xmax_del = raw_xmax_del
    xmax_mut = raw_xmax_mut

    if not plot_histogram_outliers:
        density_total = hdensity.sum()
        xmax_ins = _clip_to_percentile(x_bins_ins, y_values_ins, 0.99, total=density_total)
        xmax_del = _clip_to_percentile(x_bins_del, y_values_del, 0.99, total=density_total)
        xmax_mut = _clip_to_percentile(x_bins_mut, y_values_mut, 0.99, total=density_total)

    xmax_ins = max(15, xmax_ins)
    xmax_del = max(15, xmax_del)
    xmax_mut = max(15, xmax_mut)

    clipped_string = _build_freq_del_ins_clipped_string(
        xmax_ins, xmax_del, xmax_mut,
        raw_xmax_ins, raw_xmax_del, raw_xmax_mut,
    )

    num_refs = len(ctx.ref_names)
    plot_root = _make_plot_root(
        ctx,
        '3b.' + _ref_plot_name(ctx) + 'Insertion_deletion_substitutions_size_hist',
    )

    return {
        'ref': ref,
        'counts_total': counts_total,
        'plot_root': plot_root,
        'plot_titles': {
            'ins': plot_title_with_ref_name('Insertions', ref_name, num_refs),
            'del': plot_title_with_ref_name('Deletions', ref_name, num_refs),
            'mut': plot_title_with_ref_name('Substitutions', ref_name, num_refs),
        },
        'xmax_ins': xmax_ins,
        'xmax_del': xmax_del,
        'xmax_mut': xmax_mut,
        'save_also_png': ctx.save_png,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'ref_name': ref_name,
        'clipped_string': clipped_string,
    }


def prep_amplicon_modifications(ctx: PlotContext):
    """Prepare kwargs for plot_amplicon_modifications (plot_4a).

    Requires ``ctx.ref_name``.

    Pattern A — packaging with y_max computation and title generation.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    num_refs = len(ctx.ref_names)

    all_indelsub_count_vector = ctx.all_indelsub_count_vectors[ref_name]
    include_idxs_list = ref['include_idxs']
    cut_points = ref['sgRNA_cut_points']
    plot_cut_points = ref['sgRNA_plot_cut_points']
    sgRNA_intervals = ref['sgRNA_intervals']
    ref_len = ref['sequence_length']

    y_max = max(all_indelsub_count_vector) * 1.1

    plot_root = _make_plot_root(
        ctx,
        '4a.' + _ref_plot_name(ctx) + 'Combined_insertion_deletion_substitution_locations',
    )

    return {
        'all_indelsub_count_vectors': all_indelsub_count_vector,
        'include_idxs_list': include_idxs_list,
        'cut_points': cut_points,
        'plot_cut_points': plot_cut_points,
        'sgRNA_intervals': sgRNA_intervals,
        'n_total': ctx.N_TOTAL,
        'n_this_category': ctx.counts_total[ref_name],
        'ref_name': ref_name,
        'num_refs': num_refs,
        'ref_len': ref_len,
        'y_max': y_max,
        'plot_titles': {
            'combined': plot_title_with_ref_name(
                'Combined Insertions/Deletions/Substitutions', ref_name, num_refs,
            ),
            'main': plot_title_with_ref_name(
                'Mutation position distribution', ref_name, num_refs,
            ),
        },
        'plot_root': plot_root,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'save_also_png': ctx.save_png,
    }


def prep_modification_frequency(ctx: PlotContext):
    """Prepare kwargs for plot_modification_frequency (plot_4b).

    Requires ``ctx.ref_name``.

    Pattern A — packaging with title generation. ``y_max`` is computed
    from the combined indelsub count vector (same as plot_4a).
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    num_refs = len(ctx.ref_names)

    y_max = max(ctx.all_indelsub_count_vectors[ref_name]) * 1.1

    plot_root = _make_plot_root(
        ctx,
        '4b.' + _ref_plot_name(ctx) + 'Insertion_deletion_substitution_locations',
    )

    return {
        'include_idxs_list': ref['include_idxs'],
        'all_insertion_count_vectors': ctx.all_insertion_count_vectors[ref_name],
        'all_deletion_count_vectors': ctx.all_deletion_count_vectors[ref_name],
        'all_substitution_count_vectors': ctx.all_substitution_count_vectors[ref_name],
        'sgRNA_intervals': ref['sgRNA_intervals'],
        'ref_len': ref['sequence_length'],
        'ref_name': ref_name,
        'num_refs': num_refs,
        'n_total': ctx.N_TOTAL,
        'n_this_category': ctx.counts_total[ref_name],
        'cut_points': ref['sgRNA_cut_points'],
        'plot_cut_points': ref['sgRNA_plot_cut_points'],
        'y_max': y_max,
        'plot_title': plot_title_with_ref_name(
            'Mutation position distribution', ref_name, num_refs,
        ),
        'plot_root': plot_root,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'save_also_png': ctx.save_png,
    }


def prep_dsODN_piechart(ctx: PlotContext):
    """Prepare kwargs for plot_class_dsODN_piechart (plot_1d).

    No scope fields required.

    Computes labels and sizes from df_alleles 'contains dsODN' column.
    """
    df_alleles = ctx.df_alleles
    N_TOTAL = ctx.N_TOTAL
    plot_root = _make_plot_root(ctx, '1d.Detection_of_dsODN')

    n_contain = df_alleles[df_alleles['contains dsODN'] == True]['#Reads'].sum()
    n_not_contain = df_alleles[df_alleles['contains dsODN'] == False]['#Reads'].sum()
    labels = [
        'Contains dsODN\n(' + str(n_contain) + ' reads)',
        'No dsODN\n(' + str(n_not_contain) + ' reads)',
    ]
    sizes = [
        100 * n_contain / float(N_TOTAL),
        100 * n_not_contain / float(N_TOTAL),
    ]
    return {
        'sizes': sizes,
        'labels': labels,
        'plot_root': plot_root,
        'save_also_png': ctx.save_png,
    }


def prep_nucleotide_quilt(ctx: PlotContext):
    """Prepare kwargs for plot_nucleotide_quilt (plot_2a).

    Requires ``ctx.ref_name``.

    Builds ``nuc_pct_df`` and ``mod_pct_df`` from per-position count vectors.

    The returned dict is passed directly to ``plot_nucleotide_quilt``.
    CORE should also read ``nuc_pct_df`` and ``mod_pct_df`` back from the
    returned dict to pass as inputs to ``prep_nucleotide_quilt_around_sgRNA``.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    ref_seq = ref['sequence']
    tot = float(ctx.counts_total[ref_name])

    _df_nuc_freq_all, df_nuc_pct_all = _build_nuc_freq_df(ctx)

    counts_total = ctx.counts_total[ref_name]

    mod_pcts = []
    mod_pcts.append(np.concatenate((['Insertions'], np.array(ctx.all_insertion_count_vectors[ref_name]).astype(float) / tot)))
    mod_pcts.append(np.concatenate((['Insertions_Left'], np.array(ctx.all_insertion_left_count_vectors[ref_name]).astype(float) / tot)))
    mod_pcts.append(np.concatenate((['Deletions'], np.array(ctx.all_deletion_count_vectors[ref_name]).astype(float) / tot)))
    mod_pcts.append(np.concatenate((['Substitutions'], np.array(ctx.all_substitution_count_vectors[ref_name]).astype(float) / tot)))
    mod_pcts.append(np.concatenate((['All_modifications'], np.array(ctx.all_indelsub_count_vectors[ref_name]).astype(float) / tot)))
    mod_pcts.append(np.concatenate((['Total'], [counts_total] * len(ref_seq))))
    colnames = ['Modification'] + list(ref_seq)
    modification_percentage_summary_df = _to_numeric_ignore_columns(
        pd.DataFrame(mod_pcts, columns=colnames), {'Modification'},
    )

    nuc_df_for_plot = df_nuc_pct_all.reset_index().rename(columns={'index': 'Nucleotide'})
    nuc_df_for_plot.insert(0, 'Batch', ref_name)
    mod_df_for_plot = modification_percentage_summary_df.copy()
    mod_df_for_plot.insert(0, 'Batch', ref_name)

    plot_root = _make_plot_root(
        ctx, '2a.' + _ref_plot_name(ctx) + 'Nucleotide_percentage_quilt',
    )

    return {
        'nuc_pct_df': nuc_df_for_plot,
        'mod_pct_df': mod_df_for_plot,
        'plot_root': plot_root,
        'save_also_png': ctx.save_png,
        'sgRNA_intervals': ref['sgRNA_intervals'],
        'sgRNA_names': ref['sgRNA_names'],
        'sgRNA_mismatches': ref['sgRNA_mismatches'],
        'sgRNA_sequences': ref['sgRNA_sequences'],
        'quantification_window_idxs': ref['include_idxs'],
        'custom_colors': ctx.custom_config.get('colors', {}),
    }


def prep_nucleotide_quilt_around_sgRNA(ctx: PlotContext):
    """Prepare kwargs for plot_nucleotide_quilt around one sgRNA (plot_2b).

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.

    Internally calls :func:`prep_nucleotide_quilt` to obtain the full-amplicon
    DataFrames, then slices them to the window around the cut point and
    adjusts sgRNA interval coordinates to the local frame.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    sgRNA_ind = ctx.sgRNA_ind
    cut_point = ref['sgRNA_cut_points'][sgRNA_ind]
    plot_half_window = max(1, ctx.args.plot_window_size)
    ref_len = ref['sequence_length']
    sgRNA_intervals = ref['sgRNA_intervals']
    include_idxs_list = ref['include_idxs']

    # Build full-amplicon DataFrames
    quilt_data = prep_nucleotide_quilt(ctx)
    nuc_df_for_plot = quilt_data['nuc_pct_df']
    mod_df_for_plot = quilt_data['mod_pct_df']

    # Columns 0 and 1 are 'Batch' and 'Nucleotide'/'Modification' labels;
    # sequence position p maps to column index p + 2.
    new_sel_cols_start = max(2, cut_point - plot_half_window + 1)
    new_sel_cols_end = min(ref_len, cut_point + plot_half_window + 1)
    sel_cols = [0, 1] + list(range(new_sel_cols_start + 2, new_sel_cols_end + 2))

    new_sgRNA_intervals = [
        (int_start - new_sel_cols_start, int_end - new_sel_cols_start)
        for (int_start, int_end) in sgRNA_intervals
    ]
    new_include_idx = [x - new_sel_cols_start for x in include_idxs_list]

    plot_root = _make_plot_root(
        ctx,
        '2b.' + _ref_plot_name(ctx) + 'Nucleotide_percentage_quilt_around_' + _sgRNA_label(ctx),
    )

    return {
        'nuc_pct_df': nuc_df_for_plot.iloc[:, sel_cols],
        'mod_pct_df': mod_df_for_plot.iloc[:, sel_cols],
        'plot_root': plot_root,
        'save_also_png': ctx.save_png,
        'sgRNA_intervals': new_sgRNA_intervals,
        'sgRNA_names': ref['sgRNA_names'],
        'sgRNA_mismatches': ref['sgRNA_mismatches'],
        'sgRNA_sequences': ref['sgRNA_sequences'],
        'quantification_window_idxs': new_include_idx,
        'custom_colors': ctx.custom_config.get('colors', {}),
    }


def prep_hdr_nucleotide_quilt(ctx: PlotContext):
    """Prepare kwargs for plot_nucleotide_quilt for HDR comparison (plot_4g).

    No scope fields required (iterates over all references with reads).

    Builds multi-batch ``nuc_pct_df`` and ``mod_pct_df`` across all references
    with non-zero read counts, aligned to the first reference's coordinate
    system.
    """
    ref_names_for_hdr = [r for r in ctx.ref_names if ctx.counts_total[r] > 0]
    ref0 = ref_names_for_hdr[0]
    ref_seq = ctx.refs[ref0]['sequence']
    seq_len = ctx.refs[ref0]['sequence_length']

    nuc_pcts = []
    for rn in ref_names_for_hdr:
        tot = float(ctx.counts_total[rn])
        for nuc in ['A', 'C', 'G', 'T', 'N', '-']:
            nuc_pcts.append(np.concatenate(
                ([rn, nuc],
                 np.array(ctx.ref1_all_base_count_vectors[rn + '_' + nuc]).astype(float) / tot),
            ))
    colnames = ['Batch', 'Nucleotide'] + list(ref_seq)
    nuc_pct_df = _to_numeric_ignore_columns(
        pd.DataFrame(nuc_pcts, columns=colnames), {'Batch', 'Nucleotide'},
    )

    mod_pcts = []
    for rn in ref_names_for_hdr:
        tot = float(ctx.counts_total[rn])
        mod_pcts.append(np.concatenate(([rn, 'Insertions'], np.array(ctx.ref1_all_insertion_count_vectors[rn]).astype(float) / tot)))
        mod_pcts.append(np.concatenate(([rn, 'Insertions_Left'], np.array(ctx.ref1_all_insertion_left_count_vectors[rn]).astype(float) / tot)))
        mod_pcts.append(np.concatenate(([rn, 'Deletions'], np.array(ctx.ref1_all_deletion_count_vectors[rn]).astype(float) / tot)))
        mod_pcts.append(np.concatenate(([rn, 'Substitutions'], np.array(ctx.ref1_all_substitution_count_vectors[rn]).astype(float) / tot)))
        mod_pcts.append(np.concatenate(([rn, 'All_modifications'], np.array(ctx.ref1_all_indelsub_count_vectors[rn]).astype(float) / tot)))
        mod_pcts.append(np.concatenate(([rn, 'Total'], [ctx.counts_total[rn]] * seq_len)))
    colnames = ['Batch', 'Modification'] + list(ref_seq)
    mod_pct_df = _to_numeric_ignore_columns(
        pd.DataFrame(mod_pcts, columns=colnames), {'Batch', 'Modification'},
    )

    plot_root = _make_plot_root(ctx, '4g.HDR_nucleotide_percentage_quilt')

    return {
        'nuc_pct_df': nuc_pct_df,
        'mod_pct_df': mod_pct_df,
        'plot_root': plot_root,
        'save_also_png': ctx.save_png,
        'sgRNA_intervals': ctx.refs[ref0]['sgRNA_intervals'],
        'sgRNA_names': ctx.refs[ref0]['sgRNA_names'],
        'sgRNA_mismatches': ctx.refs[ref0]['sgRNA_mismatches'],
        'sgRNA_sequences': ctx.refs[ref0]['sgRNA_sequences'],
        'quantification_window_idxs': [],  # windows may differ between amplicons
        'custom_colors': ctx.custom_config.get('colors', {}),
    }


def prep_pe_nucleotide_quilt(ctx: PlotContext):
    """Prepare kwargs for plot_nucleotide_quilt for PE comparison (plot_11a).

    No scope fields required.

    Same computation as ``prep_hdr_nucleotide_quilt`` but uses the first
    reference's ``include_idxs`` for ``quantification_window_idxs`` (PE
    amplicons share the same quantification window, unlike HDR).
    """
    result = prep_hdr_nucleotide_quilt(ctx)
    result['quantification_window_idxs'] = ctx.refs[ctx.ref_names[0]]['include_idxs']
    result['plot_root'] = _make_plot_root(
        ctx, '11a.Prime_editing_nucleotide_percentage_quilt',
    )
    return result


def prep_pe_nucleotide_quilt_around_sgRNA(ctx: PlotContext):
    """Prepare kwargs for plot_nucleotide_quilt around one sgRNA for PE (plot_11b).

    Requires ``ctx.sgRNA_ind``. Uses the first reference for sgRNA metadata.

    Internally calls :func:`prep_pe_nucleotide_quilt` to obtain the multi-batch
    DataFrames, then applies the same window-slicing as
    :func:`prep_nucleotide_quilt_around_sgRNA`.
    """
    ref0_name = ctx.ref_names[0]
    ref0 = ctx.refs[ref0_name]
    sgRNA_ind = ctx.sgRNA_ind
    cut_point = ref0['sgRNA_cut_points'][sgRNA_ind]
    plot_half_window = max(1, ctx.args.plot_window_size)
    ref_len = ref0['sequence_length']
    sgRNA_intervals = ref0['sgRNA_intervals']
    include_idxs_list = ref0['include_idxs']

    # Build full-amplicon PE DataFrames
    pe_data = prep_pe_nucleotide_quilt(ctx)
    nuc_df_for_plot = pe_data['nuc_pct_df']
    mod_df_for_plot = pe_data['mod_pct_df']

    # Window slicing (same logic as prep_nucleotide_quilt_around_sgRNA)
    new_sel_cols_start = max(2, cut_point - plot_half_window + 1)
    new_sel_cols_end = min(ref_len, cut_point + plot_half_window + 1)
    sel_cols = [0, 1] + list(range(new_sel_cols_start + 2, new_sel_cols_end + 2))

    new_sgRNA_intervals = [
        (int_start - new_sel_cols_start, int_end - new_sel_cols_start)
        for (int_start, int_end) in sgRNA_intervals
    ]
    new_include_idx = [x - new_sel_cols_start for x in include_idxs_list]

    # Build sgRNA label from first reference
    sgRNA = ref0['sgRNA_orig_sequences'][sgRNA_ind]
    sgRNA_name = ref0['sgRNA_names'][sgRNA_ind]
    from CRISPResso2 import CRISPRessoShared
    label = "sgRNA_" + sgRNA
    if sgRNA_name != "":
        label = sgRNA_name
    label = CRISPRessoShared.slugify(label)

    plot_root = _make_plot_root(
        ctx, '11b.Nucleotide_percentage_quilt_around_' + label,
    )

    return {
        'nuc_pct_df': nuc_df_for_plot.iloc[:, sel_cols],
        'mod_pct_df': mod_df_for_plot.iloc[:, sel_cols],
        'plot_root': plot_root,
        'save_also_png': ctx.save_png,
        'sgRNA_intervals': new_sgRNA_intervals,
        'sgRNA_names': pe_data['sgRNA_names'],
        'sgRNA_mismatches': pe_data['sgRNA_mismatches'],
        'sgRNA_sequences': pe_data['sgRNA_sequences'],
        'quantification_window_idxs': new_include_idx,
        'custom_colors': ctx.custom_config.get('colors', {}),
    }


def prep_global_frameshift_data(ctx: PlotContext):
    """Aggregate frameshift/splice/inframe counts across all coding-seq refs.

    No scope fields required.

    Used by plots 5a (global frameshift pie), 6a (global frameshift profiles),
    and 8a (global splice sites). For HDR refs, unmodified reads are added to
    the non-modified-non-frameshift bucket.

    Returns a dict with all aggregated values needed by the three plot
    functions.
    """
    global_MODIFIED_FRAMESHIFT = 0
    global_MODIFIED_NON_FRAMESHIFT = 0
    global_NON_MODIFIED_NON_FRAMESHIFT = 0
    global_SPLICING_SITES_MODIFIED = 0

    global_hists_frameshift = Counter()
    global_hists_frameshift[0] = 0
    global_hists_inframe = Counter()
    global_hists_inframe[0] = 0

    global_count_total = 0
    global_count_modified = 0
    global_count_unmodified = 0

    for ref_name in ctx.ref_names:
        if ctx.refs[ref_name]['contains_coding_seq']:
            global_MODIFIED_FRAMESHIFT += ctx.counts_modified_frameshift[ref_name]
            global_MODIFIED_NON_FRAMESHIFT += ctx.counts_modified_non_frameshift[ref_name]
            global_NON_MODIFIED_NON_FRAMESHIFT += ctx.counts_non_modified_non_frameshift[ref_name]
            global_SPLICING_SITES_MODIFIED += ctx.counts_splicing_sites_modified[ref_name]

            if ref_name == "HDR":
                # for HDR, add all unmodified reads to those that have
                # modifications not in exons
                global_NON_MODIFIED_NON_FRAMESHIFT += ctx.counts_unmodified[ref_name]

            for exon_len, count in ctx.hists_frameshift[ref_name].items():
                global_hists_frameshift[exon_len] += count
            for exon_len, count in ctx.hists_inframe[ref_name].items():
                global_hists_inframe[exon_len] += count

            global_count_total += ctx.counts_total[ref_name]
            global_count_modified += ctx.counts_modified[ref_name]
            global_count_unmodified += ctx.counts_unmodified[ref_name]

    return {
        'global_modified_frameshift': global_MODIFIED_FRAMESHIFT,
        'global_modified_non_frameshift': global_MODIFIED_NON_FRAMESHIFT,
        'global_non_modified_non_frameshift': global_NON_MODIFIED_NON_FRAMESHIFT,
        'global_splicing_sites_modified': global_SPLICING_SITES_MODIFIED,
        'global_hists_frameshift': global_hists_frameshift,
        'global_hists_inframe': global_hists_inframe,
        'global_count_total': global_count_total,
        'global_count_modified': global_count_modified,
        'global_count_unmodified': global_count_unmodified,
    }


def prep_global_modifications_reference(ctx: PlotContext):
    """Prepare kwargs for plot_global_modifications_reference (plot_4e/4f).

    Requires ``ctx.ref_name`` (either the primary amplicon or ``"HDR"``).

    The plot root and title are determined automatically based on whether
    ``ctx.ref_name`` is the first reference (4e) or ``"HDR"`` (4f).
    """
    ref_name = ctx.ref_name
    ref0 = ctx.ref_names[0]

    if ref_name == ref0:
        plot_root = _make_plot_root(ctx, '4e.' + ref0 + '.Global_mutations_in_all_reads')
        plot_title = 'Mutation position distribution in all reads with reference to %s' % ref0
    else:
        plot_root = _make_plot_root(
            ctx,
            '4f.' + ref0 + '.Global_mutations_in_HDR_reads_with_reference_to_' + ref0,
        )
        plot_title = 'Mutation position distribution in %s reads with reference to %s' % (ref_name, ref0)

    return {
        'ref1_all_insertion_count_vectors': ctx.ref1_all_insertion_count_vectors[ref_name],
        'ref1_all_deletion_count_vectors': ctx.ref1_all_deletion_count_vectors[ref_name],
        'ref1_all_substitution_count_vectors': ctx.ref1_all_substitution_count_vectors[ref_name],
        'ref1': ctx.refs[ref0],
        'include_idxs_list': _ref(ctx)['include_idxs'],
        'n_total': ctx.N_TOTAL,
        'ref_len': _ref(ctx)['sequence_length'],
        'ref_name': ref0,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'save_also_png': ctx.save_png,
        'plot_title': plot_title,
        'plot_root': plot_root,
    }


def prep_log_nuc_freqs(ctx: PlotContext):
    """Prepare kwargs for plot_log_nuc_freqs (plot_10d).

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.

    Computes ``plot_quant_window_idxs`` — the indices in the plotting window
    that overlap with the quantification window.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    sgRNA_ind = ctx.sgRNA_ind

    include_idxs_list = ref['include_idxs']
    plot_idxs = ref['sgRNA_plot_idxs'][sgRNA_ind]
    ref_len = ref['sequence_length']
    tot_aln_reads = ctx.counts_total[ref_name]

    # Build nucleotide frequency DataFrame and slice to plot window
    df_nuc_freq_all, _df_nuc_pct_all = _build_nuc_freq_df(ctx)
    df_nuc_freq = df_nuc_freq_all.iloc[:, plot_idxs]

    is_window = np.zeros(ref_len)
    for include_idx in include_idxs_list:
        is_window[include_idx] = 1

    plot_quant_window_idxs = []
    for plot_ind, loc in enumerate(plot_idxs):
        if is_window[loc]:
            plot_quant_window_idxs.append(plot_ind - 2)

    num_refs = len(ctx.ref_names)
    sgRNA_leg = _sgRNA_legend(ctx)

    plot_root = _make_plot_root(
        ctx,
        '10d.' + _ref_plot_name(ctx) + 'Log2_nucleotide_frequency_around_' + _sgRNA_label(ctx),
    )

    return {
        'df_nuc_freq': df_nuc_freq,
        'tot_aln_reads': tot_aln_reads,
        'plot_title': plot_title_with_ref_name(
            'Log2 Nucleotide Frequencies Around the ' + sgRNA_leg,
            ref_name,
            num_refs,
        ),
        'plot_root': plot_root,
        'save_also_png': ctx.save_png,
        'quantification_window_idxs': plot_quant_window_idxs,
    }


def prep_amino_acid_table(ctx: PlotContext):
    """Prepare amino acid table data for plot_9a and CSV export.

    Requires ``ctx.ref_name``, ``ctx.sgRNA_ind``, and ``ctx.coding_seq_ind``.

    Converts a nucleotide coding sequence to amino acids, computes the
    amino acid cut point, and builds the amino acid allele DataFrame via
    ``CRISPRessoShared.get_amino_acid_dataframe``.

    Returns a dict with:

    - ``coding_seq_amino_acids``: amino acid string for the coding sequence
    - ``amino_acid_cut_point``: cut point in amino acid coordinates
    - ``df_to_plot``: DataFrame of amino acid alleles (for CSV and plot)
    """
    from CRISPResso2 import CRISPRessoShared

    ref_name = ctx.ref_name
    ref = _ref(ctx)
    coding_seq_ind = ctx.coding_seq_ind

    coding_seqs = ctx.run_data['running_info']['coding_seqs']
    coding_seq = coding_seqs[coding_seq_ind]
    cut_point = ref['sgRNA_cut_points'][ctx.sgRNA_ind]
    exon_positions_start = ref['exon_positions'][0]
    df_alleles_for_ref = ctx.df_alleles.loc[ctx.df_alleles['Reference_Name'] == ref_name]
    exon_interval_start = ref['exon_intervals'][coding_seq_ind][0]

    _ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    blosum_path = os.path.join(_ROOT, "BLOSUM62")

    coding_seq_amino_acids = CRISPRessoShared.get_amino_acids_from_nucs(coding_seq)
    amino_acid_cut_point = (cut_point - exon_positions_start + 1) // 3
    df_to_plot = CRISPRessoShared.get_amino_acid_dataframe(
        df_alleles_for_ref,
        exon_interval_start,
        len(coding_seq_amino_acids),
        blosum_path,
        amino_acid_cut_point,
    )
    return {
        'coding_seq_amino_acids': coding_seq_amino_acids,
        'amino_acid_cut_point': amino_acid_cut_point,
        'df_to_plot': df_to_plot,
    }


def _prep_windowed_alleles(
    df_alleles_around_cut,
    cut_point,
    window_left,
    window_right,
    ref_sequence,
    sgRNA_intervals,
    count_total,
    allele_plot_pcts_only_for_assigned_reference,
    expand_allele_plots_by_quantification,
):
    """Shared logic for ``prep_alleles_around_cut`` and ``prep_base_edit_quilt``.

    .. warning::

        **Mutates** *df_alleles_around_cut* in place when
        *allele_plot_pcts_only_for_assigned_reference* is True (adds a
        ``%AllReads`` column and overwrites ``%Reads``).  CORE relies on
        this mutation — it reads the modified DataFrame back from the
        returned dict for CSV writing.

    Applies:

    1. Optional percentage adjustment when
       *allele_plot_pcts_only_for_assigned_reference* is True
    2. Reference sequence slicing for the window
    3. Optional groupby collapse (when *expand_allele_plots_by_quantification*
       is False)
    4. sgRNA interval coordinate adjustment to the local window frame

    Returns ``(df_alleles_around_cut, df_to_plot, ref_seq_around_cut,
    new_sgRNA_intervals, new_sel_cols_start)``.
    """
    if allele_plot_pcts_only_for_assigned_reference:
        df_alleles_around_cut['%AllReads'] = df_alleles_around_cut['%Reads']
        df_alleles_around_cut['%Reads'] = df_alleles_around_cut['#Reads'] / count_total * 100

    ref_seq_around_cut = ref_sequence[
        cut_point - window_left + 1:cut_point + window_right + 1
    ]

    df_to_plot = df_alleles_around_cut
    if not expand_allele_plots_by_quantification:
        df_to_plot = df_alleles_around_cut.groupby(
            ['Aligned_Sequence', 'Reference_Sequence'],
        ).sum().reset_index().set_index('Aligned_Sequence')
        df_to_plot.sort_values(
            by=['#Reads', 'Aligned_Sequence', 'Reference_Sequence'],
            inplace=True,
            ascending=[False, True, True],
        )

    new_sgRNA_intervals = []
    new_sel_cols_start = cut_point - window_left
    for (int_start, int_end) in sgRNA_intervals:
        new_sgRNA_intervals.append(
            (int_start - new_sel_cols_start - 1, int_end - new_sel_cols_start - 1),
        )

    return (
        df_alleles_around_cut,
        df_to_plot,
        ref_seq_around_cut,
        new_sgRNA_intervals,
        new_sel_cols_start,
    )


def prep_alleles_around_cut(ctx: PlotContext):
    """Prepare alleles-around-cut data for plot_9 and CSV export.

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.

    Computes asymmetric window sizes, calls
    ``CRISPRessoShared.get_dataframe_around_cut_asymmetrical``, then applies
    percentage adjustment, optional groupby collapse, and sgRNA interval
    recomputation via ``_prep_windowed_alleles``.

    .. warning::

        **Mutates** the intermediate ``df_alleles_around_cut`` in place when
        ``args.allele_plot_pcts_only_for_assigned_reference`` is True.
        See ``_prep_windowed_alleles`` for details.

    Returns a dict with:

    - ``df_alleles_around_cut``: the (possibly mutated) alleles DataFrame
    - ``df_to_plot``: after optional groupby collapse (for ``prep_alleles_table``)
    - ``ref_seq_around_cut``: reference sequence in the window
    - ``new_sgRNA_intervals``: sgRNA intervals in local coordinates
    - ``new_cut_point``: cut point in local coordinates, or None if the cut
      point is not inside any sgRNA interval
    - ``window_truncated``: True if the window was truncated due to proximity
      to amplicon edges
    - ``plot_half_window_left``: computed left window size
    - ``plot_half_window_right``: computed right window size
    """
    from CRISPResso2 import CRISPRessoShared

    ref_name = ctx.ref_name
    ref = _ref(ctx)
    sgRNA_ind = ctx.sgRNA_ind

    cut_point = ref['sgRNA_cut_points'][sgRNA_ind]
    ref_len = ref['sequence_length']
    plot_window_size = ctx.args.plot_window_size

    plot_half_window_left, plot_half_window_right, window_truncated = \
        _compute_half_windows(cut_point, plot_window_size, ref_len)

    df_alleles_around_cut = CRISPRessoShared.get_dataframe_around_cut_asymmetrical(
        ctx.df_alleles.loc[ctx.df_alleles['Reference_Name'] == ref_name],
        cut_point,
        plot_half_window_left,
        plot_half_window_right,
    )

    (
        df_alleles_around_cut,
        df_to_plot,
        ref_seq_around_cut,
        new_sgRNA_intervals,
        new_sel_cols_start,
    ) = _prep_windowed_alleles(
        df_alleles_around_cut=df_alleles_around_cut,
        cut_point=cut_point,
        window_left=plot_half_window_left,
        window_right=plot_half_window_right,
        ref_sequence=ref['sequence'],
        sgRNA_intervals=ref['sgRNA_intervals'],
        count_total=ctx.counts_total[ref_name],
        allele_plot_pcts_only_for_assigned_reference=ctx.args.allele_plot_pcts_only_for_assigned_reference,
        expand_allele_plots_by_quantification=ctx.args.expand_allele_plots_by_quantification,
    )

    new_cut_point = None
    for (int_start, int_end) in ref['sgRNA_intervals']:
        if int_start <= cut_point <= int_end:
            new_cut_point = cut_point - new_sel_cols_start - 1

    return {
        'df_alleles_around_cut': df_alleles_around_cut,
        'df_to_plot': df_to_plot,
        'ref_seq_around_cut': ref_seq_around_cut,
        'new_sgRNA_intervals': new_sgRNA_intervals,
        'new_cut_point': new_cut_point,
        'window_truncated': window_truncated,
        'plot_half_window_left': plot_half_window_left,
        'plot_half_window_right': plot_half_window_right,
    }


def prep_base_edit_quilt(ctx: PlotContext):
    """Prepare base edit quilt data for plot_10h and CSV export.

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.

    Internally calls ``CRISPRessoShared.get_base_edit_dataframe_around_cut``
    to produce the sliced DataFrame, then applies the same windowing logic
    as ``prep_alleles_around_cut`` (via ``_prep_windowed_alleles``).

    Uses a symmetric window (``args.plot_window_size``) and computes
    ``x_labels`` — the 1-indexed positions of the conversion nucleotide in
    the full reference sequence.

    .. warning::

        **Mutates** the intermediate ``df_alleles_around_cut`` in place when
        ``args.allele_plot_pcts_only_for_assigned_reference`` is True.
        See ``_prep_windowed_alleles`` for details.

    Returns a dict with:

    - ``df_alleles_around_cut``: the (possibly mutated) alleles DataFrame
    - ``df_to_plot``: after optional groupby collapse (for ``prep_alleles_table``)
    - ``ref_seq_around_cut``: reference sequence in the window
    - ``new_sgRNA_intervals``: sgRNA intervals in local coordinates
    - ``x_labels``: 1-indexed positions of ``conversion_nuc_from`` in the
      full reference
    """
    from CRISPResso2 import CRISPRessoShared

    ref_name = ctx.ref_name
    ref = _ref(ctx)

    cut_point = ref['sgRNA_cut_points'][ctx.sgRNA_ind]
    plot_half_window = max(1, ctx.args.plot_window_size)
    conversion_nuc_from = ctx.args.conversion_nuc_from

    df_alleles_around_cut = CRISPRessoShared.get_base_edit_dataframe_around_cut(
        ctx.df_alleles.loc[ctx.df_alleles['Reference_Name'] == ref_name],
        conversion_nuc_from,
    )

    (
        df_alleles_around_cut,
        df_to_plot,
        ref_seq_around_cut,
        new_sgRNA_intervals,
        _new_sel_cols_start,
    ) = _prep_windowed_alleles(
        df_alleles_around_cut=df_alleles_around_cut,
        cut_point=cut_point,
        window_left=plot_half_window,
        window_right=plot_half_window,
        ref_sequence=ref['sequence'],
        sgRNA_intervals=ref['sgRNA_intervals'],
        count_total=ctx.counts_total[ref_name],
        allele_plot_pcts_only_for_assigned_reference=ctx.args.allele_plot_pcts_only_for_assigned_reference,
        expand_allele_plots_by_quantification=ctx.args.expand_allele_plots_by_quantification,
    )

    x_labels = [
        ind for ind, a in enumerate(ref['sequence'], start=1)
        if a == conversion_nuc_from
    ]

    return {
        'df_alleles_around_cut': df_alleles_around_cut,
        'df_to_plot': df_to_plot,
        'ref_seq_around_cut': ref_seq_around_cut,
        'new_sgRNA_intervals': new_sgRNA_intervals,
        'x_labels': x_labels,
    }


def prep_alternate_allele_counts(sub_base_vectors, ref_name, ref_sequence):
    """Aggregate per-position substitution counts by reference nucleotide.

    For each reference base, sums the substitution counts across all
    positions where that base appears. Returns a nested dict
    ``{ref_nuc: {obs_nuc: count}}``.

    This is a low-level utility (not PlotContext-based) used by
    ``count_alternate_alleles`` in CORE.
    """
    alph = ['A', 'C', 'G', 'T', 'N']
    alt_nuc_counts = {a: {b: 0 for b in alph} for a in alph}
    for pos, ref_base in enumerate(ref_sequence):
        ref_base = ref_base.upper()
        if ref_base not in alt_nuc_counts:
            continue
        for nuc in alph:
            alt_nuc_counts[ref_base][nuc] += int(sub_base_vectors[ref_name + '_' + nuc][pos])
    return alt_nuc_counts


def prep_class_piechart_and_barplot(class_counts_order, class_counts,
                                    ref_names, expected_hdr_amplicon_seq,
                                    N_TOTAL):
    """Compute labels and sizes for class piechart/barplot (plot_1b/1c).

    This is a standalone utility (not PlotContext-based) that transforms
    raw class counts into display-ready labels and percentage sizes.

    Returns dict with 'labels', 'sizes', 'N_TOTAL'.
    """
    labels = []
    sizes = []
    for class_name in class_counts_order:
        if expected_hdr_amplicon_seq != "" and class_name == ref_names[0] + "_MODIFIED":
            labels.append("NHEJ\n(" + str(class_counts[class_name]) + " reads)")
        elif expected_hdr_amplicon_seq != "" and class_name == "HDR_MODIFIED":
            labels.append("Imperfect HDR\n(" + str(class_counts[class_name]) + " reads)")
        elif expected_hdr_amplicon_seq != "" and class_name == "HDR_UNMODIFIED":
            labels.append("HDR\n(" + str(class_counts[class_name]) + " reads)")
        else:
            display_class_name = class_name
            if len(ref_names) == 1:
                display_class_name = display_class_name.replace('Reference_', '')
            labels.append(display_class_name + "\n(" + str(class_counts[class_name]) + " reads)")
        sizes.append(100 * class_counts[class_name] / float(N_TOTAL))
    return {'labels': labels, 'sizes': sizes, 'N_TOTAL': N_TOTAL}


def prep_base_edit_upset(ref_seq, df_alleles, ref_name, sgRNA_interval,
                         args):
    """Prepare data for base edit upset plot (plot_10i).

    Finds the target sequence, aligns it to reference, identifies
    base changes, and computes combination counts.

    This is a standalone utility (not PlotContext-based) that encapsulates
    the inline logic from CRISPRessoCORE for plot_10i.

    Parameters
    ----------
    ref_seq : str
        Reference amplicon sequence.
    df_alleles : pd.DataFrame
        Allele table (full, not filtered by ref_name).
    ref_name : str
        Name of the reference (wild-type) amplicon.
    sgRNA_interval : tuple
        (start, end) of the sgRNA interval for this sgRNA.
    args : argparse.Namespace
        Must have: base_editor_target_ref_skip_allele_count,
        base_editor_consider_changes_outside_qw,
        needleman_wunsch_gap_open, needleman_wunsch_gap_extend,
        needleman_wunsch_aln_matrix_loc, quantification_window_coordinates.

    Returns
    -------
    dict or None
        Dict with 'bp_substitutions_arr' and 'binary_allele_counts',
        or None if no target sequence found, quantification_window_coordinates
        is set, or no substitutions detected.
    """
    from CRISPResso2.CRISPRessoCORE import (
        get_base_edit_target_sequence,
        get_bp_substitutions,
        get_refpos_values,
        get_upset_plot_counts,
    )

    target_seq = get_base_edit_target_sequence(
        ref_seq, df_alleles,
        args.base_editor_target_ref_skip_allele_count,
    )

    if not target_seq:
        return None
    if args.quantification_window_coordinates is not None:
        return None

    # Lazy import — CRISPResso2Align is a compiled C extension
    import CRISPResso2Align

    # Align target to reference
    aln_matrix_loc = args.needleman_wunsch_aln_matrix_loc
    if aln_matrix_loc == 'EDNAFULL':
        aln_matrix = CRISPResso2Align.make_matrix()
    else:
        import os
        if not os.path.exists(aln_matrix_loc):
            return None
        aln_matrix = CRISPResso2Align.read_matrix(aln_matrix_loc)

    # gap_incentive comes from the ref data; use zeros as default
    gap_incentive = np.zeros(len(ref_seq) + 1, dtype=np.int32)

    aln_target_seq, aln_ref_seq, _aln_score = CRISPResso2Align.global_align(
        target_seq,
        ref_seq,
        matrix=aln_matrix,
        gap_incentive=gap_incentive,
        gap_open=args.needleman_wunsch_gap_open,
        gap_extend=args.needleman_wunsch_gap_extend,
    )

    # Determine which positions to include
    if args.base_editor_consider_changes_outside_qw:
        ref_positions_to_include = list(range(len(ref_seq)))
    else:
        this_start, this_stop = sgRNA_interval
        ref_positions_to_include = list(range(this_start, this_stop + 1))

    ref_changes_dict = get_refpos_values(aln_ref_seq, aln_target_seq)
    bp_substitutions_arr = get_bp_substitutions(
        ref_changes_dict, ref_seq, ref_positions_to_include,
    )

    if len(bp_substitutions_arr) == 0:
        return None

    counts_dict = get_upset_plot_counts(
        df_alleles, bp_substitutions_arr, ref_name,
    )

    return {
        'bp_substitutions_arr': bp_substitutions_arr,
        'binary_allele_counts': counts_dict['binary_allele_counts'],
    }


def prep_conversion_at_sel_nucs(ctx: PlotContext):
    """Compute from_nuc_indices and selected-nucleotide DataFrames for plots 10e/10f/10g.

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.

    Finds the column positions matching ``conversion_nuc_from`` and returns
    sliced, re-labeled DataFrames used for both CSV writes and plot inputs.

    Returns a dict with:
    - ``from_nuc_indices``: list of int column positions matching the conversion nucleotide
    - ``just_sel_nuc_pcts``: DataFrame sliced to those columns, columns renamed
    - ``just_sel_nuc_freqs``: same slicing applied to frequency DataFrame
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    sgRNA_ind = ctx.sgRNA_ind
    conversion_nuc_from = ctx.args.conversion_nuc_from

    plot_idxs = ref['sgRNA_plot_idxs'][sgRNA_ind]

    # Build nucleotide freq/pct DataFrames and slice to plot window
    df_nuc_freq_all, df_nuc_pct_all = _build_nuc_freq_df(ctx)
    plot_nuc_pcts = df_nuc_pct_all.iloc[:, plot_idxs]
    plot_nuc_freqs = df_nuc_freq_all.iloc[:, plot_idxs]

    from_nuc_indices = [
        pos for pos, char in enumerate(list(plot_nuc_pcts.columns.values))
        if char == conversion_nuc_from
    ]

    just_sel_nuc_pcts = plot_nuc_pcts.iloc[:, from_nuc_indices].copy()
    just_sel_nuc_pcts.columns = [
        conversion_nuc_from + str(pos + 1) for pos in from_nuc_indices
    ]

    just_sel_nuc_freqs = plot_nuc_freqs.iloc[:, from_nuc_indices].copy()
    just_sel_nuc_freqs.columns = [
        conversion_nuc_from + str(pos + 1) for pos in from_nuc_indices
    ]

    return {
        'from_nuc_indices': from_nuc_indices,
        'just_sel_nuc_pcts': just_sel_nuc_pcts,
        'just_sel_nuc_freqs': just_sel_nuc_freqs,
    }
