"""Extracted data preparation functions for CRISPResso2 plots.

Each function takes raw analysis data and returns the kwargs dict
that the corresponding CRISPRessoPlot function expects. These are
pure functions — no file I/O, no side effects.

Only functions with non-trivial computation are extracted here.
Trivial pass-through plots (where the dict is just packaging existing
variables) stay inline in CORE — there's nothing to reuse.

CORE calls these to build plot inputs. CRISPRessoPro can call them
from PlotContext to generate plots independently.
"""

from collections import Counter

import numpy as np
import pandas as pd


def _to_numeric_ignore_columns(df, ignore_columns):
    """Convert DataFrame columns to numeric, ignoring specified columns."""
    for col in df.columns:
        if col not in ignore_columns:
            df[col] = df[col].apply(pd.to_numeric, errors='raise')
    return df


def _plot_title_with_ref_name(title, ref_name, num_refs):
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


def prep_indel_size_distribution(hdensity, hlengths, center_index,
                                 n_this_category, ref_name, ref_names,
                                 plot_root, save_also_png,
                                 plot_histogram_outliers):
    """Prepare kwargs for plot_indel_size_distribution (plot_3a).

    Computes xmin/xmax by clipping to the 99th percentile of the
    density distribution (unless ``plot_histogram_outliers`` is True),
    then clamping to at least ±15.
    """
    xmin = min(hlengths)
    xmax = max(hlengths)

    if not plot_histogram_outliers:
        xmax = _clip_to_percentile(hlengths, hdensity, 0.99)
        xmin = _clip_to_percentile(hlengths[::-1], hdensity[::-1], 0.99)

    xmin = min(xmin, -15)
    xmax = max(xmax, 15)

    num_refs = len(ref_names)

    return {
        'hdensity': hdensity,
        'hlengths': hlengths,
        'center_index': center_index,
        'n_this_category': n_this_category,
        'xmin': xmin,
        'xmax': xmax,
        'title': _plot_title_with_ref_name(
            'Indel size distribution', ref_name, num_refs,
        ),
        'plot_root': plot_root,
        'save_also_png': save_also_png,
        'ref_name': ref_name,
    }


def prep_frequency_deletions_insertions(x_bins_ins, y_values_ins,
                                        x_bins_del, y_values_del,
                                        x_bins_mut, y_values_mut,
                                        hdensity, ref, counts_total,
                                        ref_name, ref_names, plot_root,
                                        save_also_png, custom_colors,
                                        plot_histogram_outliers):
    """Prepare kwargs for plot_frequency_deletions_insertions (plot_3b).

    Clips each histogram's xmax to the 99th percentile of the overall
    density (unless ``plot_histogram_outliers`` is True), clamped to 15.
    """
    xmax_ins = max(x_bins_ins)
    xmax_del = max(x_bins_del)
    xmax_mut = max(x_bins_mut)

    if not plot_histogram_outliers:
        density_total = hdensity.sum()
        xmax_ins = _clip_to_percentile(x_bins_ins, y_values_ins, 0.99, total=density_total)
        xmax_del = _clip_to_percentile(x_bins_del, y_values_del, 0.99, total=density_total)
        xmax_mut = _clip_to_percentile(x_bins_mut, y_values_mut, 0.99, total=density_total)

    xmax_ins = max(15, xmax_ins)
    xmax_del = max(15, xmax_del)
    xmax_mut = max(15, xmax_mut)

    num_refs = len(ref_names)

    return {
        'ref': ref,
        'counts_total': counts_total,
        'plot_path': plot_root,
        'plot_titles': {
            'ins': _plot_title_with_ref_name('Insertions', ref_name, num_refs),
            'del': _plot_title_with_ref_name('Deletions', ref_name, num_refs),
            'mut': _plot_title_with_ref_name('Substitutions', ref_name, num_refs),
        },
        'xmax_ins': xmax_ins,
        'xmax_del': xmax_del,
        'xmax_mut': xmax_mut,
        'save_also_png': save_also_png,
        'custom_colors': custom_colors,
        'ref_name': ref_name,
    }


def prep_amplicon_modifications(all_indelsub_count_vector, include_idxs_list,
                                cut_points, plot_cut_points, sgRNA_intervals,
                                N_TOTAL, n_this_category, ref_name, ref_names,
                                ref_len, plot_root, custom_colors,
                                save_also_png):
    """Prepare kwargs for plot_amplicon_modifications (plot_4a).

    Pattern A — packaging with y_max computation and title generation.
    """
    num_refs = len(ref_names)
    y_max = max(all_indelsub_count_vector) * 1.1

    return {
        'all_indelsub_count_vectors': all_indelsub_count_vector,
        'include_idxs_list': include_idxs_list,
        'cut_points': cut_points,
        'plot_cut_points': plot_cut_points,
        'sgRNA_intervals': sgRNA_intervals,
        'n_total': N_TOTAL,
        'n_this_category': n_this_category,
        'ref_name': ref_name,
        'num_refs': num_refs,
        'ref_len': ref_len,
        'y_max': y_max,
        'plot_titles': {
            'combined': _plot_title_with_ref_name(
                'Combined Insertions/Deletions/Substitutions', ref_name, num_refs,
            ),
            'main': _plot_title_with_ref_name(
                'Mutation position distribution', ref_name, num_refs,
            ),
        },
        'plot_root': plot_root,
        'custom_colors': custom_colors,
        'save_also_png': save_also_png,
    }


def prep_modification_frequency(include_idxs_list, all_insertion_count_vector,
                                all_deletion_count_vector,
                                all_substitution_count_vector,
                                sgRNA_intervals, ref_len, ref_name, ref_names,
                                N_TOTAL, n_this_category, cut_points,
                                plot_cut_points, y_max, plot_root,
                                custom_colors, save_also_png):
    """Prepare kwargs for plot_modification_frequency (plot_4b).

    Pattern A — packaging with title generation.
    """
    num_refs = len(ref_names)

    return {
        'include_idxs_list': include_idxs_list,
        'all_insertion_count_vectors': all_insertion_count_vector,
        'all_deletion_count_vectors': all_deletion_count_vector,
        'all_substitution_count_vectors': all_substitution_count_vector,
        'sgRNA_intervals': sgRNA_intervals,
        'ref_len': ref_len,
        'ref_name': ref_name,
        'num_refs': num_refs,
        'n_total': N_TOTAL,
        'n_this_category': n_this_category,
        'cut_points': cut_points,
        'plot_cut_points': plot_cut_points,
        'y_max': y_max,
        'plot_title': _plot_title_with_ref_name(
            'Mutation position distribution', ref_name, num_refs,
        ),
        'plot_root': plot_root,
        'custom_colors': custom_colors,
        'save_also_png': save_also_png,
    }


def prep_dsODN_piechart(df_alleles, N_TOTAL, plot_root, save_also_png):
    """Prepare kwargs for plot_class_dsODN_piechart (plot_1d).

    Computes labels and sizes from df_alleles 'contains dsODN' column.
    """
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
        'save_also_png': save_also_png,
    }


def prep_nucleotide_quilt(
    all_base_count_vectors,
    all_insertion_count_vectors,
    all_insertion_left_count_vectors,
    all_deletion_count_vectors,
    all_substitution_count_vectors,
    all_indelsub_count_vectors,
    counts_total,
    ref_name,
    ref_seq,
    include_idxs_list,
    sgRNA_intervals,
    sgRNA_names,
    sgRNA_mismatches,
    sgRNA_sequences,
    plot_root,
    save_also_png,
    custom_colors,
):
    """Prepare kwargs for plot_nucleotide_quilt (plot_2a).

    Builds ``nuc_pct_df`` and ``mod_pct_df`` from per-position count vectors.

    The returned dict is passed directly to ``plot_nucleotide_quilt``.
    CORE should also read ``nuc_pct_df`` and ``mod_pct_df`` back from the
    returned dict to pass as inputs to ``prep_nucleotide_quilt_around_sgRNA``.
    """
    tot = float(counts_total)

    df_nuc_freq_all = pd.DataFrame([
        all_base_count_vectors[ref_name + '_A'],
        all_base_count_vectors[ref_name + '_C'],
        all_base_count_vectors[ref_name + '_G'],
        all_base_count_vectors[ref_name + '_T'],
        all_base_count_vectors[ref_name + '_N'],
        all_base_count_vectors[ref_name + '_-'],
    ])
    df_nuc_freq_all.index = ['A', 'C', 'G', 'T', 'N', '-']
    df_nuc_freq_all.columns = list(ref_seq)
    df_nuc_pct_all = df_nuc_freq_all.divide(tot)

    mod_pcts = []
    mod_pcts.append(np.concatenate((['Insertions'], np.array(all_insertion_count_vectors).astype(float) / tot)))
    mod_pcts.append(np.concatenate((['Insertions_Left'], np.array(all_insertion_left_count_vectors).astype(float) / tot)))
    mod_pcts.append(np.concatenate((['Deletions'], np.array(all_deletion_count_vectors).astype(float) / tot)))
    mod_pcts.append(np.concatenate((['Substitutions'], np.array(all_substitution_count_vectors).astype(float) / tot)))
    mod_pcts.append(np.concatenate((['All_modifications'], np.array(all_indelsub_count_vectors).astype(float) / tot)))
    mod_pcts.append(np.concatenate((['Total'], [counts_total] * len(ref_seq))))
    colnames = ['Modification'] + list(ref_seq)
    modification_percentage_summary_df = _to_numeric_ignore_columns(
        pd.DataFrame(mod_pcts, columns=colnames), {'Modification'},
    )

    nuc_df_for_plot = df_nuc_pct_all.reset_index().rename(columns={'index': 'Nucleotide'})
    nuc_df_for_plot.insert(0, 'Batch', ref_name)
    mod_df_for_plot = modification_percentage_summary_df.copy()
    mod_df_for_plot.insert(0, 'Batch', ref_name)

    return {
        'nuc_pct_df': nuc_df_for_plot,
        'mod_pct_df': mod_df_for_plot,
        'fig_filename_root': plot_root,
        'save_also_png': save_also_png,
        'sgRNA_intervals': sgRNA_intervals,
        'sgRNA_names': sgRNA_names,
        'sgRNA_mismatches': sgRNA_mismatches,
        'sgRNA_sequences': sgRNA_sequences,
        'quantification_window_idxs': include_idxs_list,
        'custom_colors': custom_colors,
    }


def prep_nucleotide_quilt_around_sgRNA(
    nuc_df_for_plot,
    mod_df_for_plot,
    cut_point,
    plot_half_window,
    ref_len,
    sgRNA_intervals,
    include_idxs_list,
    sgRNA_names,
    sgRNA_mismatches,
    sgRNA_sequences,
    plot_root,
    save_also_png,
    custom_colors,
):
    """Prepare kwargs for plot_nucleotide_quilt around one sgRNA (plot_2b).

    Slices ``nuc_df_for_plot`` and ``mod_df_for_plot`` to the window around
    ``cut_point`` and adjusts sgRNA interval coordinates to the local frame.

    ``nuc_df_for_plot`` and ``mod_df_for_plot`` come from the ``nuc_pct_df``
    and ``mod_pct_df`` keys of ``prep_nucleotide_quilt``'s return value.
    """
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

    return {
        'nuc_pct_df': nuc_df_for_plot.iloc[:, sel_cols],
        'mod_pct_df': mod_df_for_plot.iloc[:, sel_cols],
        'fig_filename_root': plot_root,
        'save_also_png': save_also_png,
        'sgRNA_intervals': new_sgRNA_intervals,
        'sgRNA_names': sgRNA_names,
        'sgRNA_mismatches': sgRNA_mismatches,
        'sgRNA_sequences': sgRNA_sequences,
        'quantification_window_idxs': new_include_idx,
        'custom_colors': custom_colors,
    }


def prep_hdr_nucleotide_quilt(
    ref_names_for_hdr,
    counts_total,
    refs,
    ref1_all_base_count_vectors,
    ref1_all_insertion_count_vectors,
    ref1_all_insertion_left_count_vectors,
    ref1_all_deletion_count_vectors,
    ref1_all_substitution_count_vectors,
    ref1_all_indelsub_count_vectors,
    custom_colors,
    save_also_png,
):
    """Prepare kwargs for plot_nucleotide_quilt for HDR comparison (plot_4g).

    Builds multi-batch ``nuc_pct_df`` and ``mod_pct_df`` across all references
    with non-zero read counts, aligned to the first reference's coordinate
    system.

    Parameters
    ----------
    ref_names_for_hdr : list of str
        Reference names with counts > 0.
    counts_total : dict
        ref_name → total read count.
    refs : dict
        ref_name → dict with 'sequence', 'sequence_length', sgRNA metadata.
    ref1_all_base_count_vectors : dict
        Keys like ``"FANC_A"``, ``"HDR_C"`` → numpy array of per-position
        base counts aligned to ref1's coordinate system.
    ref1_all_insertion_count_vectors, ref1_all_insertion_left_count_vectors,
    ref1_all_deletion_count_vectors, ref1_all_substitution_count_vectors,
    ref1_all_indelsub_count_vectors : dict
        ref_name → numpy array of per-position modification counts aligned
        to ref1's coordinate system.
    """
    ref0 = ref_names_for_hdr[0]
    ref_seq = refs[ref0]['sequence']
    seq_len = refs[ref0]['sequence_length']

    nuc_pcts = []
    for ref_name in ref_names_for_hdr:
        tot = float(counts_total[ref_name])
        for nuc in ['A', 'C', 'G', 'T', 'N', '-']:
            nuc_pcts.append(np.concatenate(
                ([ref_name, nuc],
                 np.array(ref1_all_base_count_vectors[ref_name + '_' + nuc]).astype(float) / tot),
            ))
    colnames = ['Batch', 'Nucleotide'] + list(ref_seq)
    nuc_pct_df = _to_numeric_ignore_columns(
        pd.DataFrame(nuc_pcts, columns=colnames), {'Batch', 'Nucleotide'},
    )

    mod_pcts = []
    for ref_name in ref_names_for_hdr:
        tot = float(counts_total[ref_name])
        mod_pcts.append(np.concatenate(([ref_name, 'Insertions'], np.array(ref1_all_insertion_count_vectors[ref_name]).astype(float) / tot)))
        mod_pcts.append(np.concatenate(([ref_name, 'Insertions_Left'], np.array(ref1_all_insertion_left_count_vectors[ref_name]).astype(float) / tot)))
        mod_pcts.append(np.concatenate(([ref_name, 'Deletions'], np.array(ref1_all_deletion_count_vectors[ref_name]).astype(float) / tot)))
        mod_pcts.append(np.concatenate(([ref_name, 'Substitutions'], np.array(ref1_all_substitution_count_vectors[ref_name]).astype(float) / tot)))
        mod_pcts.append(np.concatenate(([ref_name, 'All_modifications'], np.array(ref1_all_indelsub_count_vectors[ref_name]).astype(float) / tot)))
        mod_pcts.append(np.concatenate(([ref_name, 'Total'], [counts_total[ref_name]] * seq_len)))
    colnames = ['Batch', 'Modification'] + list(ref_seq)
    mod_pct_df = _to_numeric_ignore_columns(
        pd.DataFrame(mod_pcts, columns=colnames), {'Batch', 'Modification'},
    )

    return {
        'nuc_pct_df': nuc_pct_df,
        'mod_pct_df': mod_pct_df,
        'fig_filename_root': '',  # caller sets this
        'save_also_png': save_also_png,
        'sgRNA_intervals': refs[ref0]['sgRNA_intervals'],
        'sgRNA_names': refs[ref0]['sgRNA_names'],
        'sgRNA_mismatches': refs[ref0]['sgRNA_mismatches'],
        'sgRNA_sequences': refs[ref0]['sgRNA_sequences'],
        'quantification_window_idxs': [],  # windows may differ between amplicons
        'custom_colors': custom_colors,
    }


def prep_global_frameshift_data(
    ref_names,
    refs,
    counts_modified_frameshift,
    counts_modified_non_frameshift,
    counts_non_modified_non_frameshift,
    counts_splicing_sites_modified,
    counts_total,
    counts_modified,
    counts_unmodified,
    hists_frameshift,
    hists_inframe,
):
    """Aggregate frameshift/splice/inframe counts across all coding-seq refs.

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

    for ref_name in ref_names:
        if refs[ref_name]['contains_coding_seq']:
            global_MODIFIED_FRAMESHIFT += counts_modified_frameshift[ref_name]
            global_MODIFIED_NON_FRAMESHIFT += counts_modified_non_frameshift[ref_name]
            global_NON_MODIFIED_NON_FRAMESHIFT += counts_non_modified_non_frameshift[ref_name]
            global_SPLICING_SITES_MODIFIED += counts_splicing_sites_modified[ref_name]

            if ref_name == "HDR":
                # for HDR, add all unmodified reads to those that have
                # modifications not in exons
                global_NON_MODIFIED_NON_FRAMESHIFT += counts_unmodified[ref_name]

            for exon_len, count in hists_frameshift[ref_name].items():
                global_hists_frameshift[exon_len] += count
            for exon_len, count in hists_inframe[ref_name].items():
                global_hists_inframe[exon_len] += count

            global_count_total += counts_total[ref_name]
            global_count_modified += counts_modified[ref_name]
            global_count_unmodified += counts_unmodified[ref_name]

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


def prep_global_modifications_reference(
    ref_name,
    ref_names,
    ref1_all_insertion_count_vectors,
    ref1_all_deletion_count_vectors,
    ref1_all_substitution_count_vectors,
    ref1,
    include_idxs_list,
    N_TOTAL,
    ref_len,
    custom_colors,
    save_also_png,
    _jp,
):
    """Prepare kwargs for plot_global_modifications_reference (plot_4e/4f).

    Selects plot_title and plot_root based on whether ref_name is the
    primary reference (4e) or the HDR reference (4f).

    Parameters
    ----------
    _jp : callable
        Path-joining function that prepends the output directory.
    """
    ref0 = ref_names[0]

    if ref_name == ref0:
        plot_root = _jp('4e.' + ref0 + '.Global_mutations_in_all_reads')
        plot_title = (
            'Mutation position distribution in all reads with reference to %s'
            % ref0
        )
    else:  # ref_name == "HDR"
        plot_root = _jp(
            '4f.' + ref0
            + '.Global_mutations_in_HDR_reads_with_reference_to_' + ref0
        )
        plot_title = (
            'Mutation position distribution in %s reads with reference to %s'
            % (ref_name, ref0)
        )

    return {
        'ref1_all_insertion_count_vectors': ref1_all_insertion_count_vectors,
        'ref1_all_deletion_count_vectors': ref1_all_deletion_count_vectors,
        'ref1_all_substitution_count_vectors': ref1_all_substitution_count_vectors,
        'ref1': ref1,
        'include_idxs_list': include_idxs_list,
        'n_total': N_TOTAL,
        'ref_len': ref_len,
        'ref_name': ref0,
        'custom_colors': custom_colors,
        'save_also_png': save_also_png,
        'plot_title': plot_title,
        'plot_root': plot_root,
    }


def prep_log_nuc_freqs(
    df_nuc_freq,
    tot_aln_reads,
    include_idxs_list,
    plot_idxs,
    ref_len,
    sgRNA_legend,
    ref_name,
    ref_names,
    fig_filename_root,
    save_also_png,
):
    """Prepare kwargs for plot_log_nuc_freqs (plot_10d).

    Computes ``plot_quant_window_idxs`` — the indices in the plotting window
    that overlap with the quantification window.

    Parameters
    ----------
    plot_idxs : list of int
        Indices into the full-amplicon sequence selected for the plot window.
    """
    is_window = np.zeros(ref_len)
    for include_idx in include_idxs_list:
        is_window[include_idx] = 1

    plot_quant_window_idxs = []
    for plot_ind, loc in enumerate(plot_idxs):
        if is_window[loc]:
            plot_quant_window_idxs.append(plot_ind - 2)

    num_refs = len(ref_names)

    return {
        'df_nuc_freq': df_nuc_freq,
        'tot_aln_reads': tot_aln_reads,
        'plot_title': _plot_title_with_ref_name(
            'Log2 Nucleotide Frequencies Around the ' + sgRNA_legend,
            ref_name,
            num_refs,
        ),
        'fig_filename_root': fig_filename_root,
        'save_also_png': save_also_png,
        'quantification_window_idxs': plot_quant_window_idxs,
    }


def prep_conversion_at_sel_nucs(plot_nuc_pcts, plot_nuc_freqs, conversion_nuc_from):
    """Compute from_nuc_indices and selected-nucleotide DataFrames for plots 10e/10f/10g.

    Finds the column positions matching ``conversion_nuc_from`` and returns
    sliced, re-labeled DataFrames used for both CSV writes and plot inputs.

    Returns a dict with:
    - ``from_nuc_indices``: list of int column positions matching the conversion nucleotide
    - ``just_sel_nuc_pcts``: DataFrame sliced to those columns, columns renamed
    - ``just_sel_nuc_freqs``: same slicing applied to frequency DataFrame
    """
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
