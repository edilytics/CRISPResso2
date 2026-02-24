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


def _plot_title_with_ref_name(title, ref_name, num_refs):
    """Add ref_name suffix to plot title when there are multiple references."""
    if num_refs > 1:
        return title + ": " + ref_name
    return title


def _clip_to_percentile(values, counts, percentile=0.99):
    """Find the value at which cumulative counts reach the given percentile.

    Scans ``values`` and ``counts`` in order, returning the value where
    the running sum of counts first exceeds ``percentile * total``.
    If the cutoff is never reached, returns the last value.
    """
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
