"""Tests for CRISPResso2.plots.data_prep — extracted data preparation functions.

Each prep function takes raw analysis data and returns the kwargs dict
that the corresponding CRISPRessoPlot function expects. Only functions
with non-trivial computation are extracted and tested here.
"""

import numpy as np
from inline_snapshot import snapshot

from CRISPResso2.plots.data_prep import (
    prep_amplicon_modifications,
    prep_frequency_deletions_insertions,
    prep_indel_size_distribution,
    prep_modification_frequency,
)


def _to_serializable(obj):
    """Recursively convert numpy types to plain Python for snapshotting."""
    if isinstance(obj, dict):
        return {k: _to_serializable(v) for k, v in sorted(obj.items())}
    if isinstance(obj, (list, tuple)):
        return type(obj)(_to_serializable(v) for v in obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    return obj


class TestPrepAmpliconModifications:
    """prep_amplicon_modifications (plot_4a) — packages vectors + metadata."""

    def test_single_ref(self):
        result = prep_amplicon_modifications(
            all_indelsub_count_vector=np.array([0, 1, 2, 3, 4]),
            include_idxs_list=[1, 2, 3],
            cut_points=[2],
            plot_cut_points=[2],
            sgRNA_intervals=[(1, 3)],
            N_TOTAL=100,
            n_this_category=80,
            ref_name='ref1',
            ref_names=['ref1'],
            ref_len=5,
            plot_root='/tmp/4a.Combined',
            custom_colors={},
            save_also_png=True,
        )
        assert _to_serializable(result) == snapshot({'all_indelsub_count_vectors': [0, 1, 2, 3, 4], 'custom_colors': {}, 'cut_points': [2], 'include_idxs_list': [1, 2, 3], 'n_this_category': 80, 'n_total': 100, 'num_refs': 1, 'plot_cut_points': [2], 'plot_root': '/tmp/4a.Combined', 'plot_titles': {'combined': 'Combined Insertions/Deletions/Substitutions', 'main': 'Mutation position distribution'}, 'ref_len': 5, 'ref_name': 'ref1', 'save_also_png': True, 'sgRNA_intervals': [(1, 3)], 'y_max': 4.4})

    def test_multi_ref(self):
        result = prep_amplicon_modifications(
            all_indelsub_count_vector=np.array([1]),
            include_idxs_list=[0],
            cut_points=[],
            plot_cut_points=[],
            sgRNA_intervals=[],
            N_TOTAL=10,
            n_this_category=10,
            ref_name='myref',
            ref_names=['myref', 'other'],
            ref_len=1,
            plot_root='/tmp/test',
            custom_colors={},
            save_also_png=False,
        )
        assert _to_serializable(result) == snapshot({'all_indelsub_count_vectors': [1], 'custom_colors': {}, 'cut_points': [], 'include_idxs_list': [0], 'n_this_category': 10, 'n_total': 10, 'num_refs': 2, 'plot_cut_points': [], 'plot_root': '/tmp/test', 'plot_titles': {'combined': 'Combined Insertions/Deletions/Substitutions: myref', 'main': 'Mutation position distribution: myref'}, 'ref_len': 1, 'ref_name': 'myref', 'save_also_png': False, 'sgRNA_intervals': [], 'y_max': 1.1})


class TestPrepModificationFrequency:
    """prep_modification_frequency (plot_4b) — packages per-type vectors."""

    def test_single_ref(self):
        result = prep_modification_frequency(
            include_idxs_list=[0, 1, 2],
            all_insertion_count_vector=np.array([0, 1, 0]),
            all_deletion_count_vector=np.array([1, 0, 0]),
            all_substitution_count_vector=np.array([0, 0, 1]),
            sgRNA_intervals=[(0, 2)],
            ref_len=3,
            ref_name='ref1',
            ref_names=['ref1'],
            N_TOTAL=50,
            n_this_category=40,
            cut_points=[1],
            plot_cut_points=[1],
            y_max=5.5,
            plot_root='/tmp/4b.Modification',
            custom_colors={},
            save_also_png=True,
        )
        assert _to_serializable(result) == snapshot({'all_deletion_count_vectors': [1, 0, 0], 'all_insertion_count_vectors': [0, 1, 0], 'all_substitution_count_vectors': [0, 0, 1], 'custom_colors': {}, 'cut_points': [1], 'include_idxs_list': [0, 1, 2], 'n_this_category': 40, 'n_total': 50, 'num_refs': 1, 'plot_cut_points': [1], 'plot_root': '/tmp/4b.Modification', 'plot_title': 'Mutation position distribution', 'ref_len': 3, 'ref_name': 'ref1', 'save_also_png': True, 'sgRNA_intervals': [(0, 2)], 'y_max': 5.5})

    def test_multi_ref(self):
        result = prep_modification_frequency(
            include_idxs_list=[],
            all_insertion_count_vector=np.array([]),
            all_deletion_count_vector=np.array([]),
            all_substitution_count_vector=np.array([]),
            sgRNA_intervals=[],
            ref_len=0,
            ref_name='FANC',
            ref_names=['FANC', 'HDR'],
            N_TOTAL=0,
            n_this_category=0,
            cut_points=[],
            plot_cut_points=[],
            y_max=0,
            plot_root='/tmp/t',
            custom_colors={},
            save_also_png=False,
        )
        assert _to_serializable(result) == snapshot({'all_deletion_count_vectors': [], 'all_insertion_count_vectors': [], 'all_substitution_count_vectors': [], 'custom_colors': {}, 'cut_points': [], 'include_idxs_list': [], 'n_this_category': 0, 'n_total': 0, 'num_refs': 2, 'plot_cut_points': [], 'plot_root': '/tmp/t', 'plot_title': 'Mutation position distribution: FANC', 'ref_len': 0, 'ref_name': 'FANC', 'save_also_png': False, 'sgRNA_intervals': [], 'y_max': 0})


class TestPrepIndelSizeDistribution:
    """prep_indel_size_distribution (plot_3a) — 99th percentile clipping."""

    def test_clipping_without_outliers(self):
        """99th percentile clips xmax/xmin, then clamped to ±15."""
        result = prep_indel_size_distribution(
            hdensity=np.array([1, 5, 94, 1]),
            hlengths=np.array([-30, -5, 0, 30]),
            center_index=2,
            n_this_category=101,
            ref_name='r',
            ref_names=['r'],
            plot_root='/tmp/t',
            save_also_png=False,
            plot_histogram_outliers=False,
        )
        assert _to_serializable(result) == snapshot({'center_index': 2, 'hdensity': [1, 5, 94, 1], 'hlengths': [-30, -5, 0, 30], 'n_this_category': 101, 'plot_root': '/tmp/t', 'ref_name': 'r', 'save_also_png': False, 'title': 'Indel size distribution', 'xmax': 15, 'xmin': -15})

    def test_no_clipping_with_outliers(self):
        """When plot_histogram_outliers=True, use raw min/max."""
        result = prep_indel_size_distribution(
            hdensity=np.array([1, 5, 90, 3, 1]),
            hlengths=np.array([-50, -10, 0, 10, 50]),
            center_index=2,
            n_this_category=100,
            ref_name='r',
            ref_names=['r'],
            plot_root='/tmp/t',
            save_also_png=False,
            plot_histogram_outliers=True,
        )
        assert _to_serializable(result) == snapshot({'center_index': 2, 'hdensity': [1, 5, 90, 3, 1], 'hlengths': [-50, -10, 0, 10, 50], 'n_this_category': 100, 'plot_root': '/tmp/t', 'ref_name': 'r', 'save_also_png': False, 'title': 'Indel size distribution', 'xmax': 50, 'xmin': -50})

    def test_xmin_xmax_clamped_to_15(self):
        """Even with narrow data, xmin/xmax should be at least ±15."""
        result = prep_indel_size_distribution(
            hdensity=np.array([50, 50]),
            hlengths=np.array([-2, 2]),
            center_index=0,
            n_this_category=100,
            ref_name='r',
            ref_names=['r'],
            plot_root='/tmp/t',
            save_also_png=False,
            plot_histogram_outliers=True,
        )
        assert _to_serializable(result) == snapshot({'center_index': 0, 'hdensity': [50, 50], 'hlengths': [-2, 2], 'n_this_category': 100, 'plot_root': '/tmp/t', 'ref_name': 'r', 'save_also_png': False, 'title': 'Indel size distribution', 'xmax': 15, 'xmin': -15})

    def test_multi_ref_title(self):
        result = prep_indel_size_distribution(
            hdensity=np.array([100]),
            hlengths=np.array([0]),
            center_index=0,
            n_this_category=100,
            ref_name='FANC',
            ref_names=['FANC', 'HDR'],
            plot_root='/tmp/t',
            save_also_png=False,
            plot_histogram_outliers=False,
        )
        assert _to_serializable(result) == snapshot({'center_index': 0, 'hdensity': [100], 'hlengths': [0], 'n_this_category': 100, 'plot_root': '/tmp/t', 'ref_name': 'FANC', 'save_also_png': False, 'title': 'Indel size distribution: FANC', 'xmax': 15, 'xmin': -15})


class TestPrepFrequencyDeletionsInsertions:
    """prep_frequency_deletions_insertions (plot_3b) — triple clipping."""

    def test_clipping_without_outliers(self):
        result = prep_frequency_deletions_insertions(
            x_bins_ins=np.array([0, 1, 100]),
            y_values_ins=np.array([95, 4, 1]),
            x_bins_del=np.array([0, 1, 200]),
            y_values_del=np.array([95, 4, 1]),
            x_bins_mut=np.array([0, 1, 300]),
            y_values_mut=np.array([95, 4, 1]),
            hdensity=np.array([100]),
            ref={},
            counts_total=100,
            ref_name='r',
            ref_names=['r'],
            plot_root='/tmp/t',
            save_also_png=False,
            custom_colors={},
            plot_histogram_outliers=False,
        )
        assert _to_serializable(result) == snapshot({'counts_total': 100, 'custom_colors': {}, 'plot_path': '/tmp/t', 'plot_titles': {'ins': 'Insertions', 'del': 'Deletions', 'mut': 'Substitutions'}, 'ref': {}, 'ref_name': 'r', 'save_also_png': False, 'xmax_del': 200, 'xmax_ins': 100, 'xmax_mut': 300})

    def test_no_clipping(self):
        result = prep_frequency_deletions_insertions(
            x_bins_ins=np.array([0, 50]),
            y_values_ins=np.array([50, 50]),
            x_bins_del=np.array([0, 60]),
            y_values_del=np.array([50, 50]),
            x_bins_mut=np.array([0, 70]),
            y_values_mut=np.array([50, 50]),
            hdensity=np.array([100]),
            ref={},
            counts_total=100,
            ref_name='r',
            ref_names=['r'],
            plot_root='/tmp/t',
            save_also_png=False,
            custom_colors={},
            plot_histogram_outliers=True,
        )
        assert _to_serializable(result) == snapshot({'counts_total': 100, 'custom_colors': {}, 'plot_path': '/tmp/t', 'plot_titles': {'ins': 'Insertions', 'del': 'Deletions', 'mut': 'Substitutions'}, 'ref': {}, 'ref_name': 'r', 'save_also_png': False, 'xmax_del': 60, 'xmax_ins': 50, 'xmax_mut': 70})

    def test_multi_ref_titles(self):
        result = prep_frequency_deletions_insertions(
            x_bins_ins=np.array([0]),
            y_values_ins=np.array([1]),
            x_bins_del=np.array([0]),
            y_values_del=np.array([1]),
            x_bins_mut=np.array([0]),
            y_values_mut=np.array([1]),
            hdensity=np.array([1]),
            ref={},
            counts_total=1,
            ref_name='FANC',
            ref_names=['FANC', 'HDR'],
            plot_root='/tmp/t',
            save_also_png=False,
            custom_colors={},
            plot_histogram_outliers=False,
        )
        assert _to_serializable(result) == snapshot({'counts_total': 1, 'custom_colors': {}, 'plot_path': '/tmp/t', 'plot_titles': {'ins': 'Insertions: FANC', 'del': 'Deletions: FANC', 'mut': 'Substitutions: FANC'}, 'ref': {}, 'ref_name': 'FANC', 'save_also_png': False, 'xmax_del': 15, 'xmax_ins': 15, 'xmax_mut': 15})
