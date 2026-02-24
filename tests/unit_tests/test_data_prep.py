"""Tests for CRISPResso2.plots.data_prep — extracted data preparation functions.

Each prep function takes raw analysis data and returns the kwargs dict
that the corresponding CRISPRessoPlot function expects. Only functions
with non-trivial computation are extracted and tested here.
"""

import numpy as np
import pytest
from inline_snapshot import snapshot

from CRISPResso2.plots.data_prep import (
    prep_amplicon_modifications,
    prep_frequency_deletions_insertions,
    prep_indel_size_distribution,
    prep_modification_frequency,
)


class TestPrepAmpliconModifications:
    """prep_amplicon_modifications (plot_4a) — packages vectors + metadata."""

    def test_output_shape(self):
        indelsub = np.array([0, 1, 2, 3, 4])
        result = prep_amplicon_modifications(
            all_indelsub_count_vector=indelsub,
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
        # Vector is passed through by reference
        assert result['all_indelsub_count_vectors'] is indelsub
        # Computed values
        assert result['y_max'] == snapshot(4.4)
        assert result['num_refs'] == snapshot(1)
        assert sorted(result.keys()) == snapshot(
            [
                "all_indelsub_count_vectors",
                "custom_colors",
                "cut_points",
                "include_idxs_list",
                "n_this_category",
                "n_total",
                "num_refs",
                "plot_cut_points",
                "plot_root",
                "plot_titles",
                "ref_len",
                "ref_name",
                "save_also_png",
                "sgRNA_intervals",
                "y_max",
            ]
        )

    def test_y_max_computed_from_vector(self):
        indelsub = np.array([10, 0, 5])
        result = prep_amplicon_modifications(
            all_indelsub_count_vector=indelsub,
            include_idxs_list=[0],
            cut_points=[],
            plot_cut_points=[],
            sgRNA_intervals=[],
            N_TOTAL=50,
            n_this_category=50,
            ref_name='amp',
            ref_names=['amp'],
            ref_len=3,
            plot_root='/tmp/test',
            custom_colors={},
            save_also_png=False,
        )
        assert result['y_max'] == snapshot(11.0)

    def test_plot_titles_single_ref(self):
        result = prep_amplicon_modifications(
            all_indelsub_count_vector=np.array([1]),
            include_idxs_list=[0],
            cut_points=[],
            plot_cut_points=[],
            sgRNA_intervals=[],
            N_TOTAL=10,
            n_this_category=10,
            ref_name='myref',
            ref_names=['myref'],
            ref_len=1,
            plot_root='/tmp/test',
            custom_colors={},
            save_also_png=False,
        )
        assert result['plot_titles'] == snapshot(
            {
                "combined": "Combined Insertions/Deletions/Substitutions",
                "main": "Mutation position distribution",
            }
        )

    def test_plot_titles_multi_ref(self):
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
        assert result['plot_titles'] == snapshot(
            {
                "combined": "Combined Insertions/Deletions/Substitutions: myref",
                "main": "Mutation position distribution: myref",
            }
        )


class TestPrepModificationFrequency:
    """prep_modification_frequency (plot_4b) — packages per-type vectors."""

    def test_output_shape(self):
        ins = np.array([0, 1, 0])
        dels = np.array([1, 0, 0])
        subs = np.array([0, 0, 1])
        result = prep_modification_frequency(
            include_idxs_list=[0, 1, 2],
            all_insertion_count_vector=ins,
            all_deletion_count_vector=dels,
            all_substitution_count_vector=subs,
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
        # Vectors passed through by reference
        assert result['all_insertion_count_vectors'] is ins
        assert result['all_deletion_count_vectors'] is dels
        assert result['all_substitution_count_vectors'] is subs
        # Computed values
        assert result['num_refs'] == snapshot(1)
        assert result['plot_title'] == snapshot("Mutation position distribution")
        assert sorted(result.keys()) == snapshot(
            [
                "all_deletion_count_vectors",
                "all_insertion_count_vectors",
                "all_substitution_count_vectors",
                "custom_colors",
                "cut_points",
                "include_idxs_list",
                "n_this_category",
                "n_total",
                "num_refs",
                "plot_cut_points",
                "plot_root",
                "plot_title",
                "ref_len",
                "ref_name",
                "save_also_png",
                "sgRNA_intervals",
                "y_max",
            ]
        )

    def test_plot_title_multi_ref(self):
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
        assert result['plot_title'] == snapshot(
            "Mutation position distribution: FANC"
        )


class TestPrepIndelSizeDistribution:
    """prep_indel_size_distribution (plot_3a) — 99th percentile clipping."""

    def test_output_shape(self):
        hdensity = np.array([5, 10, 80, 3, 2])
        hlengths = np.array([-20, -5, 0, 5, 50])
        result = prep_indel_size_distribution(
            hdensity=hdensity,
            hlengths=hlengths,
            center_index=2,
            n_this_category=100,
            ref_name='ref1',
            ref_names=['ref1'],
            plot_root='/tmp/3a.Indel',
            save_also_png=True,
            plot_histogram_outliers=False,
        )
        assert result['hdensity'] is hdensity
        assert result['hlengths'] is hlengths
        assert sorted(result.keys()) == snapshot(['center_index', 'hdensity', 'hlengths', 'n_this_category', 'plot_root', 'ref_name', 'save_also_png', 'title', 'xmax', 'xmin'])

    def test_clipping_without_outliers(self):
        """99th percentile should clip xmax and xmin."""
        # 100 counts: 98 at 0, 1 at -30, 1 at +30
        # 99% cutoff = 99. Scanning from left: at 0 we hit 99 (5+94), so xmax=0
        # But clamped to at least 15
        hdensity = np.array([1, 5, 94, 1])
        hlengths = np.array([-30, -5, 0, 30])
        result = prep_indel_size_distribution(
            hdensity=hdensity,
            hlengths=hlengths,
            center_index=2,
            n_this_category=101,
            ref_name='r',
            ref_names=['r'],
            plot_root='/tmp/t',
            save_also_png=False,
            plot_histogram_outliers=False,
        )
        # xmax clipped from 30 to something <= 15 (clamped), xmin similarly
        assert result['xmin'] == snapshot(-15)
        assert result['xmax'] == snapshot(15)

    def test_no_clipping_with_outliers(self):
        """When plot_histogram_outliers=True, use raw min/max."""
        hdensity = np.array([1, 5, 90, 3, 1])
        hlengths = np.array([-50, -10, 0, 10, 50])
        result = prep_indel_size_distribution(
            hdensity=hdensity,
            hlengths=hlengths,
            center_index=2,
            n_this_category=100,
            ref_name='r',
            ref_names=['r'],
            plot_root='/tmp/t',
            save_also_png=False,
            plot_histogram_outliers=True,
        )
        assert result['xmin'] == snapshot(-50)
        assert result['xmax'] == snapshot(50)

    def test_xmin_xmax_clamped_to_15(self):
        """Even with narrow data, xmin/xmax should be at least ±15."""
        hdensity = np.array([50, 50])
        hlengths = np.array([-2, 2])
        result = prep_indel_size_distribution(
            hdensity=hdensity,
            hlengths=hlengths,
            center_index=0,
            n_this_category=100,
            ref_name='r',
            ref_names=['r'],
            plot_root='/tmp/t',
            save_also_png=False,
            plot_histogram_outliers=True,
        )
        assert result['xmin'] == snapshot(-15)
        assert result['xmax'] == snapshot(15)

    def test_title_with_multi_ref(self):
        hdensity = np.array([100])
        hlengths = np.array([0])
        result = prep_indel_size_distribution(
            hdensity=hdensity,
            hlengths=hlengths,
            center_index=0,
            n_this_category=100,
            ref_name='FANC',
            ref_names=['FANC', 'HDR'],
            plot_root='/tmp/t',
            save_also_png=False,
            plot_histogram_outliers=False,
        )
        assert result['title'] == snapshot('Indel size distribution: FANC')


class TestPrepFrequencyDeletionsInsertions:
    """prep_frequency_deletions_insertions (plot_3b) — triple clipping."""

    def test_output_shape(self):
        result = prep_frequency_deletions_insertions(
            x_bins_ins=np.array([0, 1, 2, 3, 50]),
            y_values_ins=np.array([80, 10, 5, 3, 2]),
            x_bins_del=np.array([0, 1, 2, 3, 40]),
            y_values_del=np.array([75, 15, 5, 3, 2]),
            x_bins_mut=np.array([0, 1, 2, 3, 30]),
            y_values_mut=np.array([70, 20, 5, 3, 2]),
            hdensity=np.array([100]),
            ref={"some": "data"},
            counts_total=100,
            ref_name='ref1',
            ref_names=['ref1'],
            plot_root='/tmp/3b.Hist',
            save_also_png=True,
            custom_colors={},
            plot_histogram_outliers=False,
        )
        assert sorted(result.keys()) == snapshot(['counts_total', 'custom_colors', 'plot_path', 'plot_titles', 'ref', 'ref_name', 'save_also_png', 'xmax_del', 'xmax_ins', 'xmax_mut'])

    def test_clipping_without_outliers(self):
        """Each xmax should be clipped to 99th percentile, clamped to 15."""
        result = prep_frequency_deletions_insertions(
            x_bins_ins=np.array([0, 1, 100]),
            y_values_ins=np.array([95, 4, 1]),
            x_bins_del=np.array([0, 1, 200]),
            y_values_del=np.array([95, 4, 1]),
            x_bins_mut=np.array([0, 1, 300]),
            y_values_mut=np.array([95, 4, 1]),
            hdensity=np.array([100]),  # sum=100, cutoff=99
            ref={},
            counts_total=100,
            ref_name='r',
            ref_names=['r'],
            plot_root='/tmp/t',
            save_also_png=False,
            custom_colors={},
            plot_histogram_outliers=False,
        )
        # 95+4=99 → cutoff at index 1, but clamped to 15
        assert result['xmax_ins'] == snapshot(100)
        assert result['xmax_del'] == snapshot(200)
        assert result['xmax_mut'] == snapshot(300)

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
        assert result['xmax_ins'] == snapshot(50)
        assert result['xmax_del'] == snapshot(60)
        assert result['xmax_mut'] == snapshot(70)

    def test_plot_titles_multi_ref(self):
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
        assert result['plot_titles'] == snapshot({'ins': 'Insertions: FANC', 'del': 'Deletions: FANC', 'mut': 'Substitutions: FANC'})
