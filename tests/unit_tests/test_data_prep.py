"""Tests for CRISPResso2.plots.data_prep — extracted data preparation functions.

Each prep function takes raw analysis data and returns the kwargs dict
that the corresponding CRISPRessoPlot function expects.
"""

import numpy as np
import pytest


class TestPrepReadBarplot:
    """prep_read_barplot (plot_1a) — trivial packaging."""

    def test_returns_correct_keys(self):
        from CRISPResso2.plots.data_prep import prep_read_barplot

        result = prep_read_barplot(
            N_READS_INPUT=1000,
            N_READS_AFTER_PREPROCESSING=900,
            N_TOTAL=800,
            fig_filename_root='/tmp/1a.Read_barplot',
            save_png=True,
        )
        assert result == {
            'N_READS_INPUT': 1000,
            'N_READS_AFTER_PREPROCESSING': 900,
            'N_TOTAL': 800,
            'fig_filename_root': '/tmp/1a.Read_barplot',
            'save_png': True,
        }

    def test_zero_reads(self):
        from CRISPResso2.plots.data_prep import prep_read_barplot

        result = prep_read_barplot(
            N_READS_INPUT=0,
            N_READS_AFTER_PREPROCESSING=0,
            N_TOTAL=0,
            fig_filename_root='/tmp/test',
            save_png=False,
        )
        assert result['N_READS_INPUT'] == 0
        assert result['N_TOTAL'] == 0


class TestPrepClassPiechartBarplot:
    """prep_class_piechart_barplot (plot_1b/1c) — trivial packaging."""

    def test_returns_correct_keys(self):
        from CRISPResso2.plots.data_prep import prep_class_piechart_barplot

        class_counts = {'Reference': 100, 'NHEJ': 50}
        result = prep_class_piechart_barplot(
            class_counts_order=['Reference', 'NHEJ'],
            class_counts=class_counts,
            ref_names=['ref1'],
            expected_hdr_amplicon_seq='',
            N_TOTAL=150,
            piechart_plot_root='/tmp/1b.pie',
            barplot_plot_root='/tmp/1c.bar',
            custom_colors={'ref1': '#ff0000'},
            save_png=True,
        )
        assert result['class_counts'] is class_counts
        assert result['N_TOTAL'] == 150
        assert result['piechart_plot_root'] == '/tmp/1b.pie'
        assert result['barplot_plot_root'] == '/tmp/1c.bar'
        assert result['save_png'] is True

    def test_all_keys_present(self):
        from CRISPResso2.plots.data_prep import prep_class_piechart_barplot

        result = prep_class_piechart_barplot(
            class_counts_order=[],
            class_counts={},
            ref_names=[],
            expected_hdr_amplicon_seq='ACGT',
            N_TOTAL=0,
            piechart_plot_root='/tmp/pie',
            barplot_plot_root='/tmp/bar',
            custom_colors={},
            save_png=False,
        )
        expected_keys = {
            'class_counts_order', 'class_counts', 'ref_names',
            'expected_hdr_amplicon_seq', 'N_TOTAL',
            'piechart_plot_root', 'barplot_plot_root',
            'custom_colors', 'save_png',
        }
        assert set(result.keys()) == expected_keys


class TestPrepAlleleHomology:
    """prep_allele_homology (plot_1e) — trivial packaging."""

    def test_returns_correct_keys(self):
        from CRISPResso2.plots.data_prep import prep_allele_homology

        scores = [85, 90, 95]
        counts = [10, 20, 30]
        result = prep_allele_homology(
            fig_root='/tmp/1e.Allele_homology',
            homology_scores=scores,
            counts=counts,
            min_homology=60,
            save_also_png=True,
        )
        assert result['fig_root'] == '/tmp/1e.Allele_homology'
        assert result['homology_scores'] is scores
        assert result['counts'] is counts
        assert result['min_homology'] == 60
        assert result['save_also_png'] is True


class TestPrepAmpliconModifications:
    """prep_amplicon_modifications (plot_4a) — packages vectors + metadata."""

    def test_returns_correct_keys(self):
        from CRISPResso2.plots.data_prep import prep_amplicon_modifications

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
        assert result['all_indelsub_count_vectors'] is indelsub
        assert result['n_total'] == 100
        assert result['n_this_category'] == 80
        assert result['ref_name'] == 'ref1'
        assert result['num_refs'] == 1
        assert result['ref_len'] == 5
        # y_max should be max(indelsub) * 1.1
        assert result['y_max'] == pytest.approx(4 * 1.1)
        assert result['plot_root'] == '/tmp/4a.Combined'

    def test_y_max_computed_from_vector(self):
        from CRISPResso2.plots.data_prep import prep_amplicon_modifications

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
        assert result['y_max'] == pytest.approx(10 * 1.1)

    def test_plot_titles_present(self):
        from CRISPResso2.plots.data_prep import prep_amplicon_modifications

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
        assert 'plot_titles' in result
        assert 'combined' in result['plot_titles']
        assert 'main' in result['plot_titles']


class TestPrepModificationFrequency:
    """prep_modification_frequency (plot_4b) — packages per-type vectors."""

    def test_returns_correct_keys(self):
        from CRISPResso2.plots.data_prep import prep_modification_frequency

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
        assert result['all_insertion_count_vectors'] is ins
        assert result['all_deletion_count_vectors'] is dels
        assert result['all_substitution_count_vectors'] is subs
        assert result['n_total'] == 50
        assert result['num_refs'] == 1
        assert result['y_max'] == 5.5
        assert result['plot_root'] == '/tmp/4b.Modification'

    def test_all_expected_keys(self):
        from CRISPResso2.plots.data_prep import prep_modification_frequency

        result = prep_modification_frequency(
            include_idxs_list=[],
            all_insertion_count_vector=np.array([]),
            all_deletion_count_vector=np.array([]),
            all_substitution_count_vector=np.array([]),
            sgRNA_intervals=[],
            ref_len=0,
            ref_name='r',
            ref_names=['r'],
            N_TOTAL=0,
            n_this_category=0,
            cut_points=[],
            plot_cut_points=[],
            y_max=0,
            plot_root='/tmp/t',
            custom_colors={},
            save_also_png=False,
        )
        expected_keys = {
            'include_idxs_list', 'all_insertion_count_vectors',
            'all_deletion_count_vectors', 'all_substitution_count_vectors',
            'sgRNA_intervals', 'ref_len', 'ref_name', 'num_refs',
            'n_total', 'n_this_category', 'cut_points', 'plot_cut_points',
            'y_max', 'plot_title', 'plot_root', 'custom_colors', 'save_also_png',
        }
        assert set(result.keys()) == expected_keys
