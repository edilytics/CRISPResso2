"""Tests for CRISPResso2.plots.data_prep — extracted data preparation functions.

Each prep function takes a PlotContext as its only argument and returns
the kwargs dict that the corresponding CRISPRessoPlot function expects.
"""

import numpy as np
import pandas as pd
from collections import Counter
from types import SimpleNamespace

from CRISPResso2.plots.data_prep import (
    _prep_windowed_alleles,
    amino_acids_to_numbers,
    prep_alleles_around_cut,
    prep_alleles_table,
    prep_alleles_table_compare,
    prep_alternate_allele_counts,
    prep_amino_acid_table,
    prep_amino_acid_table_for_plot,
    prep_amplicon_modifications,
    prep_class_piechart_and_barplot,
    prep_base_edit_quilt,
    prep_conversion_at_sel_nucs,
    prep_dsODN_piechart,
    prep_frequency_deletions_insertions,
    prep_global_frameshift_data,
    prep_global_modifications_reference,
    prep_hdr_nucleotide_quilt,
    prep_indel_size_distribution,
    prep_modification_frequency,
    prep_nucleotide_quilt,
    prep_pe_nucleotide_quilt,
    prep_pe_nucleotide_quilt_around_sgRNA,
    _to_numeric_ignore_columns,
)
from CRISPResso2.plots.plot_context import PlotContext


# =============================================================================
# Helper: Build minimal PlotContext for tests
# =============================================================================


def _make_ctx(**overrides):
    """Build a minimal PlotContext, merging *overrides* into defaults."""
    defaults = dict(
        args=SimpleNamespace(
            plot_histogram_outliers=False,
            plot_window_size=20,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
            conversion_nuc_from='C',
            expected_hdr_amplicon_seq='',
            base_editor_output=False,
            coding_seq='',
        ),
        run_data={'running_info': {}},
        refs={},
        ref_names=[],
        counts_total={},
        counts_modified={},
        counts_unmodified={},
        counts_discarded={},
        counts_insertion={},
        counts_deletion={},
        counts_substitution={},
        class_counts={},
        N_TOTAL=0,
        df_alleles=None,
        all_insertion_count_vectors={},
        all_insertion_left_count_vectors={},
        all_deletion_count_vectors={},
        all_substitution_count_vectors={},
        all_indelsub_count_vectors={},
        all_substitution_base_vectors={},
        all_base_count_vectors={},
        insertion_count_vectors={},
        deletion_count_vectors={},
        substitution_count_vectors={},
        insertion_length_vectors={},
        deletion_length_vectors={},
        hists_frameshift={},
        hists_inframe={},
        counts_modified_frameshift={},
        counts_modified_non_frameshift={},
        counts_non_modified_non_frameshift={},
        counts_splicing_sites_modified={},
    )
    defaults.update(overrides)
    return PlotContext(**defaults)


def _ref_dict(**overrides):
    """Build a minimal per-reference dict, merging *overrides* into defaults."""
    defaults = dict(
        sequence='ACGT',
        sequence_length=4,
        ref_plot_name='',
        include_idxs=[0, 1, 2, 3],
        sgRNA_cut_points=[2],
        sgRNA_plot_cut_points=[2],
        sgRNA_intervals=[(1, 3)],
        sgRNA_names=[''],
        sgRNA_mismatches=[],
        sgRNA_sequences=['ACGT'],
        sgRNA_orig_sequences=['ACGT'],
        sgRNA_plot_idxs=[[0, 1, 2, 3]],
        hdensity=np.array([100]),
        hlengths=np.array([0]),
        center_index=0,
        x_bins_ins=np.array([0, 1]),
        y_values_ins=np.array([90, 10]),
        x_bins_del=np.array([0, 1]),
        y_values_del=np.array([80, 20]),
        x_bins_mut=np.array([0, 1]),
        y_values_mut=np.array([85, 15]),
        contains_coding_seq=False,
        exon_intervals=[],
        exon_positions=[],
        exon_len_mods=[],
    )
    defaults.update(overrides)
    return defaults


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


# =============================================================================
# Tests for _to_numeric_ignore_columns
# =============================================================================


def test_to_numeric_ignore_columns_basic():
    """Test _to_numeric_ignore_columns with basic dataframe."""
    df = pd.DataFrame({
        "name": ["a", "b", "c"],
        "value": ["1", "2", "3"],
        "count": ["10", "20", "30"]
    })
    result = _to_numeric_ignore_columns(df, {"name"})
    assert result["value"].dtype in [int, float, "int64", "float64"]
    assert result["count"].dtype in [int, float, "int64", "float64"]
    assert result["name"].dtype == object


def test_to_numeric_ignore_columns_multiple_ignore():
    """Test _to_numeric_ignore_columns ignoring multiple columns."""
    df = pd.DataFrame({
        "name": ["a", "b"],
        "label": ["x", "y"],
        "value": ["1", "2"]
    })
    result = _to_numeric_ignore_columns(df, {"name", "label"})
    assert result["name"].dtype == object
    assert result["label"].dtype == object
    assert result["value"].dtype in [int, float, "int64", "float64"]


def test_to_numeric_ignore_columns_empty_ignore():
    """Test _to_numeric_ignore_columns with empty ignore set."""
    df = pd.DataFrame({
        "a": ["1", "2"],
        "b": ["3", "4"]
    })
    result = _to_numeric_ignore_columns(df, set())
    assert result["a"].dtype in [int, float, "int64", "float64"]
    assert result["b"].dtype in [int, float, "int64", "float64"]


def test_to_numeric_ignore_columns_all_numeric():
    """Test _to_numeric_ignore_columns with all numeric columns."""
    df = pd.DataFrame({
        "a": ["1", "2", "3"],
        "b": ["4.5", "5.5", "6.5"]
    })

    result = _to_numeric_ignore_columns(df, set())

    assert result["a"].dtype in [int, float, "int64", "float64"]
    assert result["b"].dtype in [int, float, "int64", "float64"]


def test_to_numeric_ignore_columns_preserve_strings():
    """Test _to_numeric_ignore_columns preserves string columns in ignore list."""
    df = pd.DataFrame({
        "name": ["abc", "def", "ghi"],
        "value": ["1", "2", "3"]
    })

    result = _to_numeric_ignore_columns(df, {"name"})

    assert result["name"].dtype == object
    assert list(result["name"]) == ["abc", "def", "ghi"]


def test_to_numeric_ignore_columns_float_strings():
    """Test _to_numeric_ignore_columns with float strings."""
    df = pd.DataFrame({
        "val": ["1.5", "2.5", "3.5"]
    })

    result = _to_numeric_ignore_columns(df, set())

    assert result["val"].dtype == float


def test_to_numeric_ignore_columns_mixed_types():
    """Test _to_numeric_ignore_columns with mixed data types."""
    df = pd.DataFrame({
        "name": ["a", "b", "c"],
        "int_val": ["1", "2", "3"],
        "float_val": ["1.1", "2.2", "3.3"]
    })

    result = _to_numeric_ignore_columns(df, {"name"})

    assert result["name"].dtype == object
    assert result["int_val"].dtype in [int, float, "int64", "float64"]
    assert result["float_val"].dtype in [float, "float64"]


# =============================================================================
# Tests: prep_amplicon_modifications (plot 4a)
# =============================================================================


class TestPrepAmpliconModifications:

    def test_single_ref(self):
        ctx = _make_ctx(
            ref_names=['ref1'],
            refs={'ref1': _ref_dict(
                sequence='ACGTX',
                sequence_length=5,
                include_idxs=[1, 2, 3],
                sgRNA_cut_points=[2],
                sgRNA_plot_cut_points=[2],
                sgRNA_intervals=[(1, 3)],
            )},
            all_indelsub_count_vectors={'ref1': np.array([0, 1, 2, 3, 4])},
            counts_total={'ref1': 80},
            N_TOTAL=100,
        )
        ctx.ref_name = 'ref1'
        result = prep_amplicon_modifications(ctx)

        assert result['y_max'] == 4 * 1.1
        assert result['n_total'] == 100
        assert result['n_this_category'] == 80
        assert result['ref_len'] == 5
        assert result['plot_titles']['main'] == 'Mutation position distribution'

    def test_multi_ref_title(self):
        ctx = _make_ctx(
            ref_names=['myref', 'other'],
            refs={
                'myref': _ref_dict(),
                'other': _ref_dict(),
            },
            all_indelsub_count_vectors={'myref': np.array([1])},
            counts_total={'myref': 10},
            N_TOTAL=10,
        )
        ctx.ref_name = 'myref'
        result = prep_amplicon_modifications(ctx)
        assert 'myref' in result['plot_titles']['main']


# =============================================================================
# Tests: prep_modification_frequency (plot 4b)
# =============================================================================


class TestPrepModificationFrequency:

    def test_y_max_from_indelsub(self):
        """y_max is computed internally from all_indelsub_count_vectors."""
        ctx = _make_ctx(
            ref_names=['ref1'],
            refs={'ref1': _ref_dict(
                sequence='ACG',
                sequence_length=3,
                include_idxs=[0, 1, 2],
                sgRNA_intervals=[(0, 2)],
            )},
            all_insertion_count_vectors={'ref1': np.array([0, 1, 0])},
            all_deletion_count_vectors={'ref1': np.array([1, 0, 0])},
            all_substitution_count_vectors={'ref1': np.array([0, 0, 1])},
            all_indelsub_count_vectors={'ref1': np.array([1, 1, 1])},
            counts_total={'ref1': 40},
            N_TOTAL=50,
        )
        ctx.ref_name = 'ref1'
        result = prep_modification_frequency(ctx)
        assert result['y_max'] == 1 * 1.1
        assert result['n_total'] == 50
        assert result['n_this_category'] == 40


# =============================================================================
# Tests: prep_indel_size_distribution (plot 3a)
# =============================================================================


class TestPrepIndelSizeDistribution:

    def test_clipping_without_outliers(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                hdensity=np.array([1, 5, 94, 1]),
                hlengths=np.array([-30, -5, 0, 30]),
                center_index=2,
            )},
            counts_total={'r': 101},
            args=SimpleNamespace(plot_histogram_outliers=False),
        )
        ctx.ref_name = 'r'
        result = prep_indel_size_distribution(ctx)
        assert result['xmin'] == -15
        assert result['xmax'] == 15
        assert result['title'] == 'Indel size distribution'

    def test_no_clipping_with_outliers(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                hdensity=np.array([1, 5, 90, 3, 1]),
                hlengths=np.array([-50, -10, 0, 10, 50]),
                center_index=2,
            )},
            counts_total={'r': 100},
            args=SimpleNamespace(plot_histogram_outliers=True),
        )
        ctx.ref_name = 'r'
        result = prep_indel_size_distribution(ctx)
        assert result['xmin'] == -50
        assert result['xmax'] == 50

    def test_xmin_xmax_clamped_to_15(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                hdensity=np.array([50, 50]),
                hlengths=np.array([-2, 2]),
                center_index=0,
            )},
            counts_total={'r': 100},
            args=SimpleNamespace(plot_histogram_outliers=True),
        )
        ctx.ref_name = 'r'
        result = prep_indel_size_distribution(ctx)
        assert result['xmin'] == -15
        assert result['xmax'] == 15

    def test_multi_ref_title(self):
        ctx = _make_ctx(
            ref_names=['FANC', 'HDR'],
            refs={
                'FANC': _ref_dict(hdensity=np.array([100]), hlengths=np.array([0])),
                'HDR': _ref_dict(),
            },
            counts_total={'FANC': 100},
            args=SimpleNamespace(plot_histogram_outliers=False),
        )
        ctx.ref_name = 'FANC'
        result = prep_indel_size_distribution(ctx)
        assert result['title'] == 'Indel size distribution: FANC'


# =============================================================================
# Tests: prep_frequency_deletions_insertions (plot 3b)
# =============================================================================


class TestPrepFrequencyDeletionsInsertions:

    def test_clipping(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                x_bins_ins=np.arange(101),
                y_values_ins=np.zeros(101),
                x_bins_del=np.arange(101),
                y_values_del=np.zeros(101),
                x_bins_mut=np.arange(101),
                y_values_mut=np.zeros(101),
                hdensity=np.zeros(100),
            )},
            counts_total={'r': 100},
            args=SimpleNamespace(plot_histogram_outliers=True),
        )
        ctx.ref_name = 'r'
        result = prep_frequency_deletions_insertions(ctx)
        assert result['xmax_ins'] == 100
        assert result['xmax_del'] == 100

    def test_xmax_clamped_to_15(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                x_bins_ins=np.arange(5),
                y_values_ins=np.array([10, 5, 3, 1, 1]),
                x_bins_del=np.arange(5),
                y_values_del=np.array([10, 5, 3, 1, 1]),
                x_bins_mut=np.arange(5),
                y_values_mut=np.array([10, 5, 3, 1, 1]),
                hdensity=np.ones(20),
            )},
            counts_total={'r': 100},
            args=SimpleNamespace(plot_histogram_outliers=False),
        )
        ctx.ref_name = 'r'
        result = prep_frequency_deletions_insertions(ctx)
        assert result['xmax_ins'] >= 15
        assert result['xmax_del'] >= 15


# =============================================================================
# Tests: prep_dsODN_piechart (plot 1d)
# =============================================================================


class TestPrepDsODNPiechart:

    def test_basic(self):
        df = pd.DataFrame({
            '#Reads': [60, 40],
            'contains dsODN': [True, False],
        })
        ctx = _make_ctx(df_alleles=df, N_TOTAL=100)
        result = prep_dsODN_piechart(ctx)
        assert len(result['labels']) == 2
        assert abs(result['sizes'][0] - 60.0) < 0.01
        assert abs(result['sizes'][1] - 40.0) < 0.01


# =============================================================================
# Tests: prep_nucleotide_quilt (plot 2a)
# =============================================================================


class TestPrepNucleotideQuilt:

    def test_basic(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence='AC',
                sequence_length=2,
                include_idxs=[0, 1],
            )},
            all_base_count_vectors={
                'r_A': np.array([80, 10]),
                'r_C': np.array([10, 80]),
                'r_G': np.array([5, 5]),
                'r_T': np.array([5, 5]),
                'r_N': np.array([0, 0]),
                'r_-': np.array([0, 0]),
            },
            all_insertion_count_vectors={'r': np.array([2, 3])},
            all_insertion_left_count_vectors={'r': np.array([1, 2])},
            all_deletion_count_vectors={'r': np.array([1, 1])},
            all_substitution_count_vectors={'r': np.array([3, 3])},
            all_indelsub_count_vectors={'r': np.array([5, 5])},
            counts_total={'r': 100},
        )
        ctx.ref_name = 'r'
        result = prep_nucleotide_quilt(ctx)

        assert 'nuc_pct_df' in result
        assert 'mod_pct_df' in result
        assert result['nuc_pct_df'].shape[0] == 6  # A,C,G,T,N,-
        assert result['mod_pct_df'].shape[0] == 6  # Ins, InsL, Del, Sub, All, Total


# =============================================================================
# Tests: prep_global_frameshift_data
# =============================================================================


class TestPrepGlobalFrameshiftData:

    def test_aggregation(self):
        ctx = _make_ctx(
            ref_names=['ref1', 'ref2'],
            refs={
                'ref1': _ref_dict(contains_coding_seq=True),
                'ref2': _ref_dict(contains_coding_seq=False),
            },
            counts_modified_frameshift={'ref1': 10, 'ref2': 5},
            counts_modified_non_frameshift={'ref1': 20, 'ref2': 3},
            counts_non_modified_non_frameshift={'ref1': 60, 'ref2': 80},
            counts_splicing_sites_modified={'ref1': 2, 'ref2': 1},
            counts_total={'ref1': 100, 'ref2': 100},
            counts_modified={'ref1': 30, 'ref2': 8},
            counts_unmodified={'ref1': 70, 'ref2': 92},
            hists_frameshift={'ref1': Counter({0: 5, 3: 10}), 'ref2': Counter({0: 80})},
            hists_inframe={'ref1': Counter({0: 60, 3: 20}), 'ref2': Counter({0: 92})},
        )
        result = prep_global_frameshift_data(ctx)
        # Only ref1 has coding seq
        assert result['global_modified_frameshift'] == 10
        assert result['global_modified_non_frameshift'] == 20
        assert result['global_count_total'] == 100

    def test_hdr_adds_unmodified(self):
        ctx = _make_ctx(
            ref_names=['HDR'],
            refs={'HDR': _ref_dict(contains_coding_seq=True)},
            counts_modified_frameshift={'HDR': 5},
            counts_modified_non_frameshift={'HDR': 10},
            counts_non_modified_non_frameshift={'HDR': 20},
            counts_splicing_sites_modified={'HDR': 1},
            counts_total={'HDR': 100},
            counts_modified={'HDR': 15},
            counts_unmodified={'HDR': 85},
            hists_frameshift={'HDR': Counter({0: 10})},
            hists_inframe={'HDR': Counter({0: 20})},
        )
        result = prep_global_frameshift_data(ctx)
        # HDR unmodified reads added to non_modified_non_frameshift
        assert result['global_non_modified_non_frameshift'] == 20 + 85


# =============================================================================
# Tests: prep_hdr_nucleotide_quilt (plot 4g)
# =============================================================================


class TestPrepHdrNucleotideQuilt:

    def test_basic(self):
        ctx = _make_ctx(
            ref_names=['FANC', 'HDR'],
            refs={
                'FANC': _ref_dict(sequence='AC', sequence_length=2),
                'HDR': _ref_dict(sequence='AC', sequence_length=2),
            },
            counts_total={'FANC': 100, 'HDR': 50},
            ref1_all_base_count_vectors={
                'FANC_A': np.array([80, 10]),
                'FANC_C': np.array([10, 80]),
                'FANC_G': np.array([5, 5]),
                'FANC_T': np.array([5, 5]),
                'FANC_N': np.array([0, 0]),
                'FANC_-': np.array([0, 0]),
                'HDR_A': np.array([40, 5]),
                'HDR_C': np.array([5, 40]),
                'HDR_G': np.array([2, 2]),
                'HDR_T': np.array([3, 3]),
                'HDR_N': np.array([0, 0]),
                'HDR_-': np.array([0, 0]),
            },
            ref1_all_insertion_count_vectors={'FANC': np.array([1, 2]), 'HDR': np.array([0, 1])},
            ref1_all_insertion_left_count_vectors={'FANC': np.array([0, 1]), 'HDR': np.array([0, 0])},
            ref1_all_deletion_count_vectors={'FANC': np.array([1, 0]), 'HDR': np.array([0, 0])},
            ref1_all_substitution_count_vectors={'FANC': np.array([2, 1]), 'HDR': np.array([1, 0])},
            ref1_all_indelsub_count_vectors={'FANC': np.array([3, 3]), 'HDR': np.array([1, 1])},
        )
        result = prep_hdr_nucleotide_quilt(ctx)
        assert 'nuc_pct_df' in result
        assert 'mod_pct_df' in result
        assert result['quantification_window_idxs'] == []  # HDR default


# =============================================================================
# Tests: prep_pe_nucleotide_quilt (plot 11a)
# =============================================================================


class TestPrepPeNucleotideQuilt:

    def test_includes_quantification_window(self):
        ctx = _make_ctx(
            ref_names=['WT', 'PE'],
            refs={
                'WT': _ref_dict(sequence='AC', sequence_length=2, include_idxs=[0, 1]),
                'PE': _ref_dict(sequence='AC', sequence_length=2),
            },
            counts_total={'WT': 100, 'PE': 50},
            ref1_all_base_count_vectors={
                'WT_A': np.array([80, 10]), 'WT_C': np.array([10, 80]),
                'WT_G': np.array([5, 5]), 'WT_T': np.array([5, 5]),
                'WT_N': np.array([0, 0]), 'WT_-': np.array([0, 0]),
                'PE_A': np.array([40, 5]), 'PE_C': np.array([5, 40]),
                'PE_G': np.array([2, 2]), 'PE_T': np.array([3, 3]),
                'PE_N': np.array([0, 0]), 'PE_-': np.array([0, 0]),
            },
            ref1_all_insertion_count_vectors={'WT': np.array([1, 2]), 'PE': np.array([0, 1])},
            ref1_all_insertion_left_count_vectors={'WT': np.array([0, 1]), 'PE': np.array([0, 0])},
            ref1_all_deletion_count_vectors={'WT': np.array([1, 0]), 'PE': np.array([0, 0])},
            ref1_all_substitution_count_vectors={'WT': np.array([2, 1]), 'PE': np.array([1, 0])},
            ref1_all_indelsub_count_vectors={'WT': np.array([3, 3]), 'PE': np.array([1, 1])},
        )
        result = prep_pe_nucleotide_quilt(ctx)
        # PE uses first ref's include_idxs
        assert result['quantification_window_idxs'] == [0, 1]


# =============================================================================
# Tests: _prep_windowed_alleles (private helper, unchanged)
# =============================================================================


class TestPrepWindowedAlleles:

    def test_basic_windowing(self):
        df = pd.DataFrame({
            'Aligned_Sequence': ['ACG', 'AXG'],
            'Reference_Sequence': ['ACG', 'ACG'],
            '#Reads': [80, 20],
            '%Reads': [80.0, 20.0],
        }).set_index('Aligned_Sequence')

        df_out, df_plot, ref_seq, intervals, start = _prep_windowed_alleles(
            df_alleles_around_cut=df,
            cut_point=1,
            window_left=1,
            window_right=1,
            ref_sequence='AACGG',
            sgRNA_intervals=[(0, 4)],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
        )
        assert ref_seq == 'AC'
        assert len(intervals) == 1


# =============================================================================
# Tests: prep_class_piechart_and_barplot (utility, not PlotContext-based)
# =============================================================================


class TestPrepClassPiechartAndBarplot:

    def test_basic_single_ref(self):
        result = prep_class_piechart_and_barplot(
            class_counts_order=['Reference_UNMODIFIED', 'Reference_MODIFIED'],
            class_counts={'Reference_UNMODIFIED': 80, 'Reference_MODIFIED': 20},
            ref_names=['Reference'],
            expected_hdr_amplicon_seq='',
            N_TOTAL=100,
        )
        assert 'UNMODIFIED' in result['labels'][0]
        assert abs(result['sizes'][0] - 80.0) < 0.01

    def test_hdr_mode_labels(self):
        result = prep_class_piechart_and_barplot(
            class_counts_order=['FANC_MODIFIED', 'HDR_MODIFIED', 'HDR_UNMODIFIED'],
            class_counts={'FANC_MODIFIED': 10, 'HDR_MODIFIED': 5, 'HDR_UNMODIFIED': 85},
            ref_names=['FANC', 'HDR'],
            expected_hdr_amplicon_seq='ACGT',
            N_TOTAL=100,
        )
        assert 'NHEJ' in result['labels'][0]
        assert 'Imperfect HDR' in result['labels'][1]
        assert 'HDR' in result['labels'][2]

    def test_return_keys(self):
        result = prep_class_piechart_and_barplot(
            class_counts_order=['R_UNMODIFIED'],
            class_counts={'R_UNMODIFIED': 50},
            ref_names=['R'],
            expected_hdr_amplicon_seq='',
            N_TOTAL=50,
        )
        assert sorted(result.keys()) == ['N_TOTAL', 'labels', 'sizes']


# =============================================================================
# Tests: prep_alternate_allele_counts (utility, not PlotContext-based)
# =============================================================================


class TestPrepAlternateAlleleCounts:

    def test_basic(self):
        sub_base_vectors = {
            'R_A': np.array([0, 0, 0]),
            'R_C': np.array([0, 90, 0]),
            'R_G': np.array([0, 0, 95]),
            'R_T': np.array([0, 5, 0]),
            'R_N': np.array([0, 0, 0]),
        }
        result = prep_alternate_allele_counts(sub_base_vectors, 'R', 'ACG')
        assert result['C']['T'] == 5
        assert result['C']['C'] == 90
        assert result['G']['G'] == 95

    def test_all_zeros(self):
        sub_base_vectors = {
            'R_A': np.array([0, 0]),
            'R_C': np.array([0, 0]),
            'R_G': np.array([0, 0]),
            'R_T': np.array([0, 0]),
            'R_N': np.array([0, 0]),
        }
        result = prep_alternate_allele_counts(sub_base_vectors, 'R', 'AC')
        for a in 'ACGTN':
            for b in 'ACGTN':
                assert result[a][b] == 0

    def test_multiple_positions_same_base(self):
        sub_base_vectors = {
            'R_A': np.array([3, 2]),
            'R_C': np.array([85, 90]),
            'R_G': np.array([0, 0]),
            'R_T': np.array([12, 8]),
            'R_N': np.array([0, 0]),
        }
        result = prep_alternate_allele_counts(sub_base_vectors, 'R', 'CC')
        assert result['C']['A'] == 5   # 3 + 2
        assert result['C']['T'] == 20  # 12 + 8
        assert result['C']['C'] == 175 # 85 + 90

    def test_list_ref_sequence(self):
        sub_base_vectors = {
            'R_A': [10, 0],
            'R_C': [0, 5],
            'R_G': [0, 0],
            'R_T': [0, 0],
            'R_N': [0, 0],
        }
        result = prep_alternate_allele_counts(sub_base_vectors, 'R', ['A', 'C'])
        assert result['A']['A'] == 10
        assert result['C']['C'] == 5


# =============================================================================
# Tests: prep_global_modifications_reference (plot 4e/4f)
# =============================================================================


class TestPrepGlobalModificationsReference:

    def test_ref0_is_4e(self):
        ctx = _make_ctx(
            ref_names=['FANC', 'HDR'],
            refs={
                'FANC': _ref_dict(sequence_length=4, include_idxs=[0, 1, 2, 3]),
                'HDR': _ref_dict(),
            },
            N_TOTAL=100,
            ref1_all_insertion_count_vectors={'FANC': np.array([1, 0, 0, 0])},
            ref1_all_deletion_count_vectors={'FANC': np.array([0, 1, 0, 0])},
            ref1_all_substitution_count_vectors={'FANC': np.array([0, 0, 1, 0])},
        )
        ctx.ref_name = 'FANC'
        result = prep_global_modifications_reference(ctx)
        assert '4e' in result['fig_filename_root']
        assert 'all reads' in result['plot_title']

    def test_hdr_is_4f(self):
        ctx = _make_ctx(
            ref_names=['FANC', 'HDR'],
            refs={
                'FANC': _ref_dict(sequence_length=4, include_idxs=[0, 1, 2, 3]),
                'HDR': _ref_dict(sequence_length=4, include_idxs=[0, 1, 2, 3]),
            },
            N_TOTAL=100,
            ref1_all_insertion_count_vectors={'HDR': np.array([1, 0, 0, 0])},
            ref1_all_deletion_count_vectors={'HDR': np.array([0, 1, 0, 0])},
            ref1_all_substitution_count_vectors={'HDR': np.array([0, 0, 1, 0])},
        )
        ctx.ref_name = 'HDR'
        result = prep_global_modifications_reference(ctx)
        assert '4f' in result['fig_filename_root']
        assert 'HDR' in result['plot_title']


# =============================================================================
# Tests: prep_conversion_at_sel_nucs (plots 10e/10f/10g)
# =============================================================================


class TestPrepConversionAtSelNucs:

    def test_basic(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence='ACCA',
                sequence_length=4,
                sgRNA_plot_idxs=[[0, 1, 2, 3]],
            )},
            all_base_count_vectors={
                'r_A': np.array([80, 5, 5, 80]),
                'r_C': np.array([10, 85, 85, 10]),
                'r_G': np.array([5, 5, 5, 5]),
                'r_T': np.array([5, 5, 5, 5]),
                'r_N': np.array([0, 0, 0, 0]),
                'r_-': np.array([0, 0, 0, 0]),
            },
            counts_total={'r': 100},
            args=SimpleNamespace(conversion_nuc_from='C'),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0
        result = prep_conversion_at_sel_nucs(ctx)

        # C appears at positions 1, 2 in 'ACCA'
        assert len(result['from_nuc_indices']) == 2
        assert result['just_sel_nuc_pcts'].shape[1] == 2
        assert result['just_sel_nuc_freqs'].shape[1] == 2


# =============================================================================
# Tests: plot root generation
# =============================================================================


class TestPlotRootGeneration:
    """Verify that fig_filename_root/fig_filename_root is generated from ctx._jp."""

    def test_indel_size_uses_jp(self):
        def fake_jp(name):
            return '/output/' + name

        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(ref_plot_name='R.')},
            counts_total={'r': 100},
            args=SimpleNamespace(plot_histogram_outliers=False),
            _jp=fake_jp,
        )
        ctx.ref_name = 'r'
        result = prep_indel_size_distribution(ctx)
        assert result['fig_filename_root'] == '/output/3a.R.Indel_size_distribution'

    def test_no_jp_falls_back_to_filename(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(ref_plot_name='')},
            counts_total={'r': 100},
            args=SimpleNamespace(plot_histogram_outliers=False),
        )
        ctx.ref_name = 'r'
        result = prep_indel_size_distribution(ctx)
        assert result['fig_filename_root'] == '3a.Indel_size_distribution'


# =============================================================================
# Tests for amino_acids_to_numbers
# =============================================================================


def test_amino_acids_to_numbers_basic():
    """Test amino_acids_to_numbers with basic sequence."""
    result = amino_acids_to_numbers("MA")
    assert result == [11, 1]


def test_amino_acids_to_numbers_stop():
    """Test amino_acids_to_numbers with stop codon."""
    result = amino_acids_to_numbers("*")
    assert len(result) == 1
    assert result[0] == 0  # Stop codon is first in list


def test_amino_acids_to_numbers_gap():
    """Test amino_acids_to_numbers with gap."""
    result = amino_acids_to_numbers("-")
    assert result == [22]


def test_amino_acids_to_numbers_all_standard():
    """Test amino_acids_to_numbers with all standard amino acids."""
    all_aa = "ACDEFGHIKLMNPQRSTVWY"
    result = amino_acids_to_numbers(all_aa)
    assert result == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    # All should be unique
    assert len(set(result)) == 20


def test_amino_acids_to_numbers_empty():
    """Test amino_acids_to_numbers with empty sequence."""
    result = amino_acids_to_numbers("")
    assert result == []


def test_amino_acids_to_numbers_with_special():
    """Test amino_acids_to_numbers with special characters."""
    result = amino_acids_to_numbers("*-")
    assert len(result) == 2
    assert result[0] == 0  # * is first
    assert result[1] == 22  # - is last


# =============================================================================
# Tests for prep_alleles_table
# =============================================================================


def test_prep_alleles_table_basic():
    """Test prep_alleles_table with basic data."""
    df = pd.DataFrame({
        '%Reads': [50.0, 30.0, 20.0],
        '#Reads': [500, 300, 200],
        'Reference_Sequence': ['ATCG', 'ATCG', 'ATCG'],
    }, index=['ATCG', 'ATGG', 'A-CG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        prep_alleles_table(df, 'ATCG', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    assert X == [[1, 2, 3, 4], [1, 2, 4, 4], [1, 0, 3, 4]]
    assert annot == [['A', 'T', 'C', 'G'], ['A', 'T', 'G', 'G'], ['A', '-', 'C', 'G']]
    assert y_labels == ['50.00% (500 reads)', '30.00% (300 reads)', '20.00% (200 reads)']
    assert is_reference == [True, False, False]


def test_prep_alleles_table_empty():
    """Test prep_alleles_table with empty dataframe after filtering."""
    df = pd.DataFrame({
        '%Reads': [0.1],
        '#Reads': [1],
        'Reference_Sequence': ['ATCG'],
    }, index=['ATCG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        prep_alleles_table(df, 'ATCG', MAX_N_ROWS=10, MIN_FREQUENCY=1.0)

    assert len(X) == 0
    assert len(annot) == 0


def test_prep_alleles_table_max_rows():
    """Test prep_alleles_table respects MAX_N_ROWS."""
    df = pd.DataFrame({
        '%Reads': [30.0, 25.0, 20.0, 15.0, 10.0],
        '#Reads': [300, 250, 200, 150, 100],
        'Reference_Sequence': ['ATCG', 'ATCG', 'ATCG', 'ATCG', 'ATCG'],
    }, index=['ATCG', 'ATGG', 'TTCG', 'ATCA', 'GGGG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        prep_alleles_table(df, 'ATCG', MAX_N_ROWS=3, MIN_FREQUENCY=0)

    assert X == [[1, 2, 3, 4], [1, 2, 4, 4], [2, 2, 3, 4]]
    assert annot == [['A', 'T', 'C', 'G'], ['A', 'T', 'G', 'G'], ['T', 'T', 'C', 'G']]
    assert y_labels == ['30.00% (300 reads)', '25.00% (250 reads)', '20.00% (200 reads)']


def test_prep_alleles_table_with_insertions():
    """Test prep_alleles_table detects insertions."""
    # Reference with gap indicates insertion in read
    df = pd.DataFrame({
        '%Reads': [50.0],
        '#Reads': [500],
        'Reference_Sequence': ['AT--CG'],
    }, index=['ATGGCG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        prep_alleles_table(df, 'ATGGCG', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    # Should detect insertion
    assert 0 in insertion_dict
    assert len(insertion_dict[0]) > 0


def test_prep_alleles_table_with_substitution():
    """Test prep_alleles_table detects substitutions."""
    df = pd.DataFrame({
        '%Reads': [50.0],
        '#Reads': [500],
        'Reference_Sequence': ['ATCG'],  # Reference
    }, index=['GTCG'])  # G at position 0 instead of A

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        prep_alleles_table(df, 'ATCG', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    assert len(X) == 1
    assert is_reference[0] is False  # Different from reference
    # Should have bold annotation for substitution
    assert len(per_element_annot_kws[0]) > 0


def test_prep_alleles_table_all_reference():
    """Test prep_alleles_table with all reference sequences."""
    df = pd.DataFrame({
        '%Reads': [100.0],
        '#Reads': [1000],
        'Reference_Sequence': ['ATCG'],
    }, index=['ATCG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        prep_alleles_table(df, 'ATCG', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    assert len(is_reference) == 1
    assert is_reference[0] is True


# =============================================================================
# Tests for prep_alleles_table_compare
# =============================================================================


def test_prep_alleles_table_compare_basic():
    """Test prep_alleles_table_compare with basic data."""
    # Create merged allele dataframe
    df = pd.DataFrame({
        '%Reads_sample1': [50.0, 30.0],
        '%Reads_sample2': [40.0, 35.0],
        '#Reads_sample1': [500, 300],
        '#Reads_sample2': [400, 350],
        'Reference_Sequence': ['ATCG', 'ATCG'],
    }, index=['ATCG', 'ATGG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws = \
        prep_alleles_table_compare(
            df, 'sample1', 'sample2', MAX_N_ROWS=10, MIN_FREQUENCY=0
        )

    assert X == [[1, 2, 3, 4], [1, 2, 4, 4]]
    assert annot == [['A', 'T', 'C', 'G'], ['A', 'T', 'G', 'G']]
    assert y_labels == [
        '50.00% (500 reads) 40.00% (400 reads) ',
        '30.00% (300 reads) 35.00% (350 reads) ',
    ]


def test_prep_alleles_table_compare_with_insertion():
    """Test prep_alleles_table_compare detects insertions."""
    df = pd.DataFrame({
        '%Reads_s1': [50.0],
        '%Reads_s2': [50.0],
        '#Reads_s1': [500],
        '#Reads_s2': [500],
        'Reference_Sequence': ['AT--CG'],  # Insertion markers
    }, index=['ATGGCG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws = \
        prep_alleles_table_compare(
            df, 's1', 's2', MAX_N_ROWS=10, MIN_FREQUENCY=0
        )

    # Should detect insertion
    assert 0 in insertion_dict
    assert len(insertion_dict[0]) > 0


# =============================================================================
# Tests for prep_amino_acid_table_for_plot
# =============================================================================


def test_prep_amino_acid_table_for_plot_basic():
    """Test prep_amino_acid_table_for_plot with basic data."""
    df = pd.DataFrame({
        '%Reads': [60.0, 40.0],
        '#Reads': [600, 400],
        'Reference_Sequence': ['MAS', 'MAS'],
        'silent_edit_inds': [[], []],
    }, index=['MAS', 'MAT'])

    X, annot, y_labels, insertion_dict, silent_edit_dict, per_element_annot_kws, is_reference, ref_seq = \
        prep_amino_acid_table_for_plot(df, 'MAS', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    assert len(X) == 2
    assert is_reference[0] is True  # First row matches reference
    assert ref_seq == 'MAS'


def test_prep_amino_acid_table_for_plot_with_silent_edits():
    """Test prep_amino_acid_table_for_plot with silent edits."""
    df = pd.DataFrame({
        '%Reads': [50.0],
        '#Reads': [500],
        'Reference_Sequence': ['MAS'],
        'silent_edit_inds': [[1]],  # Silent edit at position 1
    }, index=['MAS'])

    X, annot, y_labels, insertion_dict, silent_edit_dict, per_element_annot_kws, is_reference, ref_seq = \
        prep_amino_acid_table_for_plot(df, 'MAS', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    # Should have silent edit at row 0, position 1
    assert 0 in silent_edit_dict
    assert 1 in silent_edit_dict[0]


# =============================================================================
# Helper: Build a minimal df_alleles for integration-style prep tests
# =============================================================================


def _make_df_alleles(ref_name, aligned_seqs, ref_seq, reads, read_status=None):
    """Build a df_alleles-like DataFrame suitable for CRISPRessoShared cut functions.

    Parameters
    ----------
    ref_name : str
        Reference name to populate ``Reference_Name``.
    aligned_seqs : list[str]
        Aligned read sequences (may contain ``-`` for deletions).
    ref_seq : str
        The aligned reference sequence (same length as each aligned_seq).
    reads : list[int]
        ``#Reads`` for each allele.
    read_status : list[str] or None
        ``Read_Status`` values; defaults to ``'MODIFIED'`` for all.

    Returns a DataFrame with the columns expected by
    ``CRISPRessoShared.get_dataframe_around_cut_asymmetrical`` et al.
    """
    total_reads = sum(reads)
    n = len(aligned_seqs)
    if read_status is None:
        read_status = ['MODIFIED'] * n

    # Build ref_positions: for each position in the alignment, what is the
    # reference coordinate (0-indexed)?  Gaps in the reference get -1.
    ref_positions = []
    pos = 0
    for c in ref_seq:
        if c == '-':
            ref_positions.append(-1)
        else:
            ref_positions.append(pos)
            pos += 1

    # Compute simple n_deleted / n_inserted / n_mutated per allele
    rows = []
    for i, (aln, status) in enumerate(zip(aligned_seqs, read_status)):
        n_del = sum(1 for a, r in zip(aln, ref_seq) if a == '-' and r != '-')
        n_ins = sum(1 for a, r in zip(aln, ref_seq) if a != '-' and r == '-')
        n_mut = sum(1 for a, r in zip(aln, ref_seq) if a != '-' and r != '-' and a != r)
        rows.append({
            'Aligned_Sequence': aln,
            'Reference_Sequence': ref_seq,
            'ref_positions': ref_positions,
            'Read_Status': status,
            'n_deleted': n_del,
            'n_inserted': n_ins,
            'n_mutated': n_mut,
            '#Reads': reads[i],
            '%Reads': 100.0 * reads[i] / total_reads,
            'Reference_Name': ref_name,
        })
    return pd.DataFrame(rows)


# =============================================================================
# Tests: prep_alleles_around_cut (plot 9)
# =============================================================================


class TestPrepAllelesAroundCut:

    def test_basic(self):
        """Smoke test: single allele, unmodified, returns expected keys."""
        #                       0123456789
        ref_sequence =         'AACCGGTTAA'
        aligned_ref =          'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=[ref_sequence],
            ref_seq=aligned_ref,
            reads=[100],
            read_status=['UNMODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
                sgRNA_names=['sgRNA1'],
                sgRNA_mismatches=[],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)

        assert 'df_alleles_around_cut' in result
        assert 'ref_seq_around_cut' in result
        assert 'new_sgRNA_intervals' in result
        assert 'new_cut_point' in result
        assert 'window_truncated' in result
        assert 'plot_input' in result
        # Window of size 3 around cut_point=4 → positions 2..7 → 6 chars
        assert len(result['ref_seq_around_cut']) <= 7

    def test_with_substitution(self):
        """Allele with a substitution produces plot_input with expected content."""
        ref_sequence = 'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=['AACCGGTTAA', 'AATCGGTTAA'],
            ref_seq='AACCGGTTAA',
            reads=[80, 20],
            read_status=['UNMODIFIED', 'MODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
                sgRNA_names=['sgRNA1'],
                sgRNA_mismatches=[],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)

        assert result['plot_input'] is not None
        pi = result['plot_input']
        assert 'reference_seq' in pi
        assert 'prepped_df_alleles' in pi
        assert 'fig_filename_root' in pi
        assert '9.' in pi['fig_filename_root']
        # Two alleles in the window
        assert len(pi['y_labels']) >= 1

    def test_plot_input_none_below_threshold(self):
        """When all alleles are below min_frequency, plot_input is None."""
        ref_sequence = 'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=['AACCGGTTAA'],
            ref_seq='AACCGGTTAA',
            reads=[100],
            read_status=['UNMODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                min_frequency_alleles_around_cut_to_plot=99999,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)
        assert result['plot_input'] is None

    def test_window_truncation_near_edge(self):
        """Window near amplicon edge sets window_truncated=True."""
        ref_sequence = 'ACGT'
        df = _make_df_alleles(
            'r', aligned_seqs=['ACGT'], ref_seq='ACGT',
            reads=[100], read_status=['UNMODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=4,
                sgRNA_cut_points=[0],
                sgRNA_plot_cut_points=[0],
                sgRNA_intervals=[(0, 3)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=20,  # much larger than sequence
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)
        assert result['window_truncated'] is True


# =============================================================================
# Tests: prep_base_edit_quilt (plot 10h)
# =============================================================================


class TestPrepBaseEditQuilt:

    def test_basic(self):
        """Smoke test: base edit quilt with one C position."""
        #                 0123456789
        ref_sequence =   'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=['AACCGGTTAA', 'AATCGGTTAA'],
            ref_seq='AACCGGTTAA',
            reads=[80, 20],
            read_status=['UNMODIFIED', 'MODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                conversion_nuc_from='C',
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_base_edit_quilt(ctx)

        assert 'df_alleles_around_cut' in result
        assert 'ref_seq_around_cut' in result
        assert 'plot_input' in result

    def test_plot_input_none_below_threshold(self):
        """When all alleles are below min_frequency, plot_input is None."""
        ref_sequence = 'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=['AACCGGTTAA'],
            ref_seq='AACCGGTTAA',
            reads=[100],
            read_status=['UNMODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                conversion_nuc_from='C',
                min_frequency_alleles_around_cut_to_plot=99999,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_base_edit_quilt(ctx)
        assert result['plot_input'] is None

    def test_x_labels_are_conversion_nuc_positions(self):
        """x_labels should be 1-indexed positions of conversion_nuc_from."""
        ref_sequence = 'ACCA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=['ACCA', 'ATCA'],
            ref_seq='ACCA',
            reads=[70, 30],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=4,
                sgRNA_cut_points=[1],
                sgRNA_plot_cut_points=[1],
                sgRNA_intervals=[(0, 3)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                conversion_nuc_from='C',
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_base_edit_quilt(ctx)
        # 'ACCA' — C at positions 1 and 2 (0-indexed) → 1-indexed: [2, 3]
        pi = result['plot_input']
        if pi is not None:
            assert pi['x_labels'] == [2, 3]


# =============================================================================
# Tests: prep_amino_acid_table (plot 9a)
# =============================================================================


class TestPrepAminoAcidTable:

    def test_basic(self):
        """Smoke test: amino acid table with a simple coding seq."""
        # 9-base coding seq → 3 amino acids
        coding_seq = 'ATGGCTTAA'  # M, A, *
        # Build a ref that contains this coding seq starting at position 0
        ref_sequence = coding_seq
        ref_len = len(ref_sequence)

        df = _make_df_alleles(
            'r',
            aligned_seqs=[ref_sequence, 'ATGACTTAA'],
            ref_seq=ref_sequence,
            reads=[80, 20],
            read_status=['UNMODIFIED', 'MODIFIED'],
        )

        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=ref_len,
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
                sgRNA_names=['sgRNA1'],
                sgRNA_mismatches=[],
                contains_coding_seq=True,
                exon_positions=list(range(ref_len)),
                exon_intervals=[(0, ref_len - 1)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                annotate_wildtype_allele='',
            ),
            run_data={'running_info': {'coding_seqs': [coding_seq]}},
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0
        ctx.coding_seq_ind = 0

        result = prep_amino_acid_table(ctx)

        assert 'coding_seq_amino_acids' in result
        assert 'amino_acid_cut_point' in result
        assert 'df_to_plot' in result
        assert 'plot_input' in result
        # Coding seq ATG GCT TAA → M A *
        assert result['coding_seq_amino_acids'] == 'MA*'

    def test_plot_input_none_below_threshold(self):
        """When no rows pass threshold, plot_input is None."""
        coding_seq = 'ATGGCTTAA'
        ref_sequence = coding_seq

        df = _make_df_alleles(
            'r',
            aligned_seqs=[ref_sequence],
            ref_seq=ref_sequence,
            reads=[100],
            read_status=['UNMODIFIED'],
        )

        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
                contains_coding_seq=True,
                exon_positions=list(range(len(ref_sequence))),
                exon_intervals=[(0, len(ref_sequence) - 1)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                min_frequency_alleles_around_cut_to_plot=99999,
                max_rows_alleles_around_cut_to_plot=10,
                annotate_wildtype_allele='',
            ),
            run_data={'running_info': {'coding_seqs': [coding_seq]}},
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0
        ctx.coding_seq_ind = 0

        result = prep_amino_acid_table(ctx)
        assert result['plot_input'] is None
