"""Tests for CRISPResso2.plots.data_prep — extracted data preparation functions.

Each prep function takes raw analysis data and returns the kwargs dict
that the corresponding CRISPRessoPlot function expects. Only functions
with non-trivial computation are extracted and tested here.
"""

import numpy as np
from inline_snapshot import snapshot

from CRISPResso2.plots.data_prep import (
    _prep_windowed_alleles,
    prep_alleles_around_cut,
    prep_amino_acid_table,
    prep_amplicon_modifications,
    prep_base_edit_quilt,
    prep_dsODN_piechart,
    prep_frequency_deletions_insertions,
    prep_global_frameshift_data,
    prep_hdr_nucleotide_quilt,
    prep_indel_size_distribution,
    prep_modification_frequency,
    prep_nucleotide_quilt,
    prep_pe_nucleotide_quilt,
    prep_pe_nucleotide_quilt_around_sgRNA,
)

import pandas as pd


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


class TestPrepDsODNPiechart:
    """prep_dsODN_piechart (plot_1d) — labels/sizes from df_alleles."""

    def _make_df(self, n_contain, n_not_contain):
        import pandas as pd
        return pd.DataFrame({
            'contains dsODN': [True] * n_contain + [False] * n_not_contain,
            '#Reads': [1] * (n_contain + n_not_contain),
        })

    def test_basic(self):
        df = self._make_df(n_contain=30, n_not_contain=70)
        result = prep_dsODN_piechart(
            df_alleles=df,
            N_TOTAL=100,
            plot_root='/tmp/1d',
            save_also_png=True,
        )
        assert result == snapshot({'sizes': [np.float64(30.0), np.float64(70.0)], 'labels': ["""\
Contains dsODN
(30 reads)\
""", """\
No dsODN
(70 reads)\
"""], 'plot_root': '/tmp/1d', 'save_also_png': True})

    def test_none_contain(self):
        df = self._make_df(n_contain=0, n_not_contain=50)
        result = prep_dsODN_piechart(
            df_alleles=df,
            N_TOTAL=50,
            plot_root='/tmp/1d',
            save_also_png=False,
        )
        assert result == snapshot({'sizes': [np.float64(0.0), np.float64(100.0)], 'labels': ["""\
Contains dsODN
(0 reads)\
""", """\
No dsODN
(50 reads)\
"""], 'plot_root': '/tmp/1d', 'save_also_png': False})


class TestPrepNucleotideQuilt:
    """prep_nucleotide_quilt (plot_2a) — builds nuc/mod DataFrames."""

    def _make_inputs(self, ref_name='R', ref_seq='ACGT', counts_total=100):
        n = len(ref_seq)
        all_base = {
            ref_name + '_A': np.array([counts_total, 0, 0, 0][:n]),
            ref_name + '_C': np.array([0, counts_total, 0, 0][:n]),
            ref_name + '_G': np.array([0, 0, counts_total, 0][:n]),
            ref_name + '_T': np.array([0, 0, 0, counts_total][:n]),
            ref_name + '_N': np.zeros(n, dtype=int),
            ref_name + '_-': np.zeros(n, dtype=int),
        }
        return dict(
            all_base_count_vectors=all_base,
            all_insertion_count_vectors=np.zeros(n, dtype=int),
            all_insertion_left_count_vectors=np.zeros(n, dtype=int),
            all_deletion_count_vectors=np.zeros(n, dtype=int),
            all_substitution_count_vectors=np.zeros(n, dtype=int),
            all_indelsub_count_vectors=np.zeros(n, dtype=int),
            counts_total=counts_total,
            ref_name=ref_name,
            ref_seq=ref_seq,
            include_idxs_list=[1, 2],
            sgRNA_intervals=[(1, 2)],
            sgRNA_names=['sg1'],
            sgRNA_mismatches=[0],
            sgRNA_sequences=['ACGT'],
            plot_root='/tmp/2a',
            save_also_png=False,
            custom_colors={},
        )

    def test_basic_structure(self):
        from CRISPResso2.plots.data_prep import prep_nucleotide_quilt
        result = prep_nucleotide_quilt(**self._make_inputs())
        # Verify the DataFrame shapes — don't snapshot the full DataFrames
        assert list(result.keys()) == snapshot(['nuc_pct_df', 'mod_pct_df', 'fig_filename_root', 'save_also_png', 'sgRNA_intervals', 'sgRNA_names', 'sgRNA_mismatches', 'sgRNA_sequences', 'quantification_window_idxs', 'custom_colors'])
        assert result['nuc_pct_df'].shape == snapshot((6, 6))
        assert result['mod_pct_df'].shape == snapshot((6, 6))
        assert result['fig_filename_root'] == '/tmp/2a'
        assert result['quantification_window_idxs'] == [1, 2]

    def test_nuc_percentages_sum_to_one(self):
        from CRISPResso2.plots.data_prep import prep_nucleotide_quilt
        result = prep_nucleotide_quilt(**self._make_inputs())
        # Each column (position) of nuc_pct_df (excluding Batch/Nucleotide columns)
        # should sum to 1.0 (all reads assigned to one nucleotide in this fixture)
        nuc_df = result['nuc_pct_df']
        seq_cols = nuc_df.columns[2:]  # skip Batch and Nucleotide
        col_sums = nuc_df[seq_cols].sum(axis=0)
        np.testing.assert_allclose(col_sums.values, 1.0)


class TestPrepNucleotideQuiltAroundSgRNA:
    """prep_nucleotide_quilt_around_sgRNA (plot_2b) — window slicing."""

    def _make_dfs(self, ref_name='R', ref_seq='ACGTACGTACGT'):
        """Build minimal nuc/mod DataFrames matching the 2a output shape."""
        import pandas as pd
        n = len(ref_seq)
        # nuc_df: shape (6 nucs × (n+2 cols: Batch + Nucleotide + n positions))
        nuc_rows = []
        for nuc in ['A', 'C', 'G', 'T', 'N', '-']:
            row = [ref_name, nuc] + [0.0] * n
            nuc_rows.append(row)
        nuc_cols = ['Batch', 'Nucleotide'] + list(ref_seq)
        nuc_df = pd.DataFrame(nuc_rows, columns=nuc_cols)

        # mod_df: shape (6 mods × (n+2 cols))
        mod_rows = []
        for mod in ['Insertions', 'Insertions_Left', 'Deletions', 'Substitutions', 'All_modifications', 'Total']:
            row = [ref_name, mod] + [0.0] * n
            mod_rows.append(row)
        mod_cols = ['Batch', 'Modification'] + list(ref_seq)
        mod_df = pd.DataFrame(mod_rows, columns=mod_cols)
        return nuc_df, mod_df

    def test_window_slicing(self):
        from CRISPResso2.plots.data_prep import prep_nucleotide_quilt_around_sgRNA
        nuc_df, mod_df = self._make_dfs()
        result = prep_nucleotide_quilt_around_sgRNA(
            nuc_df_for_plot=nuc_df,
            mod_df_for_plot=mod_df,
            cut_point=6,
            plot_half_window=3,
            ref_len=12,
            sgRNA_intervals=[(4, 8)],
            include_idxs_list=[5, 6, 7],
            sgRNA_names=['sg1'],
            sgRNA_mismatches=[0],
            sgRNA_sequences=['ACGT'],
            plot_root='/tmp/2b',
            save_also_png=False,
            custom_colors={},
        )
        assert list(result.keys()) == snapshot(['nuc_pct_df', 'mod_pct_df', 'fig_filename_root', 'save_also_png', 'sgRNA_intervals', 'sgRNA_names', 'sgRNA_mismatches', 'sgRNA_sequences', 'quantification_window_idxs', 'custom_colors'])
        # Window around cut_point=6, half=3 → cols 4..9, sel_cols=[0,1,6,7,8,9,10,11]
        assert result['nuc_pct_df'].shape == snapshot((6, 8))
        assert result['sgRNA_intervals'] == snapshot([(0, 4)])
        assert result['quantification_window_idxs'] == snapshot([1, 2, 3])

    def test_clamped_near_start(self):
        """When cut is near the start, new_sel_cols_start clamps to 2."""
        from CRISPResso2.plots.data_prep import prep_nucleotide_quilt_around_sgRNA
        nuc_df, mod_df = self._make_dfs()
        result = prep_nucleotide_quilt_around_sgRNA(
            nuc_df_for_plot=nuc_df,
            mod_df_for_plot=mod_df,
            cut_point=1,
            plot_half_window=5,
            ref_len=12,
            sgRNA_intervals=[(0, 2)],
            include_idxs_list=[1],
            sgRNA_names=['sg1'],
            sgRNA_mismatches=[0],
            sgRNA_sequences=['ACGT'],
            plot_root='/tmp/2b',
            save_also_png=False,
            custom_colors={},
        )
        # new_sel_cols_start = max(2, 1-5+1) = 2
        assert result['sgRNA_intervals'] == snapshot([(-2, 0)])
        assert result['quantification_window_idxs'] == snapshot([-1])


class TestPrepGlobalModificationsReference:
    """prep_global_modifications_reference (plot_4e/4f) — passes through plot_root/plot_title."""

    def _inputs(self, ref_name, ref_names, plot_root, plot_title):
        return dict(
            ref_name=ref_name,
            ref_names=ref_names,
            ref1_all_insertion_count_vectors=np.array([1, 2, 3]),
            ref1_all_deletion_count_vectors=np.array([0, 1, 0]),
            ref1_all_substitution_count_vectors=np.array([0, 0, 1]),
            ref1={'sequence': 'ACG', 'include_idxs': [0, 1, 2]},
            include_idxs_list=[0, 1, 2],
            N_TOTAL=100,
            ref_len=3,
            custom_colors={},
            save_also_png=False,
            plot_root=plot_root,
            plot_title=plot_title,
        )

    def test_primary_ref_uses_4e(self):
        from CRISPResso2.plots.data_prep import prep_global_modifications_reference
        result = _to_serializable(prep_global_modifications_reference(
            **self._inputs(
                ref_name='FANC',
                ref_names=['FANC', 'HDR'],
                plot_root='/out/4e.FANC.Global_mutations_in_all_reads',
                plot_title='Mutation position distribution in all reads with reference to FANC',
            )
        ))
        assert result == snapshot({'custom_colors': {}, 'include_idxs_list': [0, 1, 2], 'n_total': 100, 'plot_root': '/out/4e.FANC.Global_mutations_in_all_reads', 'plot_title': 'Mutation position distribution in all reads with reference to FANC', 'ref1': {'include_idxs': [0, 1, 2], 'sequence': 'ACG'}, 'ref1_all_deletion_count_vectors': [0, 1, 0], 'ref1_all_insertion_count_vectors': [1, 2, 3], 'ref1_all_substitution_count_vectors': [0, 0, 1], 'ref_len': 3, 'ref_name': 'FANC', 'save_also_png': False})

    def test_hdr_ref_uses_4f(self):
        from CRISPResso2.plots.data_prep import prep_global_modifications_reference
        result = _to_serializable(prep_global_modifications_reference(
            **self._inputs(
                ref_name='HDR',
                ref_names=['FANC', 'HDR'],
                plot_root='/out/4f.FANC.Global_mutations_in_HDR_reads_with_reference_to_FANC',
                plot_title='Mutation position distribution in HDR reads with reference to FANC',
            )
        ))
        assert result == snapshot({'custom_colors': {}, 'include_idxs_list': [0, 1, 2], 'n_total': 100, 'plot_root': '/out/4f.FANC.Global_mutations_in_HDR_reads_with_reference_to_FANC', 'plot_title': 'Mutation position distribution in HDR reads with reference to FANC', 'ref1': {'include_idxs': [0, 1, 2], 'sequence': 'ACG'}, 'ref1_all_deletion_count_vectors': [0, 1, 0], 'ref1_all_insertion_count_vectors': [1, 2, 3], 'ref1_all_substitution_count_vectors': [0, 0, 1], 'ref_len': 3, 'ref_name': 'FANC', 'save_also_png': False})


class TestPrepLogNucFreqs:
    """prep_log_nuc_freqs (plot_10d) — computes plot_quant_window_idxs."""

    def _make_df(self, n_cols=10):
        import pandas as pd
        return pd.DataFrame(
            np.zeros((6, n_cols)),
            index=['A', 'C', 'G', 'T', 'N', '-'],
        )

    def test_basic_window_computation(self):
        from CRISPResso2.plots.data_prep import prep_log_nuc_freqs
        df = self._make_df(10)
        result = prep_log_nuc_freqs(
            df_nuc_freq=df,
            tot_aln_reads=100,
            include_idxs_list=[3, 4, 5],
            plot_idxs=[2, 3, 4, 5, 6, 7, 8],  # window: positions 2-8
            ref_len=10,
            sgRNA_legend='sgRNA ATCG',
            ref_name='FANC',
            ref_names=['FANC'],
            fig_filename_root='/tmp/10d',
            save_also_png=False,
        )
        # Only positions 3,4,5 are in quant window; their local indices
        # in plot_idxs are 1,2,3 → after -2 offset: -1,0,1
        assert result['quantification_window_idxs'] == snapshot([-1, 0, 1])
        assert result['plot_title'] == snapshot('Log2 Nucleotide Frequencies Around the sgRNA ATCG')

    def test_multi_ref_title(self):
        from CRISPResso2.plots.data_prep import prep_log_nuc_freqs
        df = self._make_df(5)
        result = prep_log_nuc_freqs(
            df_nuc_freq=df,
            tot_aln_reads=50,
            include_idxs_list=[0, 1],
            plot_idxs=[0, 1, 2],
            ref_len=5,
            sgRNA_legend='myGuide',
            ref_name='REF1',
            ref_names=['REF1', 'REF2'],
            fig_filename_root='/tmp/10d',
            save_also_png=True,
        )
        assert result['plot_title'] == snapshot('Log2 Nucleotide Frequencies Around the myGuide: REF1')


class TestPrepConversionAtSelNucs:
    """prep_conversion_at_sel_nucs (plot_10e/10f/10g) — from_nuc_indices + sliced dfs."""

    def _make_pct_df(self, seq):
        """Build a minimal nuc percentage DataFrame with column names = seq chars."""
        import pandas as pd
        n = len(seq)
        data = np.zeros((6, n))
        df = pd.DataFrame(data, index=['A', 'C', 'G', 'T', 'N', '-'], columns=list(seq))
        return df

    def test_finds_correct_indices(self):
        from CRISPResso2.plots.data_prep import prep_conversion_at_sel_nucs
        pcts = self._make_pct_df('ACGCTA')
        freqs = self._make_pct_df('ACGCTA')
        result = prep_conversion_at_sel_nucs(
            plot_nuc_pcts=pcts,
            plot_nuc_freqs=freqs,
            conversion_nuc_from='C',
        )
        assert result['from_nuc_indices'] == snapshot([1, 3])
        assert list(result['just_sel_nuc_pcts'].columns) == snapshot(['C2', 'C4'])
        assert list(result['just_sel_nuc_freqs'].columns) == snapshot(['C2', 'C4'])

    def test_no_matching_nucs(self):
        from CRISPResso2.plots.data_prep import prep_conversion_at_sel_nucs
        pcts = self._make_pct_df('AAAA')
        freqs = self._make_pct_df('AAAA')
        result = prep_conversion_at_sel_nucs(
            plot_nuc_pcts=pcts,
            plot_nuc_freqs=freqs,
            conversion_nuc_from='C',
        )
        assert result['from_nuc_indices'] == snapshot([])
        assert result['just_sel_nuc_pcts'].shape == snapshot((6, 0))


class TestPrepHdrNucleotideQuilt:
    """prep_hdr_nucleotide_quilt (plot_4g) — cross-ref nuc/mod DataFrames."""

    def _make_inputs(self, ref_seq='ACGT'):
        n = len(ref_seq)
        ref_names = ['FANC', 'HDR']
        counts_total = {'FANC': 100, 'HDR': 50}
        # Base counts: all reads match the reference at each position
        ref1_all_base_count_vectors = {}
        for rn in ref_names:
            for nuc in ['A', 'C', 'G', 'T', 'N', '-']:
                ref1_all_base_count_vectors[rn + '_' + nuc] = np.zeros(n)
            # Put all counts on the matching base at each position
            for i, base in enumerate(ref_seq):
                ref1_all_base_count_vectors[rn + '_' + base][i] = counts_total[rn]

        ref1_all_insertion_count_vectors = {rn: np.zeros(n) for rn in ref_names}
        ref1_all_insertion_left_count_vectors = {rn: np.zeros(n) for rn in ref_names}
        ref1_all_deletion_count_vectors = {rn: np.zeros(n) for rn in ref_names}
        ref1_all_substitution_count_vectors = {rn: np.zeros(n) for rn in ref_names}
        ref1_all_indelsub_count_vectors = {rn: np.zeros(n) for rn in ref_names}

        refs = {
            'FANC': {
                'sequence': ref_seq,
                'sequence_length': n,
                'sgRNA_intervals': [(1, 2)],
                'sgRNA_names': ['sg1'],
                'sgRNA_mismatches': [0],
                'sgRNA_sequences': ['ACGT'],
            },
            'HDR': {
                'sequence': ref_seq,
                'sequence_length': n,
                'sgRNA_intervals': [(1, 2)],
                'sgRNA_names': ['sg1'],
                'sgRNA_mismatches': [0],
                'sgRNA_sequences': ['ACGT'],
            },
        }

        return dict(
            ref_names_for_hdr=ref_names,
            counts_total=counts_total,
            refs=refs,
            ref1_all_base_count_vectors=ref1_all_base_count_vectors,
            ref1_all_insertion_count_vectors=ref1_all_insertion_count_vectors,
            ref1_all_insertion_left_count_vectors=ref1_all_insertion_left_count_vectors,
            ref1_all_deletion_count_vectors=ref1_all_deletion_count_vectors,
            ref1_all_substitution_count_vectors=ref1_all_substitution_count_vectors,
            ref1_all_indelsub_count_vectors=ref1_all_indelsub_count_vectors,
            custom_colors={},
            save_also_png=False,
        )

    def test_basic_structure(self):
        result = prep_hdr_nucleotide_quilt(**self._make_inputs())
        assert list(result.keys()) == snapshot(
            [
                "nuc_pct_df",
                "mod_pct_df",
                "fig_filename_root",
                "save_also_png",
                "sgRNA_intervals",
                "sgRNA_names",
                "sgRNA_mismatches",
                "sgRNA_sequences",
                "quantification_window_idxs",
                "custom_colors",
            ]
        )
        # 2 refs × 6 nucleotides = 12 rows; 4 positions + Batch + Nucleotide = 6 cols
        assert result['nuc_pct_df'].shape == snapshot((12, 6))
        # 2 refs × 6 mod types = 12 rows
        assert result['mod_pct_df'].shape == snapshot((12, 6))
        # quantification window should be empty (different amplicons may differ)
        assert result['quantification_window_idxs'] == snapshot([])

    def test_nuc_percentages_sum_to_one_per_batch(self):
        result = prep_hdr_nucleotide_quilt(**self._make_inputs())
        nuc_df = result['nuc_pct_df']
        for batch in nuc_df['Batch'].unique():
            batch_df = nuc_df[nuc_df['Batch'] == batch]
            seq_cols = batch_df.columns[2:]
            col_sums = batch_df[seq_cols].sum(axis=0)
            np.testing.assert_allclose(col_sums.values, 1.0)

    def test_single_ref(self):
        inputs = self._make_inputs()
        inputs['ref_names_for_hdr'] = ['FANC']
        result = prep_hdr_nucleotide_quilt(**inputs)
        # 1 ref × 6 nucleotides = 6 rows
        assert result['nuc_pct_df'].shape == snapshot((6, 6))
        assert result['mod_pct_df'].shape == snapshot((6, 6))

    def test_sgRNA_from_first_ref(self):
        result = prep_hdr_nucleotide_quilt(**self._make_inputs())
        assert result['sgRNA_intervals'] == snapshot([(1, 2)])
        assert result['sgRNA_names'] == snapshot(['sg1'])


class TestPrepGlobalFrameshiftData:
    """prep_global_frameshift_data (plot_5a/6a/8a) — cross-ref aggregation."""

    def _make_inputs(self):
        from collections import Counter
        ref_names = ['FANC', 'HDR']
        refs = {
            'FANC': {'contains_coding_seq': True},
            'HDR': {'contains_coding_seq': True},
        }
        counts_modified_frameshift = {'FANC': 10, 'HDR': 5}
        counts_modified_non_frameshift = {'FANC': 20, 'HDR': 8}
        counts_non_modified_non_frameshift = {'FANC': 30, 'HDR': 12}
        counts_splicing_sites_modified = {'FANC': 2, 'HDR': 1}
        counts_total = {'FANC': 100, 'HDR': 50}
        counts_modified = {'FANC': 30, 'HDR': 13}
        counts_unmodified = {'FANC': 70, 'HDR': 37}
        hists_frameshift = {
            'FANC': Counter({0: 5, 1: 3, -1: 2}),
            'HDR': Counter({0: 2, 2: 1}),
        }
        hists_inframe = {
            'FANC': Counter({0: 10, 3: 5}),
            'HDR': Counter({0: 4, -3: 2}),
        }
        return dict(
            ref_names=ref_names,
            refs=refs,
            counts_modified_frameshift=counts_modified_frameshift,
            counts_modified_non_frameshift=counts_modified_non_frameshift,
            counts_non_modified_non_frameshift=counts_non_modified_non_frameshift,
            counts_splicing_sites_modified=counts_splicing_sites_modified,
            counts_total=counts_total,
            counts_modified=counts_modified,
            counts_unmodified=counts_unmodified,
            hists_frameshift=hists_frameshift,
            hists_inframe=hists_inframe,
        )

    def test_basic_aggregation(self):
        result = prep_global_frameshift_data(**self._make_inputs())
        # FANC: not HDR, so standard aggregation
        # HDR: adds unmodified reads to non_modified_non_frameshift
        assert result['global_modified_frameshift'] == snapshot(15)  # 10 + 5
        assert result['global_modified_non_frameshift'] == snapshot(28)  # 20 + 8
        # FANC: 30, HDR: 12 + 37 (unmodified added for HDR)
        assert result['global_non_modified_non_frameshift'] == snapshot(79)
        assert result['global_splicing_sites_modified'] == snapshot(3)  # 2 + 1
        assert result['global_count_total'] == snapshot(150)
        assert result['global_count_modified'] == snapshot(43)
        assert result['global_count_unmodified'] == snapshot(107)

    def test_histogram_aggregation(self):
        from collections import Counter
        result = prep_global_frameshift_data(**self._make_inputs())
        # Frameshift histograms: FANC {0:5, 1:3, -1:2} + HDR {0:2, 2:1}
        assert result['global_hists_frameshift'] == snapshot(
            Counter({0: 7, 1: 3, -1: 2, 2: 1})
        )
        # Inframe histograms: FANC {0:10, 3:5} + HDR {0:4, -3:2}
        assert result['global_hists_inframe'] == snapshot(
            Counter({0: 14, 3: 5, -3: 2})
        )

    def test_non_coding_refs_skipped(self):
        inputs = self._make_inputs()
        inputs['refs']['FANC']['contains_coding_seq'] = False
        result = prep_global_frameshift_data(**inputs)
        # Only HDR counts
        assert result['global_modified_frameshift'] == snapshot(5)
        assert result['global_count_total'] == snapshot(50)

    def test_no_coding_refs(self):
        inputs = self._make_inputs()
        inputs['refs']['FANC']['contains_coding_seq'] = False
        inputs['refs']['HDR']['contains_coding_seq'] = False
        result = prep_global_frameshift_data(**inputs)
        assert result['global_modified_frameshift'] == snapshot(0)
        assert result['global_count_total'] == snapshot(0)


class TestPrepPeNucleotideQuilt:
    """prep_pe_nucleotide_quilt (plot_11a) — cross-ref nuc/mod DataFrames with include_idxs."""

    def _make_inputs(self, ref_seq='ACGT'):
        n = len(ref_seq)
        ref_names = ['Prime-edited', 'Scaffold-incorporated']
        counts_total = {'Prime-edited': 100, 'Scaffold-incorporated': 20}
        ref1_all_base_count_vectors = {}
        for rn in ref_names:
            for nuc in ['A', 'C', 'G', 'T', 'N', '-']:
                ref1_all_base_count_vectors[rn + '_' + nuc] = np.zeros(n)
            for i, base in enumerate(ref_seq):
                ref1_all_base_count_vectors[rn + '_' + base][i] = counts_total[rn]

        ref1_all_insertion_count_vectors = {rn: np.zeros(n) for rn in ref_names}
        ref1_all_insertion_left_count_vectors = {rn: np.zeros(n) for rn in ref_names}
        ref1_all_deletion_count_vectors = {rn: np.zeros(n) for rn in ref_names}
        ref1_all_substitution_count_vectors = {rn: np.zeros(n) for rn in ref_names}
        ref1_all_indelsub_count_vectors = {rn: np.zeros(n) for rn in ref_names}

        refs = {
            'Prime-edited': {
                'sequence': ref_seq,
                'sequence_length': n,
                'sgRNA_intervals': [(1, 2)],
                'sgRNA_names': ['sg1'],
                'sgRNA_mismatches': [0],
                'sgRNA_sequences': ['ACGT'],
                'include_idxs': [1, 2],
            },
            'Scaffold-incorporated': {
                'sequence': ref_seq,
                'sequence_length': n,
                'sgRNA_intervals': [(1, 2)],
                'sgRNA_names': ['sg1'],
                'sgRNA_mismatches': [0],
                'sgRNA_sequences': ['ACGT'],
                'include_idxs': [1, 2],
            },
        }

        return dict(
            ref_names_for_pe=ref_names,
            ref_names=ref_names,
            counts_total=counts_total,
            refs=refs,
            ref1_all_base_count_vectors=ref1_all_base_count_vectors,
            ref1_all_insertion_count_vectors=ref1_all_insertion_count_vectors,
            ref1_all_insertion_left_count_vectors=ref1_all_insertion_left_count_vectors,
            ref1_all_deletion_count_vectors=ref1_all_deletion_count_vectors,
            ref1_all_substitution_count_vectors=ref1_all_substitution_count_vectors,
            ref1_all_indelsub_count_vectors=ref1_all_indelsub_count_vectors,
            custom_colors={},
            save_also_png=False,
        )

    def test_basic_structure(self):
        result = prep_pe_nucleotide_quilt(**self._make_inputs())
        # Same shape as HDR quilt
        assert result['nuc_pct_df'].shape == snapshot((12, 6))
        assert result['mod_pct_df'].shape == snapshot((12, 6))
        # PE includes quantification window from first ref (not empty like HDR)
        assert result['quantification_window_idxs'] == snapshot([1, 2])

    def test_nuc_percentages_sum_to_one_per_batch(self):
        result = prep_pe_nucleotide_quilt(**self._make_inputs())
        nuc_df = result['nuc_pct_df']
        for batch in nuc_df['Batch'].unique():
            batch_df = nuc_df[nuc_df['Batch'] == batch]
            seq_cols = batch_df.columns[2:]
            col_sums = batch_df[seq_cols].sum(axis=0)
            np.testing.assert_allclose(col_sums.values, 1.0)


class TestPrepPeNucleotideQuiltAroundSgRNA:
    """prep_pe_nucleotide_quilt_around_sgRNA (plot_11b) — PE window slicing."""

    def _make_dfs(self, ref_names=None, ref_seq='ACGTACGTACGT'):
        """Build minimal multi-batch nuc/mod DataFrames matching 11a output shape."""
        import pandas as pd
        if ref_names is None:
            ref_names = ['Prime-edited', 'Scaffold-incorporated']
        n = len(ref_seq)
        nuc_rows = []
        for rn in ref_names:
            for nuc in ['A', 'C', 'G', 'T', 'N', '-']:
                row = [rn, nuc] + [0.0] * n
                nuc_rows.append(row)
        nuc_cols = ['Batch', 'Nucleotide'] + list(ref_seq)
        nuc_df = pd.DataFrame(nuc_rows, columns=nuc_cols)

        mod_rows = []
        for rn in ref_names:
            for mod in ['Insertions', 'Insertions_Left', 'Deletions', 'Substitutions', 'All_modifications', 'Total']:
                row = [rn, mod] + [0.0] * n
                mod_rows.append(row)
        mod_cols = ['Batch', 'Modification'] + list(ref_seq)
        mod_df = pd.DataFrame(mod_rows, columns=mod_cols)
        return nuc_df, mod_df

    def test_window_slicing(self):
        nuc_df, mod_df = self._make_dfs()
        result = prep_pe_nucleotide_quilt_around_sgRNA(
            nuc_df_for_plot=nuc_df,
            mod_df_for_plot=mod_df,
            cut_point=6,
            plot_half_window=3,
            ref_len=12,
            sgRNA_intervals=[(4, 8)],
            include_idxs_list=[5, 6, 7],
            sgRNA_names=['sg1'],
            sgRNA_mismatches=[0],
            sgRNA_sequences=['ACGT'],
            plot_root='/tmp/11b',
            save_also_png=False,
            custom_colors={},
        )
        assert list(result.keys()) == snapshot(
            [
                "nuc_pct_df",
                "mod_pct_df",
                "fig_filename_root",
                "save_also_png",
                "sgRNA_intervals",
                "sgRNA_names",
                "sgRNA_mismatches",
                "sgRNA_sequences",
                "quantification_window_idxs",
                "custom_colors",
            ]
        )
        # 2 batches × 6 nucs = 12 rows; window = 6 cols + 2 label cols = 8
        assert result['nuc_pct_df'].shape == snapshot((12, 8))
        assert result['sgRNA_intervals'] == snapshot([(0, 4)])
        assert result['quantification_window_idxs'] == snapshot([1, 2, 3])

    def test_clamped_near_start(self):
        nuc_df, mod_df = self._make_dfs()
        result = prep_pe_nucleotide_quilt_around_sgRNA(
            nuc_df_for_plot=nuc_df,
            mod_df_for_plot=mod_df,
            cut_point=1,
            plot_half_window=5,
            ref_len=12,
            sgRNA_intervals=[(0, 2)],
            include_idxs_list=[1],
            sgRNA_names=['sg1'],
            sgRNA_mismatches=[0],
            sgRNA_sequences=['ACGT'],
            plot_root='/tmp/11b',
            save_also_png=False,
            custom_colors={},
        )
        assert result['sgRNA_intervals'] == snapshot([(-2, 0)])
        assert result['quantification_window_idxs'] == snapshot([-1])


class TestPrepWindowedAlleles:
    """_prep_windowed_alleles — shared helper for plot_9 and plot_10h."""

    def _make_df(self):
        df = pd.DataFrame({
            'Reference_Sequence': ['ACGT', 'ACGT', 'ACGT'],
            'Unedited': [True, False, False],
            'n_deleted': [0, 1, 0],
            'n_inserted': [0, 0, 1],
            'n_mutated': [0, 0, 1],
            '#Reads': [60, 30, 10],
            '%Reads': [60.0, 30.0, 10.0],
        })
        df.index = pd.Index(['ACGT', 'ACCT', 'ACGN'], name='Aligned_Sequence')
        return df

    def test_returns_same_dataframe_object(self):
        """The returned df_alleles_around_cut IS the input (identity, not copy)."""
        df = self._make_df()
        result_df, _, _, _, _ = _prep_windowed_alleles(
            df_alleles_around_cut=df,
            cut_point=5,
            window_left=3,
            window_right=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
        )
        assert result_df is df

    def test_mutates_input_when_pct_adjustment_enabled(self):
        """When pct adjustment is on, the INPUT DataFrame is mutated in place."""
        df = self._make_df()
        original_pct_reads = df['%Reads'].tolist()
        _prep_windowed_alleles(
            df_alleles_around_cut=df,
            cut_point=5,
            window_left=3,
            window_right=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[],
            count_total=200,
            allele_plot_pcts_only_for_assigned_reference=True,
            expand_allele_plots_by_quantification=True,
        )
        # The original df should now have %AllReads (the old %Reads)
        assert df['%AllReads'].tolist() == original_pct_reads
        # And %Reads should be recalculated
        assert df['%Reads'].tolist() == [30.0, 15.0, 5.0]

    def test_does_not_mutate_when_pct_adjustment_disabled(self):
        """When pct adjustment is off, the DataFrame is not mutated."""
        df = self._make_df()
        original_pct_reads = df['%Reads'].tolist()
        _prep_windowed_alleles(
            df_alleles_around_cut=df,
            cut_point=5,
            window_left=3,
            window_right=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
        )
        assert '%AllReads' not in df.columns
        assert df['%Reads'].tolist() == original_pct_reads

    def test_ref_seq_slicing(self):
        _, _, ref_seq, _, _ = _prep_windowed_alleles(
            df_alleles_around_cut=self._make_df(),
            cut_point=5,
            window_left=3,
            window_right=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
        )
        # ref_sequence[5-3+1:5+3+1] = ref_sequence[3:9] = 'CGGTTA'
        assert ref_seq == 'CGGTTA'

    def test_asymmetric_window(self):
        _, _, ref_seq, _, _ = _prep_windowed_alleles(
            df_alleles_around_cut=self._make_df(),
            cut_point=5,
            window_left=2,
            window_right=4,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
        )
        # ref_sequence[5-2+1:5+4+1] = ref_sequence[4:10] = 'GGTTAC'
        assert ref_seq == 'GGTTAC'

    def test_sgRNA_interval_adjustment(self):
        _, _, _, new_intervals, new_sel_cols_start = _prep_windowed_alleles(
            df_alleles_around_cut=self._make_df(),
            cut_point=5,
            window_left=3,
            window_right=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[(3, 7)],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
        )
        # new_sel_cols_start = 5 - 3 = 2
        assert new_sel_cols_start == 2
        # new intervals = [(3-2-1, 7-2-1)] = [(0, 4)]
        assert new_intervals == [(0, 4)]

    def test_groupby_collapse(self):
        """Rows with same Aligned_Sequence+Reference_Sequence are collapsed."""
        df = pd.DataFrame({
            'Reference_Sequence': ['ACGT', 'ACGT', 'ACGT'],
            'Unedited': [True, False, False],
            'n_deleted': [0, 1, 0],
            'n_inserted': [0, 0, 1],
            'n_mutated': [0, 0, 1],
            '#Reads': [60, 30, 10],
            '%Reads': [60.0, 30.0, 10.0],
        })
        df.index = pd.Index(['ACGT', 'ACCT', 'ACGT'], name='Aligned_Sequence')
        _, df_to_plot, _, _, _ = _prep_windowed_alleles(
            df_alleles_around_cut=df,
            cut_point=5,
            window_left=3,
            window_right=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=False,
        )
        assert df_to_plot.shape[0] == 2
        assert df_to_plot['#Reads'].iloc[0] == 70


class TestPrepAllelesAroundCut:
    """prep_alleles_around_cut (plot_9) — pct adjustment + groupby + sgRNA intervals."""

    def _make_df(self):
        """Build a DataFrame matching output of get_dataframe_around_cut_asymmetrical."""
        df = pd.DataFrame({
            'Reference_Sequence': ['ACGT', 'ACGT', 'ACGT'],
            'Unedited': [True, False, False],
            'n_deleted': [0, 1, 0],
            'n_inserted': [0, 0, 1],
            'n_mutated': [0, 0, 1],
            '#Reads': [60, 30, 10],
            '%Reads': [60.0, 30.0, 10.0],
        })
        df.index = pd.Index(['ACGT', 'ACCT', 'ACGN'], name='Aligned_Sequence')
        return df

    def test_basic_structure(self):
        df = self._make_df()
        result = prep_alleles_around_cut(
            df_alleles_around_cut=df,
            cut_point=5,
            plot_half_window_left=3,
            plot_half_window_right=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[(3, 7)],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
        )
        assert 'df_alleles_around_cut' in result
        assert 'df_to_plot' in result
        assert 'ref_seq_around_cut' in result
        assert 'new_sgRNA_intervals' in result
        assert 'new_cut_point' in result
        # ref_seq_around_cut = ref_sequence[5-3+1:5+3+1] = ref_sequence[3:9] = 'CGGTTA'
        assert result['ref_seq_around_cut'] == 'CGGTTA'

    def test_returns_same_dataframe_object(self):
        """The returned df_alleles_around_cut IS the input (mutation contract)."""
        df = self._make_df()
        result = prep_alleles_around_cut(
            df_alleles_around_cut=df,
            cut_point=5,
            plot_half_window_left=3,
            plot_half_window_right=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[(3, 7)],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
        )
        assert result['df_alleles_around_cut'] is df

    def test_sgRNA_interval_recomputation(self):
        df = self._make_df()
        result = prep_alleles_around_cut(
            df_alleles_around_cut=df,
            cut_point=5,
            plot_half_window_left=3,
            plot_half_window_right=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[(3, 7)],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
        )
        # new_sel_cols_start = 5 - 3 = 2
        # new intervals = [(3-2-1, 7-2-1)] = [(0, 4)]
        assert result['new_sgRNA_intervals'] == [(0, 4)]
        # cut_point=5 is inside (3,7), so new_cut_point = 5-2-1 = 2
        assert result['new_cut_point'] == 2

    def test_cut_point_not_in_sgRNA(self):
        df = self._make_df()
        result = prep_alleles_around_cut(
            df_alleles_around_cut=df,
            cut_point=5,
            plot_half_window_left=3,
            plot_half_window_right=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[(0, 2)],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
        )
        assert result['new_cut_point'] is None


class TestPrepBaseEditQuilt:
    """prep_base_edit_quilt (plot_10h) — base edit DataFrame + x_labels."""

    def _make_df(self):
        df = pd.DataFrame({
            'Reference_Sequence': ['C-GT', 'C-GT'],
            'Unedited': [True, False],
            'n_deleted': [0, 1],
            'n_inserted': [0, 0],
            'n_mutated': [0, 0],
            '#Reads': [70, 30],
            '%Reads': [70.0, 30.0],
        })
        df.index = pd.Index(['C-GT', 'C-AT'], name='Aligned_Sequence')
        return df

    def test_basic_structure(self):
        df = self._make_df()
        result = prep_base_edit_quilt(
            df_alleles_around_cut=df,
            cut_point=5,
            plot_half_window=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[(3, 7)],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
            conversion_nuc_from='C',
        )
        assert 'df_alleles_around_cut' in result
        assert 'df_to_plot' in result
        assert 'ref_seq_around_cut' in result
        assert 'new_sgRNA_intervals' in result
        assert 'x_labels' in result
        # ref_seq_around_cut = ref_sequence[5-3+1:5+3+1] = ref_sequence[3:9] = 'CGGTTA'
        assert result['ref_seq_around_cut'] == 'CGGTTA'

    def test_returns_same_dataframe_object(self):
        """The returned df_alleles_around_cut IS the input (mutation contract)."""
        df = self._make_df()
        result = prep_base_edit_quilt(
            df_alleles_around_cut=df,
            cut_point=5,
            plot_half_window=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
            conversion_nuc_from='C',
        )
        assert result['df_alleles_around_cut'] is df

    def test_x_labels(self):
        df = self._make_df()
        result = prep_base_edit_quilt(
            df_alleles_around_cut=df,
            cut_point=5,
            plot_half_window=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
            conversion_nuc_from='C',
        )
        # C at index 2,3,9 → positions 3,4,10 (1-indexed)
        assert result['x_labels'] == [3, 4, 10]

    def test_sgRNA_interval_recomputation(self):
        df = self._make_df()
        result = prep_base_edit_quilt(
            df_alleles_around_cut=df,
            cut_point=5,
            plot_half_window=3,
            ref_sequence='AACCGGTTACGT',
            sgRNA_intervals=[(3, 7)],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
            conversion_nuc_from='C',
        )
        # new_sel_cols_start = 5 - 3 = 2
        # new intervals = [(3-2-1, 7-2-1)] = [(0, 4)]
        assert result['new_sgRNA_intervals'] == [(0, 4)]


class TestPrepAminoAcidTable:
    """prep_amino_acid_table (plot_9a) — amino acid conversion + cut point."""

    def test_amino_acid_cut_point(self):
        """Cut point is converted from nucleotide to amino acid coordinates."""
        # Empty df to avoid needing real alleles data
        df = pd.DataFrame(columns=[
            'Aligned_Sequence', 'Reference_Sequence', 'Read_Status',
            'n_deleted', 'n_inserted', 'n_mutated',
            '#Reads', '%Reads', 'ref_positions', 'Reference_Name',
        ])
        result = prep_amino_acid_table(
            coding_seq='ATGATGATG',  # 3 codons = 3 amino acids
            cut_point=10,
            exon_positions_start=2,
            df_alleles_for_ref=df,
            exon_interval_start=0,
            blosum_path='/tmp/BLOSUM62',
        )
        # amino_acid_cut_point = (10 - 2 + 1) // 3 = 3
        assert result['amino_acid_cut_point'] == 3

    def test_coding_seq_amino_acids(self):
        """Nucleotide coding sequence is converted to amino acids."""
        df = pd.DataFrame(columns=[
            'Aligned_Sequence', 'Reference_Sequence', 'Read_Status',
            'n_deleted', 'n_inserted', 'n_mutated',
            '#Reads', '%Reads', 'ref_positions', 'Reference_Name',
        ])
        result = prep_amino_acid_table(
            coding_seq='ATGATGATG',  # Met-Met-Met
            cut_point=5,
            exon_positions_start=0,
            df_alleles_for_ref=df,
            exon_interval_start=0,
            blosum_path='/tmp/BLOSUM62',
        )
        assert result['coding_seq_amino_acids'] == 'MMM'

    def test_return_keys(self):
        df = pd.DataFrame(columns=[
            'Aligned_Sequence', 'Reference_Sequence', 'Read_Status',
            'n_deleted', 'n_inserted', 'n_mutated',
            '#Reads', '%Reads', 'ref_positions', 'Reference_Name',
        ])
        result = prep_amino_acid_table(
            coding_seq='ATGATG',
            cut_point=3,
            exon_positions_start=0,
            df_alleles_for_ref=df,
            exon_interval_start=0,
            blosum_path='/tmp/BLOSUM62',
        )
        assert 'coding_seq_amino_acids' in result
        assert 'amino_acid_cut_point' in result
        assert 'df_to_plot' in result
