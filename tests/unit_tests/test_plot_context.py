"""Unit tests for PlotContext dataclass."""

import pytest
from types import SimpleNamespace


# =============================================================================
# Construction
# =============================================================================


class TestPlotContextConstruction:
    """Test PlotContext can be constructed and fields are accessible."""

    def test_minimal_construction(self):
        """PlotContext can be constructed with required fields."""
        from CRISPResso2.plots.plot_context import PlotContext

        args = SimpleNamespace(name='test')
        ctx = PlotContext(
            args=args,
            run_data={'results': {}},
            refs={'ref1': {}},
            ref_names=['ref1'],
            counts_total={'ref1': 100},
            counts_modified={'ref1': 50},
            counts_unmodified={'ref1': 50},
            counts_discarded={'ref1': 0},
            counts_insertion={'ref1': 10},
            counts_deletion={'ref1': 20},
            counts_substitution={'ref1': 5},
            class_counts={'Reference': 50, 'NHEJ': 50},
            N_TOTAL=100,
            df_alleles=None,
            all_insertion_count_vectors={'ref1': [0, 1, 2]},
            all_insertion_left_count_vectors={'ref1': [0, 0, 1]},
            all_deletion_count_vectors={'ref1': [0, 2, 0]},
            all_substitution_count_vectors={'ref1': [1, 0, 0]},
            all_indelsub_count_vectors={'ref1': [1, 2, 2]},
            all_substitution_base_vectors={'ref1': {}},
            all_base_count_vectors={'ref1': {}},
            insertion_count_vectors={'ref1': [0, 1]},
            deletion_count_vectors={'ref1': [0, 1]},
            substitution_count_vectors={'ref1': [1, 0]},
            insertion_length_vectors={'ref1': [0, 1]},
            deletion_length_vectors={'ref1': [0, 1]},
            nucleotide_frequency_summary={'ref1': None},
            nucleotide_percentage_summary={'ref1': None},
            hists_frameshift={'ref1': {}},
            hists_inframe={'ref1': {}},
            counts_modified_frameshift={'ref1': 0},
            counts_modified_non_frameshift={'ref1': 0},
            counts_non_modified_non_frameshift={'ref1': 0},
            counts_splicing_sites_modified={'ref1': 0},
        )
        assert ctx.args.name == 'test'
        assert ctx.ref_names == ['ref1']
        assert ctx.N_TOTAL == 100

    def test_defaults(self):
        """Optional fields have correct defaults."""
        from CRISPResso2.plots.plot_context import PlotContext

        ctx = PlotContext(
            args=SimpleNamespace(),
            run_data={},
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
            nucleotide_frequency_summary={},
            nucleotide_percentage_summary={},
            hists_frameshift={},
            hists_inframe={},
            counts_modified_frameshift={},
            counts_modified_non_frameshift={},
            counts_non_modified_non_frameshift={},
            counts_splicing_sites_modified={},
        )
        # HDR vectors default to empty dicts
        assert ctx.ref1_all_insertion_count_vectors == {}
        assert ctx.ref1_all_deletion_count_vectors == {}
        assert ctx.ref1_all_substitution_count_vectors == {}
        assert ctx.ref1_all_indelsub_count_vectors == {}
        assert ctx.ref1_all_insertion_left_count_vectors == {}
        assert ctx.ref1_all_base_count_vectors == {}
        # Config defaults
        assert ctx.custom_config == {}
        assert ctx._jp is None
        assert ctx.save_png is False
        assert ctx.output_directory == ""
        # Scope fields default to None
        assert ctx.ref_name is None
        assert ctx.sgRNA_ind is None

    def test_scope_fields_mutable(self):
        """ref_name and sgRNA_ind can be set after construction."""
        from CRISPResso2.plots.plot_context import PlotContext

        ctx = PlotContext(
            args=SimpleNamespace(),
            run_data={},
            refs={'amp1': {}, 'amp2': {}},
            ref_names=['amp1', 'amp2'],
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
            nucleotide_frequency_summary={},
            nucleotide_percentage_summary={},
            hists_frameshift={},
            hists_inframe={},
            counts_modified_frameshift={},
            counts_modified_non_frameshift={},
            counts_non_modified_non_frameshift={},
            counts_splicing_sites_modified={},
        )
        ctx.ref_name = 'amp1'
        ctx.sgRNA_ind = 0
        assert ctx.ref_name == 'amp1'
        assert ctx.sgRNA_ind == 0

        ctx.ref_name = 'amp2'
        ctx.sgRNA_ind = 1
        assert ctx.ref_name == 'amp2'
        assert ctx.sgRNA_ind == 1

    def test_references_not_copied(self):
        """PlotContext holds references, not copies, of data structures."""
        from CRISPResso2.plots.plot_context import PlotContext

        counts = {'ref1': 100}
        ctx = PlotContext(
            args=SimpleNamespace(),
            run_data={},
            refs={},
            ref_names=[],
            counts_total=counts,
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
            nucleotide_frequency_summary={},
            nucleotide_percentage_summary={},
            hists_frameshift={},
            hists_inframe={},
            counts_modified_frameshift={},
            counts_modified_non_frameshift={},
            counts_non_modified_non_frameshift={},
            counts_splicing_sites_modified={},
        )
        # Mutate original — PlotContext sees it (zero-copy)
        counts['ref1'] = 200
        assert ctx.counts_total['ref1'] == 200


class TestPlotContextImport:
    """Test that PlotContext is importable from the plotting package."""

    def test_import_from_module(self):
        from CRISPResso2.plots.plot_context import PlotContext
        assert PlotContext is not None

    def test_import_from_module(self):
        from CRISPResso2.plots.plot_context import PlotContext
        assert PlotContext is not None
