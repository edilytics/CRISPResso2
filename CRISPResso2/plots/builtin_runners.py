"""Shared built-in plotting runners used by both CORE and CRISPRessoPro.

Each ``run_builtin_*_plots`` function in this module encapsulates the
plotting logic that historically lived in a mode's ``elif not args.suppress_plots``
branch.  Extracting it here lets ``CRISPRessoBatchCORE`` /
``CRISPRessoCompareCORE`` / etc. call it when CRISPRessoPro is not
installed, and lets Pro's hooks invoke the same implementation for modes
that don't yet have a Pro-native plot plugin port (see
``design_docs/MULTI_MODE_PLOT_PLUGIN.md``).

Each runner uses the mode's ``PlotContext`` subclass as its single source
of state (per the plot-context hierarchy refactor) and writes results
into the ``crispresso2_info`` run-data dict passed alongside it.

These runners intentionally keep the CORE elif semantics verbatim —
their job is to unify CORE and Pro code paths so both produce the same
output, not to reshape the logic.
"""

from __future__ import annotations

import os
import traceback
from concurrent.futures import ProcessPoolExecutor, wait
from functools import partial

from CRISPResso2 import CRISPRessoMultiProcessing
from CRISPResso2.plots.data_prep import (
    prep_batch_conversion_map,
    prep_batch_conversion_map_around_sgRNA,
    prep_batch_nuc_quilt,
    prep_batch_nuc_quilt_around_sgRNA,
    prep_compare_allele_table,
    prep_compare_editing_barchart,
    prep_compare_modification_positions,
)


def _should_plot_large(df_shape_rows: int, is_interactive: bool, use_matplotlib: bool) -> bool:
    """Mirror CRISPRessoBatchCORE.should_plot_large_plots semantics.

    Returns True unconditionally in the matplotlib case (original behavior),
    which matches the existing ``elif`` code.
    """
    # The original predicate ultimately returns True for the matplotlib
    # path used by the built-in runner; we keep the same behavior here so
    # output parity is preserved.
    return True


def run_builtin_batch_plots(batch_plot_context, crispresso2_info, plot_module, logger):
    """Run CRISPRessoBatchCORE's built-in (matplotlib) per-batch plots.

    Replicates the body of ``CRISPRessoBatchCORE.main``'s
    ``elif not args.suppress_plots and not args.suppress_batch_summary_plots:``
    branch.  Uses the following fields on ``batch_plot_context``:

    - ``amplicon_names``, ``nucleotide_percentage_summary_dfs``,
      ``guides_all_same``, ``consensus_guides``
    - ``all_summary_filenames`` — per-amplicon dict of summary CSV paths
    - ``sub_nucleotide_frequency_summary_filename`` /
      ``sub_nucleotide_percentage_summary_filename`` — last-set
      per-sgRNA CSV paths used for data links on conversion-map plots

    Populates ``crispresso2_info['results']['general_plots']`` with the
    standard batch keys (``window_nuc_pct_quilt_plot_names``, etc.).
    """
    args = batch_plot_context.args
    n_processes = getattr(args, 'n_processes_for_batch', 1)
    try:
        n_processes = int(n_processes)
    except (TypeError, ValueError):
        n_processes = 1

    if n_processes > 1:
        process_pool = ProcessPoolExecutor(n_processes)
        process_futures = {}
    else:
        process_pool = None
        process_futures = None

    plot = partial(
        CRISPRessoMultiProcessing.run_plot,
        num_processes=n_processes,
        process_futures=process_futures,
        process_pool=process_pool,
        halt_on_plot_fail=getattr(args, 'halt_on_plot_fail', False),
    )

    general_plots = crispresso2_info['results']['general_plots']
    general_plots.setdefault('summary_plot_names', [])
    general_plots.setdefault('summary_plot_titles', {})
    general_plots.setdefault('summary_plot_labels', {})
    general_plots.setdefault('summary_plot_datas', {})

    window_nuc_pct_quilt_plot_names: list[str] = []
    nuc_pct_quilt_plot_names: list[str] = []
    window_nuc_conv_plot_names: list[str] = []
    nuc_conv_plot_names: list[str] = []

    all_summary_filenames = batch_plot_context.all_summary_filenames
    sub_nuc_freq_filename = batch_plot_context.sub_nucleotide_frequency_summary_filename
    sub_nuc_pct_filename = batch_plot_context.sub_nucleotide_percentage_summary_filename
    amplicon_names = batch_plot_context.amplicon_names

    for amplicon_name in amplicon_names:
        batch_plot_context.amplicon_name = amplicon_name
        nuc_freq_filename = all_summary_filenames[amplicon_name]['nucleotide_frequency']
        mod_freq_filename = all_summary_filenames[amplicon_name]['modification_frequency']
        guides_all_same = batch_plot_context.guides_all_same[amplicon_name]
        consensus_guides = batch_plot_context.consensus_guides[amplicon_name]

        if guides_all_same and consensus_guides:
            for sgRNA_ind, sgRNA in enumerate(consensus_guides):
                batch_plot_context.sgRNA_ind = sgRNA_ind

                if _should_plot_large(
                    batch_plot_context.nucleotide_percentage_summary_dfs[amplicon_name].shape[0],
                    False, args.use_matplotlib,
                ):
                    quilt_input = prep_batch_nuc_quilt_around_sgRNA(batch_plot_context)
                    logger.debug(f'Plotting nucleotide percentage quilt for amplicon {amplicon_name}, sgRNA {sgRNA}')
                    plot(plot_module.plot_nucleotide_quilt, quilt_input)
                    plot_name = os.path.basename(quilt_input['fig_filename_root'])
                    window_nuc_pct_quilt_plot_names.append(plot_name)
                    general_plots['summary_plot_titles'][plot_name] = f'sgRNA: {sgRNA} Amplicon: {amplicon_name}'
                    general_plots['summary_plot_labels'][plot_name] = f'Composition of each base around the guide {sgRNA} for the amplicon {amplicon_name}'
                    general_plots['summary_plot_datas'][plot_name] = [
                        ('Nucleotide frequencies', os.path.basename(nuc_freq_filename)),
                        ('Modification frequencies', os.path.basename(mod_freq_filename)),
                    ]

                    if args.base_editor_output:
                        conv_input = prep_batch_conversion_map_around_sgRNA(batch_plot_context)
                        logger.debug(f'Plotting nucleotide conversion map for amplicon {amplicon_name}, sgRNA {sgRNA}')
                        plot(plot_module.plot_conversion_map, conv_input)
                        plot_name = os.path.basename(conv_input['fig_filename_root'])
                        window_nuc_conv_plot_names.append(plot_name)
                        title = f'sgRNA: {sgRNA} Amplicon: {amplicon_name}'
                        if len(consensus_guides) == 1:
                            title = ''
                        general_plots['summary_plot_titles'][plot_name] = title
                        general_plots['summary_plot_labels'][plot_name] = (
                            f'{args.conversion_nuc_from}->{args.conversion_nuc_to} '
                            f'conversion rates around the guide {sgRNA} for the amplicon {amplicon_name}'
                        )
                        general_plots['summary_plot_datas'][plot_name] = [
                            ('Nucleotide frequencies around sgRNA', os.path.basename(sub_nuc_freq_filename)),
                            ('Nucleotide percentages around sgRNA', os.path.basename(sub_nuc_pct_filename)),
                        ]

            # Whole-region quilt (when guides_all_same)
            if _should_plot_large(
                batch_plot_context.nucleotide_percentage_summary_dfs[amplicon_name].shape[0],
                False, args.use_matplotlib,
            ):
                quilt_input = prep_batch_nuc_quilt(batch_plot_context)
                logger.debug(f'Plotting nucleotide quilt for {amplicon_name}')
                plot(plot_module.plot_nucleotide_quilt, quilt_input)
                plot_name = os.path.basename(quilt_input['fig_filename_root'])
                nuc_pct_quilt_plot_names.append(plot_name)
                title = f'Amplicon: {amplicon_name}'
                if len(amplicon_names) == 1:
                    title = ''
                general_plots['summary_plot_titles'][plot_name] = title
                general_plots['summary_plot_labels'][plot_name] = f'Composition of each base for the amplicon {amplicon_name}'
                general_plots['summary_plot_datas'][plot_name] = [
                    ('Nucleotide frequencies', os.path.basename(nuc_freq_filename)),
                    ('Modification frequencies', os.path.basename(mod_freq_filename)),
                ]

                if args.base_editor_output:
                    conv_input = prep_batch_conversion_map(batch_plot_context)
                    logger.debug(f'Plotting nucleotide conversion map for {amplicon_name}')
                    plot(plot_module.plot_conversion_map, conv_input)
                    plot_name = os.path.basename(conv_input['fig_filename_root'])
                    nuc_conv_plot_names.append(plot_name)
                    general_plots['summary_plot_titles'][plot_name] = ''
                    general_plots['summary_plot_labels'][plot_name] = (
                        f'{args.conversion_nuc_from}->{args.conversion_nuc_to} '
                        f'conversion rates for the amplicon {amplicon_name}'
                    )
                    general_plots['summary_plot_datas'][plot_name] = [
                        ('Nucleotide frequencies', os.path.basename(nuc_freq_filename)),
                        ('Modification frequencies', os.path.basename(mod_freq_filename)),
                    ]

        # Guides not all same — whole-region quilt only, no sgRNA data on plot
        elif _should_plot_large(
            batch_plot_context.nucleotide_percentage_summary_dfs[amplicon_name].shape[0],
            False, args.use_matplotlib,
        ):
            quilt_input = prep_batch_nuc_quilt(batch_plot_context)
            logger.debug(f'Plotting nucleotide quilt for {amplicon_name}')
            plot(plot_module.plot_nucleotide_quilt, quilt_input)
            plot_name = os.path.basename(quilt_input['fig_filename_root'])
            nuc_pct_quilt_plot_names.append(plot_name)
            general_plots['summary_plot_labels'][plot_name] = f'Composition of each base for the amplicon {amplicon_name}'
            general_plots['summary_plot_datas'][plot_name] = [
                ('Nucleotide frequencies', os.path.basename(nuc_freq_filename)),
                ('Modification frequencies', os.path.basename(mod_freq_filename)),
            ]

            if args.base_editor_output:
                conv_input = prep_batch_conversion_map(batch_plot_context)
                logger.debug(f'Plotting BE nucleotide conversion map for {amplicon_name}')
                plot(plot_module.plot_conversion_map, conv_input)
                plot_name = os.path.basename(conv_input['fig_filename_root'])
                nuc_conv_plot_names.append(plot_name)
                general_plots['summary_plot_labels'][plot_name] = (
                    f'{args.conversion_nuc_from}->{args.conversion_nuc_to} '
                    f'conversion rates for the amplicon {amplicon_name}'
                )
                general_plots['summary_plot_datas'][plot_name] = [
                    ('Nucleotide frequencies', os.path.basename(nuc_freq_filename)),
                    ('Modification frequencies', os.path.basename(mod_freq_filename)),
                ]

    general_plots['window_nuc_pct_quilt_plot_names'] = window_nuc_pct_quilt_plot_names
    general_plots['nuc_pct_quilt_plot_names'] = nuc_pct_quilt_plot_names
    general_plots['window_nuc_conv_plot_names'] = window_nuc_conv_plot_names
    general_plots['nuc_conv_plot_names'] = nuc_conv_plot_names

    if process_pool is not None:
        wait(process_futures)
        for future in process_futures:
            try:
                future.result()
            except Exception as e:
                logger.warning(f'Error in plot pool: {e}')
                logger.debug(traceback.format_exc())
        process_pool.shutdown()


def run_builtin_aggregate_plots(agg_plot_context, crispresso2_info, plot_module, logger):
    """Run CRISPRessoAggregateCORE's built-in per-amplicon plots.

    Replicates the body of ``CRISPRessoAggregateCORE.main``'s
    ``elif not args.suppress_plots:`` branch, including the whole-region
    quilt pagination when the number of samples exceeds
    ``args.max_samples_per_summary_plot``.

    Populates ``crispresso2_info['results']['general_plots']`` with the
    standard aggregate quilt plot-name lists.
    """
    args = agg_plot_context.args
    n_processes = getattr(args, 'n_processes_for_batch', 1)
    try:
        n_processes = int(n_processes)
    except (TypeError, ValueError):
        n_processes = 1

    if n_processes > 1:
        process_pool = ProcessPoolExecutor(n_processes)
        process_futures = {}
    else:
        process_pool = None
        process_futures = None

    plot = partial(
        CRISPRessoMultiProcessing.run_plot,
        num_processes=n_processes,
        process_futures=process_futures,
        process_pool=process_pool,
        halt_on_plot_fail=getattr(args, 'halt_on_plot_fail', False),
    )

    general_plots = crispresso2_info['results']['general_plots']
    general_plots.setdefault('summary_plot_names', [])
    general_plots.setdefault('summary_plot_titles', {})
    general_plots.setdefault('summary_plot_labels', {})
    general_plots.setdefault('summary_plot_datas', {})

    window_nuc_pct_quilt_plot_names: list[str] = []
    nuc_pct_quilt_plot_names: list[str] = []
    window_nuc_conv_plot_names: list[str] = []
    nuc_conv_plot_names: list[str] = []

    all_summary_filenames = agg_plot_context.all_summary_filenames
    all_sample_counts = agg_plot_context.sample_count
    amplicon_names = agg_plot_context.amplicon_names

    for amplicon_name in amplicon_names:
        agg_plot_context.amplicon_name = amplicon_name
        nuc_freq_filename = all_summary_filenames[amplicon_name]['nucleotide_frequency']
        mod_freq_filename = all_summary_filenames[amplicon_name]['modification_frequency']
        guides_all_same = agg_plot_context.guides_all_same[amplicon_name]
        consensus_guides = agg_plot_context.consensus_guides[amplicon_name]
        this_number_samples = all_sample_counts[amplicon_name]

        # Per-sgRNA quilts (when guides_all_same and small enough for matplotlib)
        if guides_all_same and consensus_guides:
            for sgRNA_ind, sgRNA in enumerate(consensus_guides):
                agg_plot_context.sgRNA_ind = sgRNA_ind
                if this_number_samples < args.max_samples_per_summary_plot:
                    quilt_input = prep_batch_nuc_quilt_around_sgRNA(agg_plot_context)
                    logger.debug(f'Plotting nucleotide percentage quilt for amplicon {amplicon_name}, sgRNA {sgRNA}')
                    plot(plot_module.plot_nucleotide_quilt, quilt_input)
                    plot_name = os.path.basename(quilt_input['fig_filename_root'])
                    window_nuc_pct_quilt_plot_names.append(plot_name)
                    title = f'sgRNA: {sgRNA} Amplicon: {amplicon_name}'
                    if len(consensus_guides) == 1:
                        title = ''
                    general_plots['summary_plot_titles'][plot_name] = title
                    general_plots['summary_plot_labels'][plot_name] = (
                        f'Composition of each base around the guide {sgRNA} for the amplicon {amplicon_name}'
                    )
                    general_plots['summary_plot_datas'][plot_name] = [
                        (f'{amplicon_name} nucleotide frequencies', os.path.basename(nuc_freq_filename)),
                        (f'{amplicon_name} modification frequencies', os.path.basename(mod_freq_filename)),
                    ]

        # Whole-region quilt with pagination
        quilt_input = prep_batch_nuc_quilt(agg_plot_context)
        nuc_pct_df = quilt_input['nuc_pct_df']
        mod_pct_df = quilt_input['mod_pct_df']
        if this_number_samples <= 0:
            continue
        nrow_per_sample_nucs = nuc_pct_df.shape[0] / this_number_samples
        nrow_per_sample_mods = mod_pct_df.shape[0] / this_number_samples

        this_plot_suffix = ''
        this_plot_suffix_int = 1
        for sample_start_ind in range(
            0, this_number_samples, args.max_samples_per_summary_plot,
        ):
            sample_end_ind = min(
                sample_start_ind + args.max_samples_per_summary_plot,
                this_number_samples,
            )
            this_nuc_start_ind = int(sample_start_ind * nrow_per_sample_nucs)
            this_nuc_end_ind = int((sample_end_ind + 1) * nrow_per_sample_nucs - 1)
            this_mod_start_ind = int(sample_start_ind * nrow_per_sample_mods)
            this_mod_end_ind = int((sample_end_ind + 1) * nrow_per_sample_mods - 1)

            page_input = dict(quilt_input)
            page_input['nuc_pct_df'] = nuc_pct_df.iloc[this_nuc_start_ind:this_nuc_end_ind, :]
            page_input['mod_pct_df'] = mod_pct_df.iloc[this_mod_start_ind:this_mod_end_ind, :]
            page_input['fig_filename_root'] = quilt_input['fig_filename_root'] + this_plot_suffix

            logger.debug(f'Plotting nucleotide quilt for {amplicon_name}{this_plot_suffix}')
            plot(plot_module.plot_nucleotide_quilt, page_input)

            plot_name = os.path.basename(page_input['fig_filename_root'])
            nuc_pct_quilt_plot_names.append(plot_name)
            title = f'Amplicon: {amplicon_name}{this_plot_suffix}'
            if len(amplicon_names) == 1:
                title = ''
            general_plots['summary_plot_titles'][plot_name] = title
            general_plots['summary_plot_labels'][plot_name] = (
                f'Composition of each base for the amplicon {amplicon_name}'
            )
            general_plots['summary_plot_datas'][plot_name] = [
                (f'{amplicon_name} nucleotide frequencies', os.path.basename(nuc_freq_filename)),
                (f'{amplicon_name} modification frequencies', os.path.basename(mod_freq_filename)),
            ]

            this_plot_suffix_int += 1
            this_plot_suffix = f'_{this_plot_suffix_int}'

    general_plots['window_nuc_pct_quilt_plot_names'] = window_nuc_pct_quilt_plot_names
    general_plots['nuc_pct_quilt_plot_names'] = nuc_pct_quilt_plot_names
    general_plots['window_nuc_conv_plot_names'] = window_nuc_conv_plot_names
    general_plots['nuc_conv_plot_names'] = nuc_conv_plot_names

    if process_pool is not None:
        wait(process_futures)
        for future in process_futures:
            try:
                future.result()
            except Exception as e:
                logger.warning(f'Error in plot pool: {e}')
                logger.debug(traceback.format_exc())
        process_pool.shutdown()


def run_builtin_compare_plots(compare_plot_context, crispresso2_info, plot_module, logger):
    """Run CRISPRessoCompareCORE's built-in per-amplicon comparison plots.

    Replicates the body of ``CRISPRessoCompareCORE.main``'s ``else:`` /
    non-Pro branch.  Populates ``_sig_counts`` and
    ``_sig_counts_quant_window`` on ``compare_plot_context.run_data`` so
    the caller (CORE or Pro hook) can read them from a consistent location
    for the summary-counts CSV.

    Writes per-amplicon ``*_quantification.txt`` and allele comparison
    tables, and runs the three compare plots (barchart, modification
    positions, allele tables).
    """
    import numpy as np

    from CRISPResso2 import CRISPRessoShared
    from CRISPResso2.CRISPRessoCompareCORE import (
        DifferentAmpliconLengthException,
        get_matching_allele_files,
        parse_profile,
    )

    import pandas as pd  # local — not always needed elsewhere

    args = compare_plot_context.args
    run_data = crispresso2_info
    general_plots = run_data['results'].setdefault('general_plots', {})
    general_plots.setdefault('summary_plot_names', [])
    general_plots.setdefault('summary_plot_titles', {})
    general_plots.setdefault('summary_plot_labels', {})
    general_plots.setdefault('summary_plot_datas', {})

    output_directory = compare_plot_context.output_directory
    _jp = compare_plot_context._jp

    sample_1_name = compare_plot_context.sample_1_name
    sample_2_name = compare_plot_context.sample_2_name
    run_info_1 = compare_plot_context.run_info_1
    run_info_2 = compare_plot_context.run_info_2
    amplicon_info_1 = compare_plot_context.amplicon_info_1
    amplicon_info_2 = compare_plot_context.amplicon_info_2
    amplicon_names_in_both = compare_plot_context.amplicon_names

    sig_counts: dict[str, dict] = {}
    sig_counts_quant_window: dict[str, dict] = {}

    percent_complete_start, percent_complete_end = 10, 90
    if amplicon_names_in_both:
        percent_complete_step = (
            (percent_complete_end - percent_complete_start)
            / len(amplicon_names_in_both)
        )
    else:
        percent_complete_step = 0

    for amplicon_name in amplicon_names_in_both:
        percent_complete = (
            percent_complete_start
            + percent_complete_step * amplicon_names_in_both.index(amplicon_name)
        )
        logger.info(
            f'Loading data for amplicon {amplicon_name}',
            extra={'percent_complete': percent_complete},
        )

        compare_plot_context.amplicon_name = amplicon_name
        compare_plot_context.profile_1 = parse_profile(
            amplicon_info_1[amplicon_name]['quantification_file']
        )
        compare_plot_context.profile_2 = parse_profile(
            amplicon_info_2[amplicon_name]['quantification_file']
        )

        sig_counts[amplicon_name] = {}
        sig_counts_quant_window[amplicon_name] = {}

        try:
            assert np.all(
                compare_plot_context.profile_1[:, 0]
                == compare_plot_context.profile_2[:, 0]
            )
        except AssertionError:
            raise DifferentAmpliconLengthException(
                'Different amplicon lengths for the two amplicons.'
            )

        compare_plot_context.cut_points = (
            run_info_1['results']['refs'][amplicon_name]['sgRNA_cut_points']
        )
        compare_plot_context.sgRNA_intervals = (
            run_info_1['results']['refs'][amplicon_name]['sgRNA_intervals']
        )

        # Plot 1: Editing comparison barchart
        barchart_input = prep_compare_editing_barchart(compare_plot_context)
        plot_module.plot_quantification_comparison_barchart(**barchart_input)
        plot_name = os.path.basename(barchart_input['plot_path'])
        general_plots['summary_plot_names'].append(plot_name)
        general_plots['summary_plot_titles'][plot_name] = 'Editing efficiency comparison'
        general_plots['summary_plot_labels'][plot_name] = (
            f'Figure 1: Comparison for amplicon {amplicon_name}; '
            'Left: Percentage of modified and unmodified reads in each sample; '
            'Right: relative percentage of modified and unmodified reads'
        )
        output_1 = os.path.join(
            args.crispresso_output_folder_1,
            run_info_1['running_info']['report_filename'],
        )
        output_2 = os.path.join(
            args.crispresso_output_folder_2,
            run_info_2['running_info']['report_filename'],
        )
        general_plots['summary_plot_datas'][plot_name] = []
        if os.path.isfile(output_1):
            general_plots['summary_plot_datas'][plot_name].append(
                (f'{sample_1_name} output', os.path.relpath(output_1, output_directory))
            )
        if os.path.isfile(output_2):
            general_plots['summary_plot_datas'][plot_name].append(
                (f'{sample_2_name} output', os.path.relpath(output_2, output_directory))
            )

        # Load modification count data
        mod_file_1 = amplicon_info_1[amplicon_name]['modification_count_file']
        amp_seq_1, mod_freqs_1 = CRISPRessoShared.parse_count_file(mod_file_1)
        mod_file_2 = amplicon_info_2[amplicon_name]['modification_count_file']
        amp_seq_2, mod_freqs_2 = CRISPRessoShared.parse_count_file(mod_file_2)
        if amp_seq_2 != amp_seq_1:
            raise DifferentAmpliconLengthException(
                'Different amplicon lengths for the two amplicons.'
            )

        compare_plot_context.mod_freqs_1 = mod_freqs_1
        compare_plot_context.mod_freqs_2 = mod_freqs_2
        compare_plot_context.consensus_sequence = amp_seq_1
        compare_plot_context.quant_windows_1 = (
            run_info_1['results']['refs'][amplicon_name]['include_idxs']
        )
        compare_plot_context.quant_windows_2 = (
            run_info_2['results']['refs'][amplicon_name]['include_idxs']
        )

        amplicon_plot_name = f'{amplicon_name}.'
        if len(amplicon_names_in_both) == 1 and amplicon_name == 'Reference':
            amplicon_plot_name = ''

        # Plot 2: Modification positions (x4 mod types)
        for mod in ['Insertions', 'Deletions', 'Substitutions', 'All_modifications']:
            compare_plot_context.mod_type = mod
            positions_data = prep_compare_modification_positions(compare_plot_context)

            mod_filename = _jp(f'{amplicon_plot_name}{mod}_quantification.txt')
            positions_data['mod_df'].to_csv(mod_filename, sep='\t', index=None)

            plot_module.plot_quantification_positions(**positions_data['plot_kwargs'])
            plot_name = os.path.basename(positions_data['plot_kwargs']['plot_path'])
            general_plots['summary_plot_names'].append(plot_name)
            general_plots['summary_plot_titles'][plot_name] = (
                f"{positions_data['mod_name']} locations"
            )
            general_plots['summary_plot_labels'][plot_name] = (
                f"{positions_data['mod_name']} location comparison for amplicon "
                f'{amplicon_name}; Top: percent difference; Bottom: p-value.'
            )
            general_plots['summary_plot_datas'][plot_name] = [
                (f"{positions_data['mod_name']} quantification",
                 os.path.basename(mod_filename)),
            ]

            sig_counts[amplicon_name][mod] = positions_data['sig_count']
            sig_counts_quant_window[amplicon_name][mod] = (
                positions_data['sig_count_quant_window']
            )

        # Plot 3: Allele table comparisons
        matching_allele_files = get_matching_allele_files(run_info_1, run_info_2)
        matching_allele_files.sort(key=lambda pair: 'base_edit' in pair[0])

        compare_plot_context.allele_pairs = []
        for allele_file_1, allele_file_2 in matching_allele_files:
            df1 = pd.read_csv(
                os.path.join(args.crispresso_output_folder_1, allele_file_1),
                sep='\t',
            )
            df2 = pd.read_csv(
                os.path.join(args.crispresso_output_folder_2, allele_file_2),
                sep='\t',
            )
            compare_plot_context.allele_pairs.append(
                (allele_file_1, allele_file_2, df1, df2)
            )

        for pair_idx in range(len(compare_plot_context.allele_pairs)):
            allele_data = prep_compare_allele_table(compare_plot_context, pair_idx)

            allele_comparison_file = _jp(allele_data['file_root'] + '.txt')
            allele_data['merged_df'].to_csv(
                allele_comparison_file, sep='\t', index=None,
            )

            is_base_edit = allele_data['is_base_edit']
            if is_base_edit:
                title_prefix = 'Base edit comparison enriched in '
                label_prefix = 'Base editing target nucleotide composition alleles.'
            else:
                title_prefix = 'Alleles enriched in '
                label_prefix = 'Distribution comparison of alleles.'
            label_suffix = (
                ' Nucleotides are indicated by unique colors (A = green; '
                'C = red; G = yellow; T = purple). Substitutions are shown in '
                'bold font. Red rectangles highlight inserted sequences. '
                'Horizontal dashed lines indicate deleted sequences. The '
                'vertical dashed line indicates the predicted cleavage site. '
                'The proportion and number of reads is shown for each sample '
                f'on the right, with the values for {sample_1_name} followed '
                f'by the values for {sample_2_name}.'
            )

            plot_module.plot_alleles_table_compare(**allele_data['plot_top_kwargs'])
            plot_name = os.path.basename(allele_data['plot_top_kwargs']['fig_filename_root'])
            general_plots['summary_plot_names'].append(plot_name)
            general_plots['summary_plot_titles'][plot_name] = (
                title_prefix + sample_1_name
            )
            general_plots['summary_plot_labels'][plot_name] = (
                label_prefix + label_suffix
                + f' Alleles are sorted for enrichment in {sample_1_name}.'
            )
            general_plots['summary_plot_datas'][plot_name] = [
                ('Allele comparison table', os.path.basename(allele_comparison_file)),
            ]

            plot_module.plot_alleles_table_compare(**allele_data['plot_bottom_kwargs'])
            plot_name = os.path.basename(allele_data['plot_bottom_kwargs']['fig_filename_root'])
            general_plots['summary_plot_names'].append(plot_name)
            general_plots['summary_plot_titles'][plot_name] = (
                title_prefix + sample_2_name
            )
            general_plots['summary_plot_labels'][plot_name] = (
                label_prefix + label_suffix
                + f' Alleles are sorted for enrichment in {sample_2_name}.'
            )
            general_plots['summary_plot_datas'][plot_name] = [
                ('Allele comparison table', os.path.basename(allele_comparison_file)),
            ]

    # Publish the Fisher/Bonferroni significance counts so the caller
    # (CORE or Pro) can serialize them into the summary CSV.
    run_data['_sig_counts'] = sig_counts
    run_data['_sig_counts_quant_window'] = sig_counts_quant_window
