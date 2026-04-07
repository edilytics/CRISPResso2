# -*- coding: utf-8 -*-
"""CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
"""
import os
from copy import deepcopy
import sys
import traceback
from CRISPResso2 import CRISPRessoShared
from CRISPResso2.CRISPRessoReports import CRISPRessoReport

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(CRISPRessoShared.LogStreamHandler())

error = logger.critical
warn = logger.warning
debug = logger.debug
info = logger.info


def check_library(library_name):
        try:
                return __import__(library_name)
        except:
                error('You need to install %s module to use CRISPRessoCompare!' % library_name)
                sys.exit(1)


def parse_profile(profile_file):
    return np.loadtxt(profile_file, skiprows=1)


# EXCEPTIONS############################

class MixedRunningModeException(Exception):
    pass


class DifferentAmpliconLengthException(Exception):
    pass
############################


np = check_library('numpy')
pd = check_library('pandas')


_ROOT = os.path.abspath(os.path.dirname(__file__))


def normalize_name(name, output_folder_1, output_folder_2):
    get_name_from_folder = lambda x: os.path.basename(os.path.abspath(x)).replace('CRISPResso_on_', '')
    if not name:
        return '{0}_VS_{1}'.format(
            get_name_from_folder(output_folder_1),
            get_name_from_folder(output_folder_2),
        )
    else:
        return name


def get_matching_allele_files(run_info_1, run_info_2):
    def get_amplicon_info(run_info):
        return {
            amplicon['sequence']: {
                'name': amplicon_name,
                'guides': amplicon['sgRNA_orig_sequences'],
                'cut_points': amplicon['sgRNA_cut_points'],
                'allele_files': amplicon['allele_frequency_files'],
            }
            for amplicon_name, amplicon in run_info['results']['refs'].items()
        }
    amplicons_1 = get_amplicon_info(run_info_1)
    amplicons_2 = get_amplicon_info(run_info_2)
    matching_allele_files = []
    for sequence_1 in amplicons_1:
        if sequence_1 in amplicons_2:
            if amplicons_1[sequence_1]['cut_points'] != amplicons_2[sequence_1]['cut_points']:
                warn(f'Report 1 has different cut points than report 2 for amplicon {amplicons_1[sequence_1]["name"]}, skipping comparison')
                continue
            guides_1 = set(amplicons_1[sequence_1]['guides'])
            guides_2 = set(amplicons_2[sequence_1]['guides'])
            if not guides_1 & guides_2:
                warn(f'Report 1 has no shared guides with report 2 for amplicon {amplicons_1[sequence_1]["name"]}, skipping comparison')
                continue
            matching_allele_files.extend((f_1, f_2) for f_1, f_2 in zip(amplicons_1[sequence_1]['allele_files'], amplicons_2[sequence_1]['allele_files']))

    return matching_allele_files


def main():
    try:
        parser = CRISPRessoShared.getCRISPRessoArgParser("Compare", parser_title='CRISPRessoCompare Parameters')

        args = parser.parse_args()

        CRISPRessoShared.set_console_log_level(logger, args.verbosity, args.debug)

        description = ['~~~CRISPRessoCompare~~~', '-Comparison of two CRISPResso analyses-']
        compare_header = r'''
 ___________________________
| __ __      __      __  __ |
|/  /  \|\/||__) /\ |__)|_  |
|\__\__/|  ||   /--\| \ |__ |
|___________________________|
        '''
        compare_header = CRISPRessoShared.get_crispresso_header(description, compare_header)
        info(compare_header)

        if args.use_matplotlib or not CRISPRessoShared.is_C2Pro_installed():
            from CRISPResso2.plots import CRISPRessoPlot
        else:
            from CRISPRessoPro import plot as CRISPRessoPlot

        if args.zip_output and not args.place_report_in_output_folder:
            warn('Invalid argument combination: If zip_output is True then place_report_in_output_folder must also be True. Setting place_report_in_output_folder to True.')
            args.place_report_in_output_folder = True
        # check that the CRISPResso output is present and fill amplicon_info
        quantification_file_1, amplicon_names_1, amplicon_info_1 = CRISPRessoShared.check_output_folder(args.crispresso_output_folder_1)
        quantification_file_2, amplicon_names_2, amplicon_info_2 = CRISPRessoShared.check_output_folder(args.crispresso_output_folder_2)

        run_info_1 = CRISPRessoShared.load_crispresso_info(args.crispresso_output_folder_1)

        run_info_2 = CRISPRessoShared.load_crispresso_info(args.crispresso_output_folder_2)

        sample_1_name = args.sample_1_name
        if args.sample_1_name is None:
            sample_1_name = "Sample 1"
            if 'running_info' in run_info_1 and 'name' in run_info_1['running_info'] and run_info_1['running_info']['name']:
                sample_1_name = run_info_1['running_info']['name']

        sample_2_name = args.sample_2_name
        if args.sample_2_name is None:
            sample_2_name = "Sample 2"
            if 'running_info' in run_info_2 and 'name' in run_info_2['running_info'] and run_info_2['running_info']['name']:
                sample_2_name = run_info_2['running_info']['name']

        if sample_1_name == sample_2_name:
            sample_2_name += '_2'

        OUTPUT_DIRECTORY = 'CRISPRessoCompare_on_{0}'.format(normalize_name(
            args.name, args.crispresso_output_folder_1, args.crispresso_output_folder_2,
        ))

        if args.output_folder:
            OUTPUT_DIRECTORY = os.path.join(os.path.abspath(args.output_folder), OUTPUT_DIRECTORY)

        _jp = lambda filename: os.path.join(OUTPUT_DIRECTORY, filename)  # handy function to put a file in the output directory
        log_filename = _jp('CRISPRessoCompare_RUNNING_LOG.txt')

        try:
            info('Creating Folder %s' % OUTPUT_DIRECTORY, {'percent_complete': 0})
            os.makedirs(OUTPUT_DIRECTORY)
            info('Done!')
        except:
            warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename = _jp('CRISPRessoCompare_RUNNING_LOG.txt')
        logger.addHandler(logging.FileHandler(log_filename))
        logger.addHandler(CRISPRessoShared.StatusHandler(os.path.join(OUTPUT_DIRECTORY, 'CRISPRessoCompare_status.json')))

        with open(log_filename, 'w+') as outfile:
            outfile.write('[Command used]:\nCRISPRessoCompare %s\n\n[Execution log]:\n' % ' '.join(sys.argv))

        crispresso2Compare_info_file = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso2Compare_info.json')
        crispresso2_info = {'running_info': {}, 'results': {'alignment_stats': {}, 'general_plots': {}}}  # keep track of all information for this run to be pickled and saved at the end of the run
        crispresso2_info['running_info']['version'] = CRISPRessoShared.__version__
        crispresso2_info['running_info']['args'] = deepcopy(args)

        crispresso2_info['running_info']['log_filename'] = os.path.basename(log_filename)

        crispresso2_info['results']['general_plots']['summary_plot_names'] = []
        crispresso2_info['results']['general_plots']['summary_plot_titles'] = {}
        crispresso2_info['results']['general_plots']['summary_plot_labels'] = {}
        crispresso2_info['results']['general_plots']['summary_plot_datas'] = {}

        save_png = True
        if args.suppress_report:
            save_png = False

        # LOAD DATA
        amplicon_names_in_both = [amplicon_name for amplicon_name in amplicon_names_1 if amplicon_name in amplicon_names_2]

        # --- Build ComparePlotContext ---
        from CRISPResso2.plots.plot_context import ComparePlotContext

        plot_context = ComparePlotContext(
            args=args,
            run_data=crispresso2_info,
            output_directory=OUTPUT_DIRECTORY,
            save_png=save_png,
            _jp=_jp,
            custom_config={},
            amplicon_names=amplicon_names_in_both,
            sample_1_name=sample_1_name,
            sample_2_name=sample_2_name,
            run_info_1=run_info_1,
            run_info_2=run_info_2,
            amplicon_info_1=amplicon_info_1,
            amplicon_info_2=amplicon_info_2,
        )

        C2PRO_INSTALLED = CRISPRessoShared.is_C2Pro_installed()
        if C2PRO_INSTALLED:
            try:
                from CRISPRessoPro import hooks as pro_hooks
                pro_hooks.on_compare_plots_complete(plot_context, logger)
            except Exception as e:
                if args.halt_on_plot_fail:
                    raise
                logger.warning(f"CRISPRessoPro plugin hook failed: {e}")
        else:
            from CRISPResso2.plots.builtin_runners import run_builtin_compare_plots
            run_builtin_compare_plots(
                plot_context, crispresso2_info, CRISPRessoPlot, logger,
            )

        # Both paths publish Fisher/Bonferroni counts via run_data;
        # read them back for the summary-counts CSV below.
        sig_counts = plot_context.run_data.get('_sig_counts', {})
        sig_counts_quant_window = plot_context.run_data.get('_sig_counts_quant_window', {})

        debug('Calculating significant base counts...', {'percent_complete': 95})
        sig_counts_filename = _jp('CRISPRessoCompare_significant_base_counts.txt')
        with open(sig_counts_filename, 'w') as fout:
            fout.write('Amplicon\tModification\tsig_base_count\tsig_base_count_quant_window\n')
            for amplicon_name in amplicon_names_in_both:
                for mod in ['Insertions', 'Deletions', 'Substitutions', 'All_modifications']:
                    val = np.nan
                    if amplicon_name in sig_counts and mod in sig_counts[amplicon_name]:
                        val = sig_counts[amplicon_name][mod]

                    val_quant_window = np.nan
                    if amplicon_name in sig_counts_quant_window and mod in sig_counts_quant_window[amplicon_name]:
                        val_quant_window = sig_counts_quant_window[amplicon_name][mod]
                    line = "%s\t%s\t%s\t%s\n" % (amplicon_name, mod, val, val_quant_window)
                    fout.write(line)
        crispresso2_info['running_info']['sig_counts_report_location'] = sig_counts_filename

        if not args.suppress_report:
            if (args.place_report_in_output_folder):
                report_name = _jp("CRISPResso2Compare_report.html")
            else:
                report_name = OUTPUT_DIRECTORY + '.html'
            if C2PRO_INSTALLED:
                from CRISPRessoPro import hooks as pro_hooks
                pro_hooks.make_compare_report(crispresso2_info, report_name, OUTPUT_DIRECTORY, _ROOT, logger, plot_context)
            else:
                CRISPRessoReport.make_compare_report_from_folder(report_name, crispresso2_info, OUTPUT_DIRECTORY, _ROOT, logger)
            crispresso2_info['running_info']['report_location'] = report_name
            crispresso2_info['running_info']['report_filename'] = os.path.basename(report_name)

        CRISPRessoShared.write_crispresso_info(crispresso2Compare_info_file, crispresso2_info)

        if args.zip_output:
            CRISPRessoShared.zip_results(OUTPUT_DIRECTORY)

        info('Analysis Complete!', {'percent_complete': 100})
        info(CRISPRessoShared.get_crispresso_footer())
        sys.exit(0)

    except Exception as e:
        debug_flag = False
        if 'args' in vars() and 'debug' in args:
            debug_flag = args.debug

        if debug_flag:
            traceback.print_exc(file=sys.stdout)

        error('\n\nERROR: %s' % e)
        sys.exit(-1)
