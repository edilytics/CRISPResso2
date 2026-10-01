"""Smoke tests for standalone scripts using the legacy allele-plot API."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace
import sys
import zipfile

import pandas as pd
import pytest

from CRISPResso2 import CRISPRessoShared


@pytest.mark.parametrize('script, extra_args, label', [
    ('plotAmbiguous', [], 'guide'),
    ('plotCustomAllelePlot', ['--plot_left', '2', '--plot_right', '2'], 'guide'),
    ('plotCustomAllelePlot', ['--plot_left', '2', '--plot_right', '2', '--plot_center', '9'], 'custom'),
])
def test_allele_plot_script_main(tmp_path, monkeypatch, script, extra_args, label):
    path = Path(__file__).resolve().parents[2] / 'scripts' / (script + '.py')
    spec = importlib.util.spec_from_file_location(script, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    reference = 'ACGT' * 5
    df = pd.DataFrame([{
        'Aligned_Sequence': reference,
        'Reference_Sequence': reference,
        'Reference_Name': 'AMBIGUOUS_r' if script == 'plotAmbiguous' else 'r',
        'Read_Status': 'UNMODIFIED',
        'ref_positions': list(range(len(reference))),
        'n_deleted': 0, 'n_inserted': 0, 'n_mutated': 0,
        '#Reads': 100, '%Reads': 100.0,
    }])
    with zipfile.ZipFile(tmp_path / 'alleles.zip', 'w') as archive:
        archive.writestr('alleles.txt', df.to_csv(sep='\t', index=False))

    info = {
        'running_info': {
            'args': SimpleNamespace(
                write_detailed_allele_table=True, plot_window_size=2,
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                annotate_wildtype_allele='****',
            ),
            'allele_frequency_table_zip_filename': 'alleles.zip',
            'allele_frequency_table_filename': 'alleles.txt',
        },
        'results': {
            'ref_names': ['r'],
            'refs': {'r': {
                'sequence': reference, 'sgRNA_sequences': ['ACGT'],
                'sgRNA_cut_points': [9], 'sgRNA_plot_cut_points': [True],
                'sgRNA_intervals': [(8, 11)], 'sgRNA_names': ['guide'],
                'sgRNA_mismatches': [[]], 'sgRNA_plot_idxs': [[8, 9, 10, 11]],
            }},
        },
    }
    monkeypatch.setattr(CRISPRessoShared, 'load_crispresso_info', lambda _: info)
    monkeypatch.setattr(sys, 'argv', [
        str(path), '-f', str(tmp_path), '-o', str(tmp_path / 'plot'),
        '--use_matplotlib', '--max_rows', '1', '--save_png', *extra_args,
    ])
    # Exercise main's backend selection, the legacy six-result prep API, and
    # the real heatmap signature without requiring marker metadata.
    module.main()
    for suffix in ('pdf', 'png'):
        assert (tmp_path / f'plot_r_{label}.{suffix}').stat().st_size > 0
