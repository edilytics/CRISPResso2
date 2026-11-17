"""Cloned coordinates across indels, truncations and reverse-complement guides, with plots off and on."""
import importlib.util
import json
import subprocess
import sys

import pytest

AMPLICON = (
    'CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGC'
    'TTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG'
)
GUIDE = 'GTCTCAAGCACTACCTACGT'
GUIDE_START = AMPLICON.index(GUIDE)
INSERTION = 'TGACTGATCA'
GUIDE_NEAR_5P = AMPLICON[3:23]
GUIDE_NEAR_3P = AMPLICON[-40:-20]
WITHOUT_GUIDE = AMPLICON[:GUIDE_START - 3] + AMPLICON[GUIDE_START + 23:]
PRO_INSTALLED = importlib.util.find_spec('CRISPRessoPro') is not None
# CRISPRessoPro's amino acid table reads the first guide, so any -c run without -g crashes there
pro_needs_a_guide = pytest.mark.xfail(PRO_INSTALLED, reason='CRISPRessoPro amino_acid_table_data_func indexes sgRNA_orig_sequences[0]')
EXON = AMPLICON[GUIDE_START - 30:GUIDE_START + 60]
EXON_START = GUIDE_START - 30


def _revcomp(sequence):
    return sequence[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def _edit(start, deleted=0, inserted=''):
    return AMPLICON[:start] + inserted + AMPLICON[start + deleted:]


def _crispresso(tmp_path, target, *args, multi_amplicon=False, plots=False, name='clone'):
    reads = tmp_path / f'{name}.fastq'
    with reads.open('w') as fh:
        for i, sequence in enumerate([AMPLICON, target] * 15):
            fh.write(f'@r{i}\n{sequence}\n+\n{"I" * len(sequence)}\n')
    amplicons = ['-a', f'{AMPLICON},{target}', '-an', 'Reference,HDR'] if multi_amplicon else ['-a', AMPLICON, '-e', target]
    output = ['--halt_on_plot_fail'] if plots else ['--suppress_plots', '--suppress_report']
    return subprocess.run(
        [
            sys.executable, '-c', 'from CRISPResso2.CRISPRessoCORE import main; main()',
            '-r1', str(reads), *amplicons, *output,
            '-o', str(tmp_path), '-n', name, *args,
        ],
        check=False,
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )


def _assert_ok(result):
    assert result.returncode == 0, f'exit {result.returncode}: {(result.stdout + result.stderr)[-400:]}'


def _refs(tmp_path, name='clone'):
    info = json.loads((tmp_path / f'CRISPResso_on_{name}' / 'CRISPResso2_info.json').read_text())
    return info['results']['refs']


def _run(tmp_path, target, *args, **kwargs):
    _assert_ok(_crispresso(tmp_path, target, *args, **kwargs))
    return _refs(tmp_path, kwargs.get('name', 'clone'))


def _window(ref):
    return ref['include_idxs']['value']


def _frameshift_reads(tmp_path):
    text = (tmp_path / 'CRISPResso_on_clone' / 'HDR.Frameshift_analysis.txt').read_text()
    return int(text.split('Frameshift mutation:')[1].split(' reads')[0])


def test_insertion_upstream_moves_coordinates_right(tmp_path):
    refs = _run(tmp_path, _edit(60, inserted=INSERTION), '-g', GUIDE)
    assert refs['HDR']['sgRNA_cut_points'] == [refs['Reference']['sgRNA_cut_points'][0] + 10]
    assert refs['HDR']['sgRNA_intervals'] == [[GUIDE_START + 10, GUIDE_START + 29]]
    assert _window(refs['HDR']) == [x + 10 for x in _window(refs['Reference'])]


def test_reverse_complement_guide_past_deletion(tmp_path):
    refs = _run(tmp_path, _edit(60, deleted=11), '-g', _revcomp(GUIDE))
    assert refs['HDR']['sgRNA_cut_points'] == [refs['Reference']['sgRNA_cut_points'][0] - 11]
    assert _window(refs['HDR']) == [x - 11 for x in _window(refs['Reference'])]


def test_deletion_interior_to_guide_keeps_cut_on_its_base(tmp_path):
    target = _edit(GUIDE_START + 6, deleted=3)
    refs = _run(tmp_path, target, '-g', GUIDE)
    reference_cut = refs['Reference']['sgRNA_cut_points'][0]
    hdr_cut = refs['HDR']['sgRNA_cut_points'][0]
    assert hdr_cut == reference_cut - 3
    assert target[hdr_cut - 5:hdr_cut + 6] == AMPLICON[reference_cut - 5:reference_cut + 6]
    assert refs['HDR']['sgRNA_intervals'] == [[GUIDE_START, GUIDE_START + 16]]


def test_cloned_extra_amplicon_without_the_guide(tmp_path):
    target = _edit(GUIDE_START + 5, inserted=INSERTION)
    refs = _run(tmp_path, target, '-g', GUIDE, multi_amplicon=True)
    assert refs['HDR']['idx_cloned_from'] == 'Reference'
    assert refs['HDR']['sgRNA_cut_points'] == [refs['Reference']['sgRNA_cut_points'][0] + 10]
    assert refs['HDR']['sgRNA_intervals'] == [[GUIDE_START, GUIDE_START + 29]]


def test_exon_length_modification_has_the_edit_sign(tmp_path):
    deleted = _run(tmp_path, _edit(GUIDE_START + 40, deleted=6), '-g', GUIDE, '-c', EXON)
    assert deleted['HDR']['exon_intervals'] == [[EXON_START, EXON_START + 84]]
    assert deleted['HDR']['exon_len_mods'] == [-6]


def test_exon_insertion_length_modification(tmp_path):
    refs = _run(tmp_path, _edit(GUIDE_START + 40, inserted='ACT'), '-g', GUIDE, '-c', EXON)
    assert refs['HDR']['exon_intervals'] == [[EXON_START, EXON_START + 93]]
    assert refs['HDR']['exon_len_mods'] == [3]


def test_qwc_clones_onto_the_same_bases(tmp_path):
    refs = _run(tmp_path, _edit(0, inserted='CGGC')[:99] + _edit(0, inserted='CGGC')[110:], '-g', GUIDE, '-qwc', '120-140')
    assert _window(refs['HDR']) == [x - 7 for x in _window(refs['Reference'])]


def test_qwc_without_guide_on_a_shorter_reference(tmp_path):
    refs = _run(tmp_path, _edit(60, deleted=11), '-qwc', '120-140')
    assert _window(refs['HDR']) == [x - 11 for x in _window(refs['Reference'])]


def test_qwc_clips_a_partly_deleted_window(tmp_path):
    refs = _run(tmp_path, _edit(125, deleted=6), '-qwc', '120-140')
    assert _window(refs['HDR']) == list(range(120, 135))


def test_guide_overlapping_bases_missing_from_the_target_3p_end(tmp_path):
    target = AMPLICON[:-30]
    refs = _run(tmp_path, target, '-g', AMPLICON[-40:-20])
    assert refs['HDR']['sgRNA_intervals'] == [[len(AMPLICON) - 40, len(target) - 1]]


@pytest.mark.parametrize('multi_amplicon', [False, True], ids=['hdr', 'extra_amplicon'])
def test_guide_overlapping_bases_missing_from_the_target_5p_end(tmp_path, multi_amplicon):
    refs = _run(tmp_path, AMPLICON[25:], '-g', AMPLICON[20:40], multi_amplicon=multi_amplicon)
    assert refs['HDR']['sgRNA_intervals'] == [[0, 14]]


def test_exon_overlapping_bases_missing_from_the_target_5p_end(tmp_path):
    refs = _run(tmp_path, AMPLICON[30:], '-g', GUIDE, '-c', AMPLICON[20:80])
    assert refs['HDR']['exon_intervals'] == [[0, 50]]
    assert refs['HDR']['exon_len_mods'] == [0]


@pro_needs_a_guide
def test_truncated_exon_is_not_a_frameshift(tmp_path):
    refs = _run(tmp_path, AMPLICON[:-40], '-c', AMPLICON[-70:-10], plots=True)
    assert refs['HDR']['exon_len_mods'] == [0]
    assert _frameshift_reads(tmp_path) == 0


def test_window_never_contains_a_negative_index(tmp_path):
    refs = _run(
        tmp_path, AMPLICON[25:], '-g', AMPLICON[28:48], '-wc', '-10', '-w', '20',
        '--exclude_bp_from_left', '0', '--exclude_bp_from_right', '0',
    )
    assert min(_window(refs['HDR'])) >= 0


def test_exon_start_inside_deletion_counts_only_exon_bases(tmp_path):
    refs = _run(tmp_path, _edit(EXON_START - 1, deleted=4), '-g', GUIDE, '-c', EXON)
    assert refs['HDR']['exon_len_mods'] == [-3]


def test_fully_deleted_exon_has_no_positions(tmp_path):
    refs = _run(tmp_path, _edit(GUIDE_START + 60, deleted=20), '-g', GUIDE, '-c', AMPLICON[GUIDE_START + 62:GUIDE_START + 78], plots=True)
    assert refs['HDR']['exon_intervals'] == [[GUIDE_START + 60, GUIDE_START + 60]]
    assert refs['HDR']['exon_positions'] == []
    assert refs['HDR']['exon_len_mods'] == [-16]


def test_exons_clone_without_a_guide(tmp_path):
    refs = _run(tmp_path, _edit(GUIDE_START + 10, deleted=6), '-c', EXON)
    assert refs['HDR']['exon_intervals'] == [[EXON_START, EXON_START + 84]]


@pytest.mark.parametrize('multi_amplicon', [False, True], ids=['hdr', 'extra_amplicon'])
def test_cut_point_in_bases_missing_from_the_target_5p_end(tmp_path, multi_amplicon):
    refs = _run(tmp_path, AMPLICON[30:], '-g', GUIDE_NEAR_5P, multi_amplicon=multi_amplicon, plots=True)
    assert refs['HDR']['sgRNA_cut_points'] == [0]


def test_cut_point_in_missing_5p_bases_without_plots(tmp_path):
    refs = _run(tmp_path, AMPLICON[30:], '-g', GUIDE_NEAR_5P)
    assert refs['HDR']['sgRNA_cut_points'] == [0]


@pytest.mark.parametrize('guide_args', [pytest.param([], marks=pro_needs_a_guide), ['-g', GUIDE]], ids=['no_guide', 'guide'])
def test_exon_in_bases_missing_from_the_target_3p_end(tmp_path, guide_args):
    target = AMPLICON[:-40]
    refs = _run(tmp_path, target, '-c', AMPLICON[-30:-10], *guide_args, plots=True)
    assert refs['HDR']['exon_intervals'] == [[len(target), len(target)]]
    assert refs['HDR']['exon_len_mods'] == [0]
    assert _frameshift_reads(tmp_path) == 0


def test_compare_accepts_an_empty_cloned_window(tmp_path):
    for name in ('one', 'two'):
        refs = _run(tmp_path, WITHOUT_GUIDE, '-g', GUIDE, name=name, plots=True)
    assert _window(refs['HDR']) == []
    result = subprocess.run(
        [
            sys.executable, '-c', 'from CRISPResso2.CRISPRessoCompareCORE import main; main()',
            str(tmp_path / 'CRISPResso_on_one'), str(tmp_path / 'CRISPResso_on_two'),
            '-o', str(tmp_path), '-n', 'compare',
        ],
        check=False,
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    _assert_ok(result)


def test_empty_cloned_window_warns(tmp_path):
    result = _crispresso(tmp_path, WITHOUT_GUIDE, '-g', GUIDE)
    _assert_ok(result)
    assert "covers no bases of 'HDR'" in result.stdout + result.stderr


def test_guide_missing_from_the_middle_of_the_target_plots(tmp_path):
    refs = _run(tmp_path, WITHOUT_GUIDE, '-g', GUIDE, plots=True)
    assert refs['HDR']['sgRNA_intervals'] == [[GUIDE_START - 3, GUIDE_START - 4]]
    assert refs['HDR']['sgRNA_mismatches'] == [[]]


def test_guide_missing_from_the_target_plots_base_editor_output(tmp_path):
    _run(tmp_path, WITHOUT_GUIDE, '-g', GUIDE, '--base_editor_output', '-wc', '-10', '-w', '10', plots=True)


def test_guide_start_at_the_target_length_plots(tmp_path):
    target = AMPLICON[:-45]
    refs = _run(tmp_path, target, '-g', GUIDE_NEAR_3P, plots=True)
    assert refs['HDR']['sgRNA_intervals'][0][0] == len(target)


def test_reverse_complement_guide_missing_from_the_target_3p_end(tmp_path):
    _run(tmp_path, AMPLICON[:-45], '-g', _revcomp(GUIDE_NEAR_3P), plots=True)


def test_plot_indexes_follow_the_cut_point(tmp_path):
    refs = _run(tmp_path, _edit(60, inserted=INSERTION), '-g', GUIDE, plots=True)
    cut = refs['HDR']['sgRNA_cut_points'][0]
    assert refs['HDR']['sgRNA_plot_idxs'][0] == list(range(cut - 19, cut + 21))


@pro_needs_a_guide
def test_exon_without_a_guide_gets_frameshift_analysis_and_no_cut_points(tmp_path):
    refs = _run(tmp_path, _edit(GUIDE_START + 10, deleted=5), '-c', EXON, plots=True)
    assert refs['HDR']['contains_coding_seq'] is True
    assert refs['HDR']['exon_len_mods'] == [-5]
    assert refs['HDR']['sgRNA_cut_points'] == []
    assert refs['HDR']['contains_guide'] is False
    assert _frameshift_reads(tmp_path) == 15


def test_quantification_window_wholly_deleted_from_the_target_plots(tmp_path):
    refs = _run(tmp_path, AMPLICON[:100] + AMPLICON[140:], '-qwc', '110-120', plots=True)
    assert _window(refs['HDR']) == []


def test_exon_at_the_start_of_a_target_with_a_5p_extension(tmp_path):
    refs = _run(tmp_path, 'CGGC' + AMPLICON, '-g', GUIDE, '-c', AMPLICON[:60])
    assert refs['HDR']['exon_intervals'] == [[4, 64]]
    assert refs['HDR']['exon_len_mods'] == [0]


def test_exon_only_template_clones_onto_an_extra_amplicon(tmp_path):
    refs = _run(tmp_path, _edit(GUIDE_START + 10, deleted=6), '-c', EXON, multi_amplicon=True)
    assert refs['HDR']['idx_cloned_from'] == 'Reference'
    assert refs['HDR']['exon_intervals'] == [[EXON_START, EXON_START + 84]]
    assert refs['HDR']['exon_len_mods'] == [-6]


def test_cut_point_past_the_target_3p_end_lands_on_its_last_base(tmp_path):
    target = AMPLICON[:-10]
    refs = _run(tmp_path, target, '-g', AMPLICON[-25:-5], '--exclude_bp_from_right', '0', plots=True)
    assert refs['HDR']['sgRNA_cut_points'] == [len(target) - 1]


def test_insertion_inside_the_window_widens_it(tmp_path):
    reference_cut = GUIDE_START + 16
    refs = _run(tmp_path, _edit(reference_cut + 3, inserted='TTT'), '-g', GUIDE, '-w', '5')
    reference_window = _window(refs['Reference'])
    assert _window(refs['HDR']) == list(range(reference_window[0], reference_window[-1] + 4))
