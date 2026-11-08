"""HDR windows and cut points cloned from a guide must land on the same bases as reference 1's."""
import json
import subprocess
import sys

AMPLICON = (
    'CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGC'
    'TTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG'
)
# 4 bp longer at the 5' end and 11 bp shorter at the deletion, so HDR position = reference position - 7 past the deletion.
HDR = (
    'CGGCCGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCTGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCTTTTCCGAGC'
    'TTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG'
)
GUIDE_PAST_THE_DELETION = 'GTCTCAAGCACTACCTACGT'
CONTEXT = 8


def _window_context(sequence, include_idxs):
    return sequence[include_idxs[0] - CONTEXT:include_idxs[-1] + CONTEXT + 1]


def test_guide_window_clones_onto_the_same_hdr_bases(tmp_path):
    reads = tmp_path / 'reads.fastq'
    with reads.open('w') as fh:
        for i, sequence in enumerate([AMPLICON, HDR] * 50):
            fh.write(f'@r{i}\n{sequence}\n+\n{"I" * len(sequence)}\n')
    subprocess.run(
        [
            sys.executable, '-c', 'from CRISPResso2.CRISPRessoCORE import main; main()',
            '-r1', str(reads), '-a', AMPLICON, '-e', HDR, '-g', GUIDE_PAST_THE_DELETION,
            '--suppress_plots', '--suppress_report', '-o', str(tmp_path), '-n', 'clone',
        ],
        check=True,
        capture_output=True,
    )
    info = json.loads((tmp_path / 'CRISPResso_on_clone' / 'CRISPResso2_info.json').read_text())
    refs = info['results']['refs']
    reference_window = refs['Reference']['include_idxs']['value']
    hdr_window = refs['HDR']['include_idxs']['value']
    assert _window_context(HDR, hdr_window) == _window_context(AMPLICON, reference_window)
    reference_cut = refs['Reference']['sgRNA_cut_points'][0]
    hdr_cut = refs['HDR']['sgRNA_cut_points'][0]
    assert HDR[hdr_cut - CONTEXT:hdr_cut + CONTEXT + 1] == AMPLICON[reference_cut - CONTEXT:reference_cut + CONTEXT + 1]
