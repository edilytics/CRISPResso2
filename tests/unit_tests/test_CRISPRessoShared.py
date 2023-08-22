import pytest
import argparse
from unittest.mock import Mock
from CRISPResso2 import CRISPResso2Align, CRISPRessoShared

ALN_MATRIX = CRISPResso2Align.read_matrix('./CRISPResso2/EDNAFULL')


def test_get_mismatches():
    mismatch_cords = CRISPRessoShared.get_mismatches(
        'ATTA',
        'ATTA',
        ALN_MATRIX,
        -5,
        -3,
    )
    assert len(mismatch_cords) == 0

    mismatch_cords = CRISPRessoShared.get_mismatches(
        'GCAGTGGGCGCGCTA',
        'CCCACTGAAGGCCC',
        ALN_MATRIX,
        -5,
        -3,
    )
    assert len(mismatch_cords) == 6


def test_check_custom_config_base(monkeypatch):
    default_config =  {
        "colors": 
        {
            'Substitution': '#0000FF',
            'Insertion': '#008000',
            'Deletion': '#FF0000',
            'A': '#7FC97F',
            'T': '#BEAED4',
            'C': '#FDC086',
            'G': '#FFFF99',
            'N': '#C8C8C8',
            '-': '#C1C1C1'
        }, 
        "guardrails": 
        {
            'min_total_reads': 10000,
            'alignedCutoff': 0.9,
            'alternateAlignment': 0.3,
            'minRatioOfModsInToOut': 0.01,
            'modificationsAtEnds': 0.01,
            'outsideWindowMaxSubRate': 0.002,
            'maxRateOfSubs': 0.3,
            'guide_len': 19,
            'amplicon_len': 50,
            'ampliconToReadLen': 1.5
        }       
    }

    args = argparse.Namespace(config_file='./tests/constTestInputs/Oceano.json')
    mock_parser = Mock()
    mock_parser.parse_args = Mock(return_value=args)
    monkeypatch.setattr(argparse, 'ArgumentParser', Mock(return_value=mock_parser))

    config = CRISPRessoShared.check_custom_config(args, False)
    assert config == default_config


def test_check_custom_config_pro(monkeypatch):
    custom_config =  {
        "colors": {
            "Substitution": "#0089b3", 
            "Insertion": "#ade8f3", 
            "Deletion": "#0197c5",
            "A": "#0179b1", 
            "C": "#00b3d8", 
            "G": "#005d8f", 
            "T": "#90e0ef", 
            "N": "#03045e", 
            "-": "#03045e"   
        }, 
        "guardrails": 
        {
            'min_total_reads': 10000,
            'alignedCutoff': 0.9,
            'alternateAlignment': 0.3,
            'minRatioOfModsInToOut': 0.01,
            'modificationsAtEnds': 0.01,
            'outsideWindowMaxSubRate': 0.002,
            'maxRateOfSubs': 0.3,
            'guide_len': 19,
            'amplicon_len': 50,
            'ampliconToReadLen': 1.5
        }       
    }
    
    args = argparse.Namespace(config_file='./tests/constTestInputs/Oceano.json')
    mock_parser = Mock()
    mock_parser.parse_args = Mock(return_value=args)
    monkeypatch.setattr(argparse, 'ArgumentParser', Mock(return_value=mock_parser))
    
    config = CRISPRessoShared.check_custom_config(args, True)
    assert config == custom_config


def test_check_custom_config_no_file(monkeypatch):
    default_config =  {
        "colors": 
        {
            'Substitution': '#0000FF',
            'Insertion': '#008000',
            'Deletion': '#FF0000',
            'A': '#7FC97F',
            'T': '#BEAED4',
            'C': '#FDC086',
            'G': '#FFFF99',
            'N': '#C8C8C8',
            '-': '#C1C1C1'
        }, 
        "guardrails": 
        {
            'min_total_reads': 10000,
            'alignedCutoff': 0.9,
            'alternateAlignment': 0.3,
            'minRatioOfModsInToOut': 0.01,
            'modificationsAtEnds': 0.01,
            'outsideWindowMaxSubRate': 0.002,
            'maxRateOfSubs': 0.3,
            'guide_len': 19,
            'amplicon_len': 50,
            'ampliconToReadLen': 1.5
        }       
    }

    args = argparse.Namespace(config_file=None)
    mock_parser = Mock()
    mock_parser.parse_args = Mock(return_value=args)
    monkeypatch.setattr(argparse, 'ArgumentParser', Mock(return_value=mock_parser))

    config = CRISPRessoShared.check_custom_config(args, True)
    assert config == default_config


def test_safety_check():
    CRISPRESSO2_INFO = {
        'version': '2.2.14', 
        'results': {
            'alignment_stats': {
                'counts_total': {
                    'FANC': 192, 
                    'HDR': 16
                    }, 
                },  
            'ref_names': ['FANC', 'HDR'], 
            'refs': {
                'FANC': {
                    'name': 'FANC',  
                    'sequence_length': 223, 
                    'sgRNA_sequences': ['GGAATCCCTTCTGCAGCACC', 'GGCCTTGCAGTGGGCGCGCTA', 'CCCACTGCAAGGCCC'], 
                },
                'HDR': {
                    'name': 'HDR', 
                    'sequence_length': 216, 
                    'sgRNA_sequences': ['GGAATCCCTTCTGCAGCACC', 'GGCCTTGCAGTGGGCGCGCTA', 'CCCACTGCAAGGCCC'],
                }
            }
        }
    }
    ALN_STATS = {
        'N_TOT_READS': 231, 
        'N_CACHED_ALN': 37, 
        'N_CACHED_NOTALN': 0, 
        'N_COMPUTED_ALN': 171, 
        'N_COMPUTED_NOTALN': 23, 
        'N_GLOBAL_SUBS': 113, 
        'N_SUBS_OUTSIDE_WINDOW': 97, 
        'N_MODS_IN_WINDOW': 16, 
        'N_MODS_OUTSIDE_WINDOW': 143.0, 
        'N_READS_IRREGULAR_ENDS': 208, 
        'READ_LENGTH': 250
    }
    guardrails = {
        'min_total_reads': 10000,
        'alignedCutoff': 0.9,
        'alternateAlignment': 0.3,
        'minRatioOfModsInToOut': 0.01,
        'modificationsAtEnds': 0.01,
        'outsideWindowMaxSubRate': 0.002,
        'maxRateOfSubs': 0.3,
        'guide_len': 19,
        'amplicon_len': 50,
        'ampliconToReadLen': 1.5
    }  
    expected_messages = [
        '<div class="alert alert-danger"><strong>Guardrail Warning!</strong> Low number of total reads: <10000</div>', 
        '<div class="alert alert-danger"><strong>Guardrail Warning!</strong> <=90.0% of expected reads were aligned to amplicon: HDR</div>', 
        '<div class="alert alert-danger"><strong>Guardrail Warning!</strong> >=30.0% more reads than expected were aligned to amplicon: FANC</div>', 
        '<div class="alert alert-danger"><strong>Guardrail Warning!</strong> >=1.0% of reads have modifications at the start or end. </div>', 
        '<div class="alert alert-danger"><strong>Guardrail Warning!</strong> >=0.2% of substitutions were outside of the quantification window. </div>', 
        '<div class="alert alert-danger"><strong>Guardrail Warning!</strong> >=30.0% of modifications were substitutions. This could potentially indicate poor sequencing quality. </div>', 
        '<div class="alert alert-danger"><strong>Guardrail Warning!</strong> guide length <19: CCCACTGCAAGGCCC</div>']
    messages = CRISPRessoShared.safety_check(CRISPRESSO2_INFO, ALN_STATS, guardrails=guardrails)
    assert messages == expected_messages

