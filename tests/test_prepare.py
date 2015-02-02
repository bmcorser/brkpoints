from nose.tools import *
from pkg_resources import resource_filename
from brkpoints.prepare import main, start_position, split_transcript


import os
import yaml
import tempfile


OUTER_TODAY = None


def setup():
    'Setup function that runs before every test'
    print "SETUP!"


def teardown():
    'Teardown function that runs after every test'
    print "TEAR DOWN!"


def test_start_position():
    'Test start_position'
    test_value = 99
    component = {'start_position': test_value}
    assert_equal(start_position(component), test_value)


def compare_transcript_details(obtained, expected):
    """ Compare transcript objects """
    for key in expected:
        assert_true(key in obtained)
        if key != 'components':
            assert_equal(expected[key], obtained[key])
    assert_true('components' in expected)
    assert_equal(len(obtained['components']),
                 len(expected['components']))
    for i in range(0, 2):
        for key in expected['components'][i]:
            assert_equal(obtained['components'][i][key],
                         expected['components'][i][key])


def test_split_transcript():
    'Test split_transcript'
    # Work area
    workdir = tempfile.mkdtemp()
    # Open handle to sequence file
    sequences_file = workdir + '/sequences.fa'
    # Create inputs
    fasta_seq = {
        'header': 'ENST01|ABL|X|1|10|1|ENSE01;ENSE02|2;8|5;9',
        'body': 'GTATCCGCTG'
        }
    output_fh = open(sequences_file, 'w')
    transcripts = {}
    # Run the splitter
    seq_counter = split_transcript(
        fasta_seq, output_fh, transcripts, 0)
    output_fh.close()
    # Comparison of output with expected
    expected = {
        'name': 'ENST01',
        'components': [
            {'name': 'ENSE01',
             'start_position': 2, 'end_position': 5,
             'type': 'exon'},
            {'name': 'ENST01:1',
             'start_position': 6, 'end_position': 7,
             'type': 'intron'},
            {'name': 'ENSE02',
             'start_position': 8, 'end_position': 9,
             'type': 'exon'},
            ],
        'chromosome_name': 'X',
        'start_position': 1,
        'end_position': 10,
        'strand': 1
        }
    assert_equal(seq_counter, 3)
    assert_true('ABL' in transcripts)
    assert_equal(len(transcripts['ABL']), 1)
    compare_transcript_details(transcripts['ABL'][0], expected)
    input_fh = open(sequences_file)
    lines = input_fh.readlines()
    input_fh.close()
    assert_equal(lines[1].rstrip(), 'TATC')
    assert_equal(lines[3].rstrip(), 'CG')
    assert_equal(lines[5].rstrip(), 'CT')
    # Reverse the transcript
    fasta_seq = {
        'header': 'ENST01|ABL|X|1|10|-1|ENSE01;ENSE02|2;8|5;9',
        'body': 'CAGCGGATAC'
        }
    output_fh = open(sequences_file, 'w')
    transcripts = {}
    # Rerun
    seq_counter = split_transcript(
        fasta_seq, output_fh, transcripts, 0)
    output_fh.close()
    # Comparison of output with expected
    expected['strand'] = -1
    assert_equal(seq_counter, 3)
    assert_true('ABL' in transcripts)
    assert_equal(len(transcripts['ABL']), 1)
    compare_transcript_details(transcripts['ABL'][0], expected)
    input_fh = open(sequences_file)
    lines = input_fh.readlines()
    input_fh.close()
    assert_equal(lines[1].rstrip(), 'GATA')
    assert_equal(lines[3].rstrip(), 'CG')
    assert_equal(lines[5].rstrip(), 'AG')


def test_main():
    'Test main'
    # Work area
    workdir = tempfile.mkdtemp()
    # Input file
    transcripts_file = \
        resource_filename(__name__, 'data/ensembl_GRCh38_reduced.fa')
    # Output files
    sequences_file = workdir + '/sequences.fa'
    details_file = workdir + '/details.yaml'
    # Run main
    main([transcripts_file,
          sequences_file,
          details_file])
    # Load the resulting transcript details_file (YAML)
    stream = open(details_file, 'r')
    transcript_details = yaml.load(stream)
    # Comparison of output with expected
    # Input file has four HGNC gene ids and 2 transcripts per gene
    # So expecting 8 transcripts
    seq_counter = 0
    base_counter = 0
    assert_equal(len(transcript_details), 4)
    for hgnc_symbol in transcript_details:
        transcripts = transcript_details[hgnc_symbol]
        assert_equal(len(transcripts), 2)
        for transcript in transcripts:
            assert_true('components' in transcript)
            seq_counter += len(transcript['components'])
            for component in transcript['components']:
                base_counter += \
                    component['end_position'] - \
                    component['start_position'] + 1
    # Compare with the number of sequences written and check
    # that the total DNA length in and out is the same
    input_fh = open(sequences_file)
    lines = input_fh.readlines()
    input_fh.close()
    cmp_seq_counter = 0
    cmp_base_counter = 0
    for line in lines:
        if line.startswith('>'):
            cmp_seq_counter += 1
        else:
            stripped = line.rstrip()
            cmp_base_counter += len(stripped)
    assert_equal(seq_counter, cmp_seq_counter)
    assert_equal(base_counter, cmp_base_counter)
