from nose.tools import *
from pkg_resources import resource_filename
from brkpoints.find import main, parse_translocations_file


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


def compare_translocations(obtained, expected):
    """ Compare translocation objects, with or without alignment """
    for key in expected:
        assert_true(key in obtained)
        if key != 'partners':
            assert_equal(expected[key], obtained[key])
    assert_true('partners' in obtained)
    assert_equal(len(obtained['partners']), 2)
    expected_partners = expected['partners']
    obtained_partners = obtained['partners']
    for i in range(0, 1):
        for key in expected_partners[i]:
            assert_true(key in obtained_partners[i])
            if key != 'aln':
                assert_equal(expected_partners[i][key],
                             obtained_partners[i][key])
    if 'aln' in expected_partners:
        for i in range(0, 1):
            expected_aln = expected_partners[i]['partners']['aln']
            obtained_aln = obtained_partners[i]['partners']['aln']
            # Check chromosome and strand
            for key in ['chromosome_name', 'strand']:
                assert_true(key in obtained_aln)
                assert_equal(expected_aln[key],
                             obtained_aln[key])


def test_parse_translocations_file():
    'Test find.parse_translocations_file'
    # Input files
    translocations_file = resource_filename(
        __name__,
        'data/allseqs_TICdb_3.3_reduced.txt')
    # Run it
    translocations = parse_translocations_file(translocations_file)
    # Compare with expected output
    assert_equal(len(translocations), 17)
    # First translocation
    expected = {
        'xref': '18594527',
        'partners': [
            {'disp': '5prime', 'hgnc_symbol': 'ACSL3'},
            {'disp': '3prime', 'hgnc_symbol': 'ETV1'},
            ],
        'junction_sequence':
            ('ctgtgtcacaccaccttagcctcttgatcgaggaagTGCCTATGATCA' +
             'GAAGCCACAAGTGGGAATGAGGC')
        }
    # Compare
    compare_translocations(translocations[0], expected)


def test_main():
    'Test find.main'
    # Work area
    workdir = tempfile.mkdtemp()
    # Input files
    translocations_file = resource_filename(
        __name__,
        'data/allseqs_TICdb_3.3_reduced.txt')
    sequences_file = resource_filename(
        __name__,
        'data/sequences.fa')
    details_file = resource_filename(
        __name__,
        'data/details.yaml')
    blast_bin_dir = '/usr/bin'
    # Output files
    output_file = workdir + '/find.yaml'
    # Run main
    main([translocations_file,
          sequences_file,
          details_file,
          workdir,
          blast_bin_dir,
          output_file])
    # Load the output file (YAML)
    stream = open(output_file, 'r')
    translocations = yaml.load(stream)
    # First translocation - now with alignment
    expected = {
        'xref': '18594527',
        'partners': [
            {'disp': '5prime', 'hgnc_symbol': 'ACSL3',
             'aln': {'breakpoint': 222900780,
                     'query_start': 1,
                     'query_end': 36,
                     'chromosome_name': '2',
                     'start_position': 222900745,
                     'end_position': 222900780,
                     'strand': 1
                     }},
            {'disp': '3prime', 'hgnc_symbol': 'ETV1',
             'aln': {'breakpoint': 13935896,
                     'query_start': 37,
                     'query_end': 71,
                     'chromosome_name': '7',
                     'start_position': 13935862,
                     'end_position': 13935896,
                     'strand': -1
                     }},
            ],
        'junction_sequence':
            ('ctgtgtcacaccaccttagcctcttgatcgaggaagTGCCTATGATCA' +
             'GAAGCCACAAGTGGGAATGAGGC')
        }
    # Compare
    compare_translocations(translocations[0], expected)
