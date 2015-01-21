"""
Test of the fragile finder
"""

import pytest
import mock

import fragile_finder

@pytest.fixture
def input_file():
    """A mock input file of translocation sites"""

@pytest.fixture
def output_file():
    """A file like object I can write to and assert against"""

@pytest.fixture
def mock_blast_search(request):
    """A fake function we will call instead of the real blast query that
    will return known data"""

@pytest.mark.usefixtures('mock_blast_search')
def test_main(input_file, output_file):
    fragile_finder.main(input_file, output_file, 100)
    output_file.seek()
    results = output_file.readlines()
    assert len(results) > 0  # make more specific
