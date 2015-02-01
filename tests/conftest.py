import pytest

def pytest_addoption(parser):
    parser.addoption("--ticdb_file", action="store", default="",
        help="This input is the tab-delimited TICdb output file.")
    parser.addoption("--distance", action="store", default="100",
        help="This input is the threshold distance in kb to determine \
        if two fragile sites are from the same fragile region.")
    parser.addoption("--output_file", action="store", default="",
        help="This input is the name of the .csv output file.")

@pytest.fixture
def ticdb_file(request):
    """A mock input file of translocation sites"""
    return request.config.getoption("--ticdb_file")
    #return 'fragile_finder_test_site.txt'

@pytest.fixture
def output_file(request):
    """A file like object I can write to and assert against"""
    #return request.config.getoption("--output_file")
    return 'ff_test.csv'

@pytest.fixture
def distance(request):
    """threshold distance between duplicate fragile sites"""
    # return request.config.getoption("--distance")
    return '100'
