from nose.tools import *
from brkpoints import main

OUTER_TODAY = None


def setup():
    'Setup function that runs before every test'
    print "SETUP!"


def teardown():
    'Teardown function that runs after every test'
    print "TEAR DOWN!"


def test_basic():
    'Test functions must start with test_*'
    print "BASIC TEST!"


def test_main():
    'Call the "main" function'
    main()
    print "MAIN TEST!"
