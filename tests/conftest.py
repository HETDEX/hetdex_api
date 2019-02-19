"""

Setup fixtures for the tests

Taken from Francesco Montesano's
pyhetdex stuff

"""

import os
import py
import pytest


@pytest.fixture(scope="session")
def datadir():
    """ Return a py.path.local object for the test data directory"""
    return py.path.local(os.path.dirname(__file__)).join('data')


