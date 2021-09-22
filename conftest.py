"""Configuration for pytest."""

import pytest
from click.testing import CliRunner


@pytest.fixture
def runner():
    """Return a CliRunner."""
    return CliRunner()


@pytest.fixture(autouse=True)
def add_standard_imports(doctest_namespace):
    """Add py3dep namespace for doctest."""
    import py3dep

    doctest_namespace["py3dep"] = py3dep
