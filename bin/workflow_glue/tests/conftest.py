#!/usr/bin/env python
"""Pytests argument definitions."""

from pathlib import Path
import sys


# The CI template invokes `pytest /host/bin/workflow_glue/tests`, so add
# `/host/bin` explicitly to keep imports stable regardless of the working
# directory pytest picks.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def pytest_addoption(parser):
    """Define command line arguments for pytest."""
    parser.addoption(
        "--test_data",
        action="store",
        default="/host/test_data"
    )
