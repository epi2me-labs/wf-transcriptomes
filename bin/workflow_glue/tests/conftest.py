#!/usr/bin/env python
"""Pytests argument definitions."""


def pytest_addoption(parser):
    """Define command line arguments for pytest."""
    parser.addoption(
        "--test_data",
        action="store",
        default="/host/test_data"
    )
