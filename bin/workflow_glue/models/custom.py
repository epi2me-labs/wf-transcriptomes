"""Extended models for the workflow."""
from dataclasses import dataclass

from workflow_glue.models.common import CheckResult


@dataclass
class CheckResult(CheckResult):
    """
    A result of some check the workflow has performed.

    This can be at a sample or workflow level.
    """

    categories = dict(
        example_check="Example check category")
