"""Test check_sample_sheet.py."""
import os

import pytest
from workflow_glue import check_sample_sheet_condition


# define a list of error messages
ERROR_MESSAGES = [
    ("sample_sheet_1.csv", "There must be only two unique conditions in the condition column of the sample sheet."),  # noqa: E501
    ("sample_sheet_2.csv", "Sample sheet has no condition column which is required for the differential expression subworkflow."),  # noqa: E501
    ("sample_sheet_3.csv", "There must be at least 2 repeats for each condition indicated in the sample sheet."),  # noqa: E501
]


@pytest.fixture
def test_data(request):
    """Define data location fixture."""
    return os.path.join(
        request.config.getoption("--test_data"),
        "workflow_glue",
        "check_sample_sheet_condition")


@pytest.mark.parametrize("sample_sheet_name,error_msg", ERROR_MESSAGES)
def test_check_sample_sheet(
        test_data, sample_sheet_name, error_msg):
    """Test the sample sheets."""
    expected_error_message = error_msg
    sample_sheet_path = f"{test_data}/{sample_sheet_name}"
    args = check_sample_sheet_condition.argparser().parse_args(
        [sample_sheet_path]
    )
    try:
        check_sample_sheet_condition.main(args)
    except SystemExit as e:
        assert str(e) == expected_error_message
