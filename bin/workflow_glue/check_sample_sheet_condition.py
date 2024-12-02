#!/usr/bin/env python
"""Check if a sample sheet is valid."""
from collections import Counter
import csv
import sys

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("checkSheetCondition")
    with open(args.sample_sheet, "r") as f:
        csv_reader = csv.DictReader(f)
        conditions_count = Counter()
        for row in csv_reader:
            if "condition" in row:
                conditions_count[row['condition']] += 1
            else:
                sys.exit(
                    "Sample sheet has no condition column "
                    "which is required for the "
                    "differential expression subworkflow.")
        if len(conditions_count.keys()) != 2:
            sys.exit(
                "There must be only two unique conditions "
                "in the condition column of the sample sheet.")
        if "control" not in conditions_count:
            sys.exit(
                "One of the condition types must be control, "
                "to indicate which samples to use as the reference.")
        if any(v < 2 for v in conditions_count.values()):
            sys.exit(
                "There must be at least 2 repeats for each "
                "condition indicated in the sample sheet.")
    logger.info(f"Checked sample sheet for condition column {args.sample_sheet}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("check_sample_sheet_condition")
    parser.add_argument("sample_sheet", help="Sample sheet to check")
    return parser
