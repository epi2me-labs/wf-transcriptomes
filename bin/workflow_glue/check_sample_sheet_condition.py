#!/usr/bin/env python
"""Check if a sample sheet is valid."""
import csv
import sys

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("checkSheetCondition")
    with open(args.sample_sheet, "r") as f:
        csv_reader = csv.DictReader(f)
        unique_controls = []
        controls_dic = {}
        for row in csv_reader:
            if 'condition' in list(row.keys()):
                unique_controls.append(row['condition'])
                if row['condition'] not in controls_dic:
                    controls_dic[row['condition']] = 1
                else:
                    controls_dic[row['condition']] += 1
            else:
                sys.exit(
                    "Sample sheet has no condition column "
                    "which is required for the "
                    "differential expression subworkflow.")
        if len(list(set(controls_dic.keys()))) != 2:
            sys.exit(
                "There must be only two unique conditions "
                "in the condition column of the sample sheet.")
        for val in list(controls_dic.values()):
            if val < 2:
                sys.exit(
                    "There must be at least 2 repeats for each "
                    "condition indicated in the sample sheet.")

    logger.info(f"Checked sample sheet for condition column {args.sample_sheet}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("check_sample_sheet_condition")
    parser.add_argument("sample_sheet", help="Sample sheet to check")
    return parser
