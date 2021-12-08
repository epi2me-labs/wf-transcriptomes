#!/usr/bin/env python
"""Generate CSV of pychopper stats."""

# -*- coding: utf-8 -*-

import argparse
import os
import sys

import pandas as pd


def parse_args(argv=sys.argv[1:]):
    """Parse arguments."""
    description = """Script to run the isoform workflow """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--data", required=True, help="")
    parser.add_argument("--output_dir", required=True, help="")
    return parser.parse_args(argv)


def generate_pychopper_stats(tsv, output):
    """Make CSV of pychopper stats."""
    classified_path = os.path.join(output, "pychopper_stats.csv")
    df = pd.read_csv(tsv, sep="\t", index_col="Name")
    classified = df.loc[df["Category"] == "Classification"]\
        .copy().reset_index().rename(columns={'Name': 'Classification'})
    classified["Percentage"] = \
        100 * classified["Value"] / classified["Value"].sum()
    tuning = df.loc[df["Category"] == "AutotuneSample"]\
        .copy().reset_index().rename(columns={'Name': 'Filter'})
    tuning.to_csv(classified_path)


def main(args):
    """Run entry point."""
    assert os.path.isfile(args.data)
    assert os.path.isdir(args.output_dir)
    generate_pychopper_stats(tsv=args.data, output=args.output_dir)


if __name__ == '__main__':
    main(args=parse_args())
