#!/usr/bin/env python
"""Merge salmon output count files."""

from functools import reduce

import numpy as np
import pandas as pd

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("merge_count_tsvs")
    parser.add_argument(
        '-j', metavar='join', help="Join type (outer).", default="outer")
    parser.add_argument(
        '-f', metavar='field',
        help="Join on this field (Reference).", default="Reference")
    parser.add_argument(
        '-o', metavar='out_tsv',
        help="Output tsv (merge_tsvs.tsv).", default="merge_tsvs.tsv")
    parser.add_argument(
        '-z', action="store_true",
        help="Fill NA values with zero.", default=False)
    parser.add_argument(
        '-tpm', type=bool, default=False,
        help="TPM instead of counts")
    parser.add_argument(
        '-tsvs', metavar='input_tsvs', nargs='*',
        help="Input tab separated files.")

    return parser


def main(args):
    """Run entry point."""
    dfs = {x: pd.read_csv(x, sep="\t") for x in args.tsvs}

    ndfs = []
    for x, df in dfs.items():
        # Transform counts to integers:
        if args.tpm:
            df = df.rename(columns={'TPM': 'Count', 'Name': 'Reference'})
        else:
            df = df.rename(columns={'NumReads': 'Count', 'Name': 'Reference'})
        df.Count = np.array(df.Count, dtype=int)
        # Take only non-zero counts:
        df = df[df.Count > 0]
        df = df[["Reference", "Count"]]
        df = df.sort_values(by=["Count"], ascending=False)
        name = x.split('.')[0]
        df = df.rename(columns={'Count': name})
        ndfs.append(df)
    dfs = ndfs

    df_merged = reduce(lambda left, right: pd.merge(
        left, right, on=args.f, how=args.j), dfs)
    if args.z:
        df_merged = df_merged.fillna(0)

    df_merged.to_csv(args.o, sep="\t", index=False)
