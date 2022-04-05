#!/usr/bin/env python
"""Generate per-transcript class sumarrarry files from gffcompare."""

import argparse
import os
import sys

import pandas as pd


def parse_args(argv=sys.argv[1:]):
    """Parse arguments."""
    description = """Script to run the isoform workflow """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--tracking", required=True, help="")
    parser.add_argument("--output_dir", required=True, help="")
    parser.add_argument("--annotation", required=False, default=None, help="")
    return parser.parse_args(argv)


def generate_tracking_summary(tracking_file, output_dir, annotations=None):
    """Write per transcript class gffcompare tracking files."""
    # results/gffcompare/gff_comparestringtie.tracking
    write_empty_tsvs = False
    tracking_headings = [
        "query_transfrag_id", "query_locus_id", "ref_gene_id",
        "class", "details"]
    nice_names = {
        '=': 'complete', 'c': 'contained', 'k': 'containment',
        'm': 'retained', 'n': 'retained (partial)', 'j': 'multi',
        'e': 'single', 'o': 'overlap', 's': 'opposite',
        'x': 'exonic', 'i': 'intron', 'y': 'contains', 'p': 'runon',
        'r': 'repeat', 'u': 'unknown'}

    if os.path.exists(annotations):
        tracking = pd.read_csv(
            tracking_file, sep="\t", names=tracking_headings[1:],
            index_col=0)

        d = pd.DataFrame(tracking['class'].value_counts()) \
            .reset_index().rename(columns={'index': 'class', 'class': 'count'})
        d['description'] = [nice_names[x] for x in d['class']]

        # write a separate table for each class
        for class_code, table in tracking.groupby('class'):
            if not write_empty_tsvs and table.empty:
                print("Skipping: No transcripts found for: {}".format(
                    class_code))
                continue

            path = tracking_file + ".{}.tsv".format(class_code)
            table.to_csv(path)
    else:
        print("Skipping classification summary as no annotation provided.")


def main(args):
    """Run entry point."""
    assert os.path.isfile(args.tracking)
    assert os.path.isdir(args.output_dir)
    if args.annotation:
        os.path.isfile(args.annotation)
    generate_tracking_summary(
        args.tracking, output_dir=args.output_dir, annotations=args.annotation)


if __name__ == '__main__':
    main(parse_args(sys.argv[1:]))
