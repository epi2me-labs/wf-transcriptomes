#!/usr/bin/env python
"""Get fastq QC reports."""

# -*- coding: utf-8 -*-

import argparse
import os
import sys

import numpy as np
from pysam import FastxFile


def parse_args(argv=sys.argv[1:]):
    """Parse args."""
    description = """Script to run the isoform workflow """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--fastq", required=True, help="")
    parser.add_argument("--output_dir", required=True, help="")
    return parser.parse_args(argv)


def run_fastq_qc(fastq_path, output):
    """Write QC info to files."""
    qualities = list()
    mean_qualities = list()
    lengths = list()

    with FastxFile(fastq_path) as fq:
        for rec in fq:
            # ONT calculation for "mean Q score"
            quals = np.fromiter(
                (ord(x) - 33 for x in rec.quality),
                dtype=int, count=len(rec.quality))
            mean_p = np.mean(np.power(10, quals / -10))
            mean_qualities.append(-10 * np.log10(mean_p))
            # all qualities
            qualities.extend(quals)
            lengths.append(len(quals))

    with open(os.path.join(output, "base_qual.txt"), 'w') as f:
        f.write("\n".join((str(q) for q in qualities)))

    with open(os.path.join(output, "read_qual.txt"), 'w') as f:
        f.write("\n".join((str(q) for q in mean_qualities)))

    with open(os.path.join(output, "lengths.txt"), 'w') as f:
        f.write("\n".join((str(_l) for _l in lengths)))


def main(args):
    """Run entry point."""
    assert os.path.isfile(args.fastq)
    assert os.path.isdir(args.output_dir)
    run_fastq_qc(fastq_path=args.fastq, output=args.output_dir)


if __name__ == '__main__':
    main(args=parse_args())
