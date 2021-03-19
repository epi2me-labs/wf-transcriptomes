#!/usr/bin/env python
"""Create a simple summary of a fastq file."""

import argparse
import glob
import itertools
import os

import numpy as np
import pysam


def mean_qual(quals):
    """Calculate mean quality of a read."""
    qual = np.fromiter(
        (ord(x) - 33 for x in quals),
        dtype=int, count=len(quals))
    mean_p = np.mean(np.power(10, qual / -10))
    return -10 * np.log10(mean_p)


def main():
    """Run entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fastq", help="Directory containing .fastq(.gz) files, or single fastq")
    parser.add_argument(
        "output", help="Output file")
    args = parser.parse_args()

    if os.path.isfile(args.fastq):
        fastqs = [args.fastq]
    elif os.path.isdir(args.fastq):
        fastqs = glob.glob(os.path.join(args.directory, "*.fastq*"))
    else:
        raise IOError("fastq argument should be directory of file.")
    reads = itertools.chain.from_iterable(
        pysam.FastxFile(fname) for fname in fastqs)

    with open(args.output, "w") as fh:
        # names as in Guppy
        fh.write("read_id\tsequence_length_template\tmean_qscore_template\n")
        for read in reads:
            fh.write("\t".join(str(x) for x in (
                read.name, len(read.sequence), mean_qual(read.quality))))
            fh.write("\n")


if __name__ == "__main__":
    main()
