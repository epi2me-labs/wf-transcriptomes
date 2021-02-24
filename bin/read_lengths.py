#!/usr/bin/env python

import argparse
import pysam


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta')
    parser.add_argument('output')
    args = parser.parse_args()

    with open(args.output, 'w') as fh:
        for rec in pysam.FastxFile(args.fasta):
            fh.write("{}\t{}\n".format(rec.name, len(rec.sequence)))

if __name__ == '__main__':
    main()
