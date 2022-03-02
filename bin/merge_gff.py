#!/usr/bin/env python

"""
Merge and fix gff files.

Merge multiple gff files into single file.
Rename gene and transcript ids to avoid attribute conflicts from
independently-created files.
"""
import argparse
from pathlib import Path
import re

from natsort import natsorted


def main(gff_files: str, outfile: str):
    """Entry point."""
    regx_id = re.compile(r'gene_id "STRG\.(\d+)"')

    start = 1
    with open(outfile, 'w') as fh:
        for gff in gff_files:

            text = Path(gff).read_text()
            if start != 1:
                # Strip headers
                text = [x for x in text.splitlines() if not x.startswith('#')]
                text = '\n'.join(text)

            ids = natsorted(set(re.findall(regx_id, text)))
            new_gene_ids = list(range(start, start + len(ids)))
            id_map = dict(zip(ids, new_gene_ids))

            for old_id, new_id in id_map.items():
                text = text.replace(
                    f'gene_id "STRG.{old_id}"',
                    f'gene_id "STRG.{new_id}"')
                text = text.replace(
                    f'transcript_id "STRG.{old_id}.',
                    f'transcript_id "STRG.{new_id}.')
            fh.write(text)
            fh.write('\n')
            start += len(ids)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_files", help="gff files to merge",
                        required=True, nargs='+')
    parser.add_argument("--out_file", help="where to save merged files",
                        required=True)
    args = parser.parse_args()
    main(args.gff_files, args.out_file)
