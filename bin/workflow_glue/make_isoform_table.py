#!/usr/bin/env python
"""Make report tables and data for plotting."""
from pathlib import Path

import numpy as np
import pandas as pd
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("Prepare report data")
    parser.add_argument(
        '--sample_id', help="Sample ID", required=True)
    parser.add_argument(
        '--gffcompare_dir',
        help="The gffcompare output directory",
        required=False,
        type=Path)

    return parser


def make_isoform_table(gffcompare_dir, sample_id):
    """Make an isoform summary table."""
    try:
        tmap_file = next(gffcompare_dir.glob('*.tmap'))
    except StopIteration:
        raise ValueError("Cannot find .tmap file in {}".format(gffcompare_dir))
    dtypes = {
        'ref_gene_id': str,
        'ref_id': str,
        'class_code': str,
        'qry_id': str,
        'num_exons': np.uint16,
        'cov': np.uint32,
        'len': np.uint32
    }
    df = pd.read_csv(
        tmap_file, sep='\t+',
        index_col=None,
        usecols=list(dtypes.keys()),
        dtype=dtypes)

    if len(df) == 0:  # No transcripts. Write a header only result file
        df = pd.DataFrame(
            columns=list(dtypes.keys()) + ['sample_id', 'parent gene iso num'])
        df.to_csv(f'{sample_id}_transcripts_table.tsv', sep='\t', index=False)
    else:
        df = df.assign(sample_id=sample_id)

        #  Make a column of number of isoforms in parent gene
        gb = df.groupby(['ref_gene_id']).count()
        gb.rename(columns={'ref_id': 'num_isoforms'}, inplace=True)

        df['parent gene iso num'] = df.apply(
            lambda x: gb.loc[(x.ref_gene_id), 'num_isoforms'], axis=1)

        # Unclassified transcripts should not be lumped together
        df.loc[df.class_code == 'u', 'parent gene iso num'] = None

        df.to_csv(f'{sample_id}_transcripts_table.tsv', sep='\t', index=False)


def main(args):
    """Entry point."""
    if args.gffcompare_dir:
        make_isoform_table(args.gffcompare_dir, args.sample_id)
