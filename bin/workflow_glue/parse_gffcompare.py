#!/usr/bin/env python
"""Make report tables and data for plotting."""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("Parse gffcompare")
    parser.add_argument(
        '--sample_id', help="Sample ID", required=True)
    parser.add_argument(
        '--gffcompare_dir',
        help="The gffcompare output directory",
        required=False,
        type=Path)
    parser.add_argument(
        '--isoform_table_out',
        help="Output path for per-isoform table",
        type=Path)
    parser.add_argument(
        '--tracking',
        help="gffcompare tracking file",
        type=Path)
    parser.add_argument(
        "--annotation",
        required=False,
        default=None, help="Reference annotation GFF file")

    return parser


def _parse_stat_line(sl):
    """Parse a stats line."""
    res = {}
    tmp = sl.split(':')[1].split('|')
    res['sensitivity'] = float(tmp[0].strip())
    res['precision'] = float(tmp[1].strip())
    return res


def _parse_matching_line(line):
    """Parse a matching line."""
    tmp = line.split(':')[1].strip()
    return int(tmp)


def _parse_mn_line(line):
    """Parse a miss or novel line."""
    res = {}
    tmp = line.split(':')[1].strip()
    tmp = tmp.split('/')
    res['value'] = int(tmp[0])
    tmp = tmp[1].split('(')
    res['value_total'] = int(tmp[0].strip())
    res['percent'] = float(tmp[1].split('%)')[0])
    return res


def _parse_total_line(line):
    """Parse a total line."""
    res = {}
    tmp = line.split(':')[1].strip()
    tmp = tmp.split('in')
    res['transcripts'] = int(tmp[0].strip())
    tmp = tmp[1].split('loci')
    res['loci'] = int(tmp[0].strip())
    tmp = int(tmp[1].split('(')[1].split(' ')[0])
    res['me_transcripts'] = tmp
    return res


def parse_gffcmp_stats(gffcompare_stats, sample_id, outpath):
    """Parse a gffcompare stats file.

    Gffcompare stats file
    :param gffcompare_stats: Path to the gffcompare stats file.
    :returns: Return as tuple of dataframes containing:
    perfromance statistics, match statistics, miss statistics,
    novel statistics, total statistics.
    :rtype: tuple
    """
    performance = []
    missed = []
    novel = []
    total = []

    with open(gffcompare_stats, 'r') as fh:
        for line in fh:
            line = line.strip()
            if len(line) == 0:
                continue

            # Parse totals:
            if line.startswith('#     Query mRNAs'):
                r = _parse_total_line(line)
                total.append([r['loci'], 'loci', 'query'])
                total.append([r['transcripts'], 'transcripts', 'query'])
                total.append([r['me_transcripts'], 'multexonic', 'query'])
            if line.startswith('# Reference mRNAs '):
                r = _parse_total_line(line)
                total.append([r['loci'], 'loci', 'reference'])
                total.append([r['transcripts'], 'transcripts', 'reference'])
                total.append([r['me_transcripts'], 'multexonic', 'reference'])

            # Parse basic statistics:
            if line.startswith('Base level'):
                st = _parse_stat_line(line)
                performance.append((st['sensitivity'], 'Sensitivity', 'Base'))
                performance.append((st['precision'], 'Precision', 'Base'))
            if line.startswith('Exon level'):
                st = _parse_stat_line(line)
                performance.append((st['sensitivity'], 'Sensitivity', 'Exon'))
                performance.append((st['precision'], 'Precision', 'Exon'))
            if line.startswith('Intron level'):
                st = _parse_stat_line(line)
                performance.append((st['sensitivity'], 'Sensitivity', 'Intron'))
                performance.append((st['precision'], 'Precision', 'Intron'))
            if line.startswith('Intron chain level'):
                st = _parse_stat_line(line)
                performance.append((st['sensitivity'], 'Sensitivity', 'Intron_chain'))
                performance.append((st['precision'], 'Precision', 'Intron_chain'))
            if line.startswith('Transcript level'):
                st = _parse_stat_line(line)
                performance.append((st['sensitivity'], 'Sensitivity', 'Transcript'))
                performance.append((st['precision'], 'Precision', 'Transcript'))
            if line.startswith('Locus level'):
                st = _parse_stat_line(line)
                performance.append((st['sensitivity'], 'Sensitivity', 'Locus'))
                performance.append((st['precision'], 'Precision', 'Locus'))

            # Parse missing statistics:
            if line.startswith('Missed exons'):
                r = _parse_mn_line(line)
                missed.append((r['value'], 'Missed', 'Exons'))
                missed.append((r['value_total'], 'total', 'Exons'))
                missed.append((r['percent'], 'Percent', 'Exons'))
            if line.startswith('Missed introns'):
                r = _parse_mn_line(line)
                missed.append((r['value'], 'Missed', 'Introns'))
                missed.append((r['value_total'], 'total', 'Introns'))
                missed.append((r['percent'], 'Percent', 'Introns'))
            if line.startswith('Missed loci'):
                r = _parse_mn_line(line)
                missed.append((r['value'], 'Missed', 'Loci'))
                missed.append((r['value_total'], 'total', 'Loci'))
                missed.append((r['percent'], 'Percent', 'Loci'))

            # Parse novel statistics:
            if line.startswith('Novel exons'):
                r = _parse_mn_line(line)
                novel.append((r['value'], 'Novel', 'Exons'))
                novel.append((r['value_total'], 'Total', 'Exons'))
                novel.append((r['percent'], 'Percent_novel', 'Exons'))
            if line.startswith('Novel introns'):
                r = _parse_mn_line(line)
                novel.append((r['value'], 'Novel', 'Introns'))
                novel.append((r['value_total'], 'Total', 'Introns'))
                novel.append((r['percent'], 'Percent_novel', 'Introns'))
            if line.startswith('Novel loci'):
                r = _parse_mn_line(line)
                novel.append((r['value'], 'Novel', 'Loci'))
                novel.append((r['value_total'], 'Total', 'Loci'))
                novel.append((r['percent'], 'Percent_novel', 'Loci'))

    def write_records(records, fn):
        pd.DataFrame.from_records(records, columns=['counts', 'type', 'source']) \
            .to_csv(outpath / fn, sep='\t')

    write_records(total, 'Totals.tsv')
    write_records(missed, 'Missed.tsv')
    write_records(performance, 'Performance.tsv')
    write_records(novel, 'Novel.tsv')


def tracking_summary(tracking_file, output_dir, annotations=None):
    """Write per transcript class gffcompare tracking files."""
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

        df = (
            pd.DataFrame(tracking['class'].value_counts())
            .reset_index()
            .rename(columns={'index': 'class', 'class': 'Count'})
        )

        df['Percent'] = round(df['Count'] * 100 / df['Count'].sum(), 2)
        df['description'] = [nice_names[x] for x in df['class']]

        df = df.sort_values('Count', ascending=True)
        df.to_csv(output_dir / 'tracking_summary.tsv', sep='\t')

    else:
        logger = get_named_logger('trackingSum')
        logger.info("Skipping classification summary as no annotation provided.")


def make_isoform_table(gffcompare_dir, sample_id, outpath):
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

    if df.empty:  # No transcripts. Write a header only result file
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

        df.to_csv(outpath, sep='\t', index=False)


def main(args):
    """Entry point."""
    if args.gffcompare_dir:  # TODO: should this every be optional?
        stats = args.gffcompare_dir / 'str_merged.stats'
        parse_gffcmp_stats(stats, args.sample_id, args.gffcompare_dir)
        make_isoform_table(args.gffcompare_dir, args.sample_id, args.isoform_table_out)
        tracking_summary(
            args.tracking, args.gffcompare_dir, args.annotation)
