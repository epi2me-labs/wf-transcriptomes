#!/usr/bin/env python
"""Create workflow report."""

import json
from pathlib import Path
import pickle

from bokeh.models import HoverTool, Range1d
from bokeh.models.tickers import AdaptiveTicker
from dominate.tags import li, ul
from dominate.util import raw
import ezcharts as ezc
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Grid, Tabs
from ezcharts.plots.categorical import barplot
from ezcharts.util import get_named_logger
import pandas as pd
from . import de_plots  # noqa: ABS101
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("--report", help="Report output file")
    parser.add_argument(
        "--metadata", default='metadata.json', required=True,
        help="sample metadata")
    parser.add_argument(
        "--stats", nargs='+',
        help="Fastcat per-read stats, ordered as per entries in --metadata.")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--wf_version", default='unknown',
        help="version of the executed workflow")
    parser.add_argument(
        "--alignment_stats", required=False, default=None, type=Path,
        help="TSV summary file of alignment statistics")
    parser.add_argument(
        "--gff_annotation", required=False, type=Path,
        help="transcriptome annotation gff file")
    parser.add_argument(
        "--gffcompare_dir", required=False, default=None, type=Path,
        help="gffcompare outout dir")
    parser.add_argument(
        "--pychop_report", required=False, default=None, type=Path,
        help="TSV summary file of pychopper statistics")
    parser.add_argument(
        "--isoform_table", required=False, type=Path,
        help="Path to directory of TSV files with isoform summaries")
    parser.add_argument(
        "--isoform_table_nrows", required=False, type=int, default=5000,
        help="Maximum rows to display in isoforms table")
    parser.add_argument(
        "--transcriptome_summary", required=False, type=Path,
        help="Path to dir containing transcriptome summary results ")
    parser.add_argument(
        "--de_report", required=False, type=Path, default=None,
        help="Differential expression report optional")
    parser.add_argument(
        "--de_stats", required=False, type=Path, default=None,
        help="Differential expression report optional")
    parser.add_argument(
        "--pval_threshold", required=False, type=float, default=0.01,
        help=(
            "pvalue theshold for inclusion of differentially expressed genes"
            " transcripts in plots "))

    return parser


def gff_compare_plots(report, gffcompare_outdirs):
    """Create various sections and plots in a WfReport.

    :param gffcompare_outdirs: List of output directories from run_gffcompare
    :return: None

    """
    # Plot overview panel:
    with report.add_section("Annotation summary", "Annotation"):

        raw("""The following plots summarize some of the output from
        <a href="https://ccb.jhu.edu/software/stringtie/gffcompare.shtml">gffcompare</a>
        """)

        raw("""<ul>
            <li><b>Totals:</b> Comparison of the number of stringtie-generated
            transcripts, multiexonic transcripts and loci between reference and
            query</li>
            <li><b>Performance:</b> How accurate are the query transcript annotations
            with respect to the reference at various levels</li>
            <li><b>Missed:</b> Features present in the reference, but absent in the
            query</li>
            <li><b>Novel</b>: Features present in the query transcripts,
            but absent in the reference</li>""")

        # Create the 4 bar plots of gffcompare summaries
        gffcomp_fnames = ['Totals.tsv', 'Performance.tsv', 'Missed.tsv', 'Novel.tsv']
        sample_tabs = Tabs()

        for dir_ in sorted(gffcompare_outdirs):
            sample_id = dir_.name  # Get sample ids from the folder name
            with sample_tabs.add_tab(sample_id):
                with Grid(columns=2):
                    for name in gffcomp_fnames:
                        df = pd.read_csv(dir_ / name, sep='\t')
                        if df.empty:
                            raw("""No gffcompare summary for this sample.
                            This could be due to  incompatible reference fasta and gff
                            files""")
                            break
                        plot = barplot(data=df, x='source', y='counts', hue='type')
                        plot._fig.title = name.replace('.tsv', '')
                        EZChart(plot, height='250px')

    # Plot overlaps panel:
    with report.add_section('Query transfrag classification', 'Classes'):

        raw("""Summaries of the classes assigned by
        <a href="https://ccb.jhu.edu/software/stringtie/gffcompare.shtml">
        gffcompare</a>, which describe the relationship between query transfrag and
        the most similar reference transcript.

        <a href="https://ccb.jhu.edu/software/stringtie/gffcompare_codes.png">
        This diagram</a> illustrates the different classes.
        """)
        tabs = Tabs()
        for dir_ in sorted(gffcompare_outdirs):
            sample_id = dir_.name  # Get sample ids from the folder name
            with tabs.add_tab(sample_id):
                with Grid(columns=2):
                    track_file = dir_ / 'tracking_summary.tsv'

                    df_tracking = (
                        pd.read_csv(track_file, sep="\t", index_col=0)
                        .sort_values('Count', ascending=True))
                    df_tracking.columns = [x.capitalize() for x in df_tracking.columns]
                    # Save class name for table
                    df_tracking.rename(columns={'Class': 'Class_tmp'}, inplace=True)
                    # Create new class name for barplot Description:Class
                    df_tracking['Class'] = \
                        df_tracking['Description'] + ":" + df_tracking['Class_tmp']
                    plot = barplot(
                        data=df_tracking, x='Count', y='Class',  orient='h')
                    # revert to original class name
                    df_tracking = (
                        df_tracking
                        .drop(columns=['Class'])
                        .rename(columns={'Class_tmp': 'Class'}))
                    plot._fig.title.text = sample_id
                    plot._fig.yaxis.major_label_text_font_size = "10pt"
                    # Should the height of the plot be based on number of plot rows?
                    EZChart(plot, height='550px')

                    DataTable.from_pandas(
                        df_tracking.sort_values('Count', ascending=False),
                        use_index=False, paging=False, searchable=False)


def pychopper_plots(report, pychop_report):
    """Make plots from pychopper output.

    :param report: ezcharts report
    :param pychop_report: path to pychopper stats file
    """
    with report.add_section("Pychopper summary statisitcs", "Pychopper"):
        raw("""The following plots summarize
        <a href="https://github.com/epi2me-labs/pychopper">pychopper</a> output""")

        ul(
            li("""Full length: Reads with primers found in correct orientation at
            both ends."""),
            li("""Rescued: subreads extracted from chimeric reads"""),
            li("""Unusable: Reads with missing or incorrect primer orientation"""),
            li(""""+/-: Orientation of reads relative to the mRNA""")
        )

        df = pd.read_csv(pychop_report, sep='\t', index_col=0)

        with Grid(columns=2):
            for sample, df in df.groupby('sample_id'):
                df = df.set_index('Name', drop=True)
                df = df.loc[['Primers_found', 'Rescue', 'Unusable', '+', '-']]
                df = df.rename(
                    index={
                        'Primers_found': 'full length',
                        'Rescue': 'Rescued'
                    },
                    columns={
                        'Value': 'n_reads'
                    }
                ).reset_index(drop=False)
                plot = barplot(
                    data=df, x='Name', y='n_reads')
                plot._fig.title.text = f'{sample}'
                plot._fig.xaxis.axis_label = ""
                EZChart(plot, height="300px")


def transcript_table(report, isoform_table, max_rows):
    """Create searchable table of transcripts.

    :param isoform_table: path to folder of isoform table files
    """
    with report.add_section("Isoforms", "Isoforms"):
        raw("""
        This table details each isoforms identified per sample.

        Table interactivity can be slow if too many isoforms are loaded.
        The number of isoform rows to load in this table can be set with
        <b>isoform_table_nrows</b>. It is currently set to {}.
        """.format(max_rows))
        tabs = Tabs()
        # drop some columns for the big table and do some filtering

        dfs = [pd.read_csv(x, sep='\t') for x in Path(isoform_table).iterdir()]
        df = pd.concat(dfs)
        for sample, df in df.groupby('sample_id'):
            with tabs.add_tab(sample):
                # Keep top n rows with most coverage
                df.sort_values('cov', ascending=False, inplace=True)
                df['cov'] = df['cov'].astype(int)
                df = df.iloc[0: max_rows, :]

                # Sort by transcripts with the highest isoform diversity
                df.sort_values('parent gene iso num', inplace=True, ascending=False)
                DataTable.from_pandas(df, use_index=False)


def transcriptome_summary(report, summaries_dir):
    """Plot transcriptome summaries.

    This section consists of four plots
    1: Isoforms per gene histogram
    2: Exons per transcript histogram
    3: Transcript lengths box plot
    4: Transcriptome summary table
    """
    with report.add_section("Transcriptome summary", 'Summary'):

        plot_height = "300px"
        tabs = Tabs()
        data = {}
        # Load all the dataframes upfront to get sample_id for sorting.
        for summ_file in summaries_dir.glob('summary_*.pkl'):
            with open(summ_file, 'rb') as fh:
                summ = pickle.load(fh)
            sample_id = summ['sample_id']
            data[sample_id] = summ

        for sample_id in sorted(data):
            summ = data[sample_id]
            with tabs.add_tab(sample_id):
                with Grid(columns=4):

                    df_isoforms_per_gene = pd.DataFrame(
                        summ['isoforms_per_gene'].items(),
                        columns=['n_isoforms', 'count'])
                    iso_per_gene_plt = ezc.histplot(
                        data=df_isoforms_per_gene['n_isoforms'],
                        weights=df_isoforms_per_gene['count'],
                        discrete=df_isoforms_per_gene['n_isoforms'].max() < 20,
                        bins=20)
                    iso_per_gene_plt._fig.title = "Isoforms per gene"
                    iso_per_gene_plt._fig.xaxis.axis_label = 'Number of isoforms'
                    iso_per_gene_plt._fig.yaxis.axis_label = 'Number of genes'
                    iso_per_gene_plt._fig.xaxis.ticker = \
                        AdaptiveTicker(min_interval=1)
                    EZChart(iso_per_gene_plt, height=plot_height)

                    df_exons_per_transcript = pd.DataFrame(
                        summ['exons_per_transcript'].items(),
                        columns=['n_exons', 'count'])
                    exons_per_tr_plt = ezc.histplot(
                        data=df_exons_per_transcript['n_exons'],
                        weights=df_exons_per_transcript['count'],
                        discrete=df_exons_per_transcript['n_exons'].max() < 20,
                        bins=20)
                    exons_per_tr_plt._fig.title = "Exons per transcript"
                    exons_per_tr_plt._fig.xaxis.axis_label = 'Number of exons'
                    exons_per_tr_plt._fig.yaxis.axis_label = 'Number of transcripts'
                    exons_per_tr_plt._fig.xaxis.ticker = \
                        AdaptiveTicker(min_interval=1)
                    EZChart(exons_per_tr_plt, height=plot_height)

                    df_tr_len = pd.DataFrame.from_dict({
                        'transcript_lengths': summ['transcript_lengths']})
                    df_tr_len['group'] = sample_id
                    tr_lens_plt = ezc.boxplot(
                        data=df_tr_len, x='group', y='transcript_lengths')
                    tr_lens_plt._fig.title = "Transcript lengths"
                    tr_lens_plt._fig.y_range = (
                        Range1d(0, df_tr_len['transcript_lengths'].max()))

                    # Disable hover tool as ezcharts sometimes produces strange
                    # results for boxplots
                    hover_tools = tr_lens_plt._fig.select(dict(type=HoverTool))
                    for hover_tool in hover_tools:
                        tr_lens_plt._fig.tools.remove(hover_tool)
                    EZChart(tr_lens_plt, height=plot_height)

                    df_table = pd.DataFrame.from_dict(summ['summaries']).T
                    df_table = df_table.reset_index(drop=False)
                    df_table.columns = ['Transcriptome summary',  '']
                    DataTable.from_pandas(
                        df_table, use_index=False, searchable=False, paging=False)


def de_section(report, de_report_dir, de_aln_stats_dir, pval_threshold):
    """Make differential transcript expression section."""
    dexseq = de_report_dir / "results_dexseq.tsv"
    dge = de_report_dir / "results_dge.tsv"
    dtu = de_report_dir / "results_dtu_stageR.tsv"
    # GFF file can have gtf or gff extension.
    # Will be the original (transcriptome_source=precomputed) or wf-assembled annotation
    annotation = next(de_report_dir.glob("*.g*f*"))
    tpm = de_report_dir / "unfiltered_tpm_transcript_counts.tsv"
    filtered = de_report_dir / "filtered_transcript_counts_with_genes.tsv"
    unfiltered = de_report_dir / "unfiltered_transcript_counts_with_genes.tsv"
    gene_counts = de_report_dir / "all_gene_counts.tsv"
    # This will also add a gene name column to the above counts tsv files
    de_plots.de_section(
        annotation=annotation,
        dexseq=dexseq,
        dge=dge,
        dtu=dtu,
        tpm=tpm,
        report=report,
        filtered=filtered,
        unfiltered=unfiltered,
        gene_counts=gene_counts,
        aln_stats_dir=de_aln_stats_dir,
        pval_threshold=pval_threshold
    )


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    logger.info('Building report')

    report = LabsReport(
        'Workflow transcriptomes report', 'wf-transcriptomes',
        args.params, args.versions, args.wf_version,
        head_resources=[*LAB_head_resources])

    with open(args.metadata) as metadata:
        sample_details = [{
            'sample': d['alias'],
            'type': d['type'],
            'barcode': d['barcode']
        } for d in json.load(metadata)]

    with report.add_section('Read summary', 'Read summary'):
        names = tuple(d['sample'] for d in sample_details)
        stats = tuple(args.stats)
        if len(stats) == 1:
            stats = stats[0]
            names = names[0]
        fastcat.SeqSummary(stats, sample_names=names)

    if args.alignment_stats is not None:
        with report.add_section('Alignment summary', 'Alignment'):
            raw("""Reference genome alignment statistics. These were generated using
              <a href="https://bioinf.shenwei.me/seqkit/">seqkit</a> )
              `seqkit bam -s` """)

            stats_dfs = []
            for stats_file in args.alignment_stats.iterdir():
                df = pd.read_csv(stats_file, sep='\\s+')
                stats_dfs.append(df)
            aln_stats_df = pd.concat(stats_dfs).sort_values(by='sample_id')
            DataTable.from_pandas(aln_stats_df, use_index=False, export=True)

    if args.pychop_report is not None:
        # report is single file concatenated for all samples
        pychopper_plots(report, next(args.pychop_report.iterdir()))

    # Do upstream
    if args.gff_annotation is not None:
        transcriptome_summary(report, args.transcriptome_summary)

    if args.gffcompare_dir is not None:
        gff_compare_plots(
            report,
            [x for x in Path(args.gffcompare_dir).iterdir()])

    if args.isoform_table is not None:
        transcript_table(report, args.isoform_table, args.isoform_table_nrows)

    if args.de_report:
        de_section(report, args.de_report, args.de_stats, pval_threshold=0.01)

    report.write(args.report)
    logger.info('Report writing finished')
