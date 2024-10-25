#!/usr/bin/env python
"""Create de report section."""

import os

from dominate.tags import h5, p
from dominate.util import raw
from ezcharts import scatterplot
from ezcharts.components.ezchart import EZChart
from ezcharts.layout.snippets import DataTable
import numpy as np
import pandas as pd


def parse_seqkit(fname):
    """Get seqkit columns."""
    cols = {
        'Read': str, 'Ref': str, 'MapQual': int, 'Acc': float, 'ReadLen': int,
        'ReadAln': int, 'ReadCov': float, 'MeanQual': float,
        'IsSec': bool, 'IsSup': bool}
    df = pd.read_csv(fname, sep="\t", dtype=cols, usecols=cols.keys())
    df['Clipped'] = df['ReadLen'] - df['ReadAln']
    df['Type'] = 'Primary'
    df.loc[df['IsSec'], 'Type'] = 'Secondary'
    df.loc[df['IsSup'], 'Type'] = 'Supplementary'
    df["fname"] = os.path.basename(fname).rstrip(".seqkit.stats")
    return df


def number_of_alignments(df, field_name):
    """Group alignments for summary table."""
    grouped = df.groupby('fname').agg(**{
        field_name: ('Read', 'size'),
    })
    return grouped.transpose()


def create_summary_table(df):
    """Create summary table."""
    all_aln = number_of_alignments(df, "Read mappings")
    primary = number_of_alignments(df.loc[df['Type'] == 'Primary'], "Primary")
    secondary = number_of_alignments(
        df.loc[df['Type'] == 'Secondary'], "Secondary")
    supplementary = number_of_alignments(
        df.loc[df['Type'] == 'Supplementary'], "Supplementary")
    avg_acc = df.loc[df['Type'] == 'Primary'].groupby(
        'fname').agg(**{"Median Qscore": ('MeanQual', 'median'), }).transpose()
    avg_mapq = df.loc[df['Type'] == 'Primary'].groupby(
        'fname').agg(**{"Median MAPQ": ('MapQual', 'median'), }).transpose()
    return pd.concat([
        all_aln, primary, secondary, supplementary,
        avg_acc, avg_mapq])


def dexseq_section(dexseq_file, id_dic, pval_thresh):
    """Add gene isoforms table and plot."""
    h5("Differential Isoform usage")
    p("""Table showing gene isoforms, ranked by adjusted
    p-value, from the DEXSeq analysis. Information shown includes the log2 fold
    change between experimental conditions, the log-scaled transcript
    abundance and the false discovery corrected p-value (FDR - Benjamini-Hochberg) .
    This table has not been filtered
    for genes that satisfy statistical or magnitudinal thresholds""")

    dexseq_results = pd.read_csv(dexseq_file, sep='\t')
    dexseq_results.index.name = "gene_id:transcript_id"
    # Replace gene id with more useful gene name where possible
    dexseq_results.index = dexseq_results.index.map(
        lambda x: str(id_dic.get(x.split(':')[0])) + ':' + str(x.split(':')[1]))

    DataTable.from_pandas(
        dexseq_results.sort_values(by='pvalue', ascending=True), use_index=True)

    p(
        """The figure below presents the MA plot from the DEXSeq analysis.
        M is the log2 ratio of isoform transcript abundance between conditions.
        A is the log2 transformed mean abundance value.
        Transcripts that satisfy the logFC and FDR-corrected
        (False discovery rate - Benjamini-Hochberg) p-value
        thresholds defined are shaded as 'Up-' or 'Down-' regulated.""")

    dexseq_results['direction'] = 'not_sig'

    dexseq_results.loc[
        (dexseq_results["Log2FC"] > 0) & (dexseq_results['pvalue'] < pval_thresh),
        'direction'] = 'up'

    dexseq_results.loc[
        (dexseq_results["Log2FC"] <= 0) & (dexseq_results['pvalue'] < pval_thresh),
        'direction'] = 'down'

    plot = scatterplot(
        data=dexseq_results, x='Log2MeanExon', y='Log2FC', hue='direction',
        palette=['#E32636', '#7E8896', '#0A22DE'],
        hue_order=['up', 'down', 'not_sig'], marker='circle')
    plot._fig.xaxis.axis_label = "A (log2 transformed mean exon read counts)"
    plot._fig.yaxis.axis_label = "M (log2 transformed differential abundance)"
    plot.legend = dict(orient='horizontal', top=30)
    plot._fig.title = "Average copy per million (CPM) vs Log-fold change (LFC)"
    EZChart(plot)


def dtu_section(dtu_file, gt_dic, ge_dic):
    """Plot dtu section."""
    dtu_results = pd.read_csv(dtu_file, sep='\t')
    dtu_results["gene_name"] = dtu_results["txID"].apply(
        lambda x: gt_dic.get(x))
    dtu_results["geneID"] = dtu_results["geneID"].apply(
        lambda x: ge_dic.get(x))
    dtu_pvals = dtu_results.sort_values(by='gene', ascending=True)
    raw("""Table showing gene and transcript identifiers
    and their FDR-corrected (False discovery rate - Benjamini-Hochberg) probabilities
    for the genes and their isoforms that have been
    identified as showing DTU using the R packages DEXSeq and StageR.
    This list has been shortened requiring that both gene and transcript
    must satisfy the p-value
    threshold""")
    DataTable.from_pandas(dtu_results.loc[dtu_pvals.index], use_index=False)

    raw("""View dtu_plots.pdf file to see plots of differential isoform usage""")


def dge_section(dge_file, ids_dic, pval_thresh):
    """Create DGE table and MA plot."""
    h5("Differential gene expression")
    dge_results = pd.read_csv(dge_file, sep='\t')
    dge_results[['logFC', 'logCPM', 'F']] = dge_results[
        ['logFC', 'logCPM', 'F']].round(2)

    p("""Table showing the genes from the edgeR analysis.
    Information shown includes the log2 fold change between
    experimental conditions, the log-scaled counts per million measure of abundance
    and the FDR-corrected p-value (False discovery rate - Benjamini-Hochberg).
    This table has not been
    filtered for genes that satisfy statistical or magnitudinal thresholds""")

    dge_results.index = dge_results.index.map(lambda x: ids_dic.get(x))
    dge_results = dge_results.sort_values('FDR', ascending=True)
    dge_results.index.name = 'Transcript'
    DataTable.from_pandas(dge_results, use_index=True)

    h5("Results of the edgeR Analysis.")

    p("""This plot visualises differences in measurements between the
    two experimental conditions. M is the log2 ratio of gene expression
    calculated between the conditions.
    A is a log2 transformed mean expression value.
    The figure below presents the MA figure from this edgeR analysis.
    Genes that satisfy the logFC and FDR-corrected
    (False discovery rate - Benjamini-Hochberg) p-value thresholds
    defined are shaded as 'Up-' or 'Down-' regulated.
    """)

    dge = pd.read_csv(dge_file, sep="\t")
    dge['sig'] = None
    dge.loc[(dge["logFC"] > 0) & (dge['PValue'] < pval_thresh), 'sig'] = 'up'
    dge.loc[(dge["logFC"] <= 0) & (dge['PValue'] < pval_thresh), 'sig'] = 'down'
    dge.loc[(dge["PValue"] >= pval_thresh), 'sig'] = 'not_sig'

    plot = scatterplot(
        data=dge, x='logCPM', y='logFC', hue='sig',
        palette=['#E32636', '#7E8896', '#0A22DE'],
        hue_order=['up', 'not_sig', 'down'], marker='circle')
    plot._fig.x_range.start = 10
    plot._fig.xaxis.axis_label = "Average log CPM"
    plot._fig.yaxis.axis_label = "Log-fold change"
    plot.legend = dict(orient='horizontal', top=30)
    # Should opacity of the symbols be lowered?
    plot._fig.title = "Average copy per million (CPM) vs Log-fold change (LFC)"
    EZChart(plot)


def salmon_table(salmon_counts):
    """Create salmon counts summary table."""
    salmon_counts = pd.read_csv(salmon_counts, sep='\t')
    salmon_counts.set_index("Reference", drop=True, append=False, inplace=True)
    salmon_size_top = salmon_counts.sum(axis=1).sort_values(ascending=False)
    salmon_counts = salmon_counts.applymap(np.int64)
    h5("Transcripts Per Million")
    p("""Table showing the annotated Transcripts Per Million
    identified by Minimap2 mapping and Salmon transcript
    detection. Displaying the top 100 transcripts with the highest
    number of mapped reads""")

    salmon_counts = salmon_counts[sorted(salmon_counts.columns)]
    DataTable.from_pandas(
        salmon_counts.loc[salmon_size_top.index].head(n=100), use_index=True)


def get_translations(gtf):
    """Create dict with gene_name and gene_references."""
    with open(gtf) as fh:
        gene_txid = {}
        gene_geid = {}
        geid_gname = {}

        def get_feature(row, feature):
            return row.split(feature)[1].split(
                ";")[0].replace('=', '').replace("\"", "").strip()

        for i in fh:
            if i.startswith("#"):
                continue
            # Different gtf/gff formats contain different attributes
            # and different formating (eg. gene_name="xyz" or gene_name "xyz")
            gene_name = None
            for var_name in ["gene_name", "gene_id", "gene"]:
                if var_name in i:
                    gene_name = get_feature(i, var_name)
                    break

            if 'ref_gene_id' in i:
                gene_reference = get_feature(i, 'ref_gene_id')
            elif 'gene_id' in i:
                gene_reference = get_feature(i, 'gene_id')
            else:
                gene_reference = gene_name
            if 'transcript_id' in i:
                transcript_id = get_feature(i, 'transcript_id')
            else:
                transcript_id = "unknown"
            if 'gene_id' in i:
                gene_id = get_feature(i, 'gene_id')
            else:
                gene_id = gene_name
            gene_txid[transcript_id] = gene_name
            gene_geid[gene_id] = gene_reference
            geid_gname[gene_reference] = gene_name
    return gene_txid, gene_geid, geid_gname


def de_section(
        stringtie, dge, dexseq, dtu,
        tpm, report, filtered, unfiltered,
        gene_counts, aln_stats_dir, pval_threshold=0.01):
    """Differential expression sections."""
    with report.add_section("Differential expression", "DE"):

        p("""This section shows differential gene expression
        and differential isoform usage. Salmon was used to
        assign reads to individual annotated isoforms defined by
        the GTF-format annotation.
        These counts were used to perform a statistical analysis to identify
        the genes and isoforms that show differences in abundance between
        the experimental conditions.
        Any novel genes or transcripts that do not have relevant gene or transcript IDs
        are prefixed with MSTRG for use in differential expression analysis.
        Find the full sequences of any transcripts in the
        final_non_redundant_transcriptome.fasta file.
        """)
        alignment_stats = pd.concat([parse_seqkit(f) for f in aln_stats_dir.iterdir()])
        alignment_summary_df = create_summary_table(alignment_stats)
        alignment_summary_df = alignment_summary_df.fillna(0).applymap(np.int64)
        h5("Alignment summary stats")
        alignment_summary_df.index.name = "statistic"
        DataTable.from_pandas(alignment_summary_df, use_index=True)

        salmon_table(tpm)
        gene_txid, gene_name, geid_gname = get_translations(stringtie)

        # Add gene names columns to counts files and write out
        # for publishing to user dir.
        df_dge = pd.read_csv(dge, sep='\t')
        df_dge.insert(0, 'gene_name', df_dge.index.map(lambda x: geid_gname.get(x)))
        df_dge.to_csv('results_dge.tsv', index=True, index_label="gene_id", sep="\t")

        # write_dge(gene_counts, geid_gname, "all_gene_counts.tsv")
        df_gene_counts = pd.read_csv(gene_counts, sep='\t')
        df_gene_counts.insert(
                0, 'gene_name', df_gene_counts.index.map(lambda x: geid_gname.get(x)))
        df_gene_counts.to_csv(
            'results_dge.tsv', index=True, index_label="gene_id", sep="\t")

        df_filtered = pd.read_csv(filtered, sep='\t')
        df_filtered.insert(1, "gene_name", df_filtered.gene_id.map(
            lambda x: geid_gname.get(x)))
        df_filtered.to_csv(
            'filtered_transcript_counts_with_genes.tsv', index=False, sep='\t')

        df_unfiltered = pd.read_csv(unfiltered, sep='\t')
        df_unfiltered.insert(1, "gene_name", df_unfiltered.gene_id.map(
            lambda x: geid_gname.get(x)))
        df_unfiltered.to_csv(
            'unfiltered_transcript_counts_with_genes.tsv', index=False, sep='\t')

        df_tpm = pd.read_csv(tpm, sep='\t')
        df_tpm.insert(1, "gene_name", df_tpm.Reference.map(
            lambda x: gene_txid.get(x)))
        df_tpm.to_csv("unfiltered_tpm_transcript_counts.tsv", index=False, sep='\t')

        # Add tables to report
        dge_section(dge, gene_name, pval_threshold)
        dexseq_section(dexseq, gene_name, pval_threshold)
        dtu_section(dtu, gene_txid, gene_name)
