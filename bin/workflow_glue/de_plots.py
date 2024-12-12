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


def dexseq_section(dexseq_file, tr_id_to_gene_name, tr_id_to_gene_id, pval_thresh):
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

    # Replace any occurrences of stringtie-generated MSTRG gene ids with
    # reference gene_ids.
    dexseq_results.index = dexseq_results.index.map(
        lambda ge_tr: str(  # lookup gene_id from transcript_id [1]
            f"{tr_id_to_gene_id.get(ge_tr.split(':')[1])}: {str(ge_tr.split(':')[1])}")
    )

    # Add gene name column.
    dexseq_results.insert(0, "gene_name", dexseq_results.index.map(
        lambda x: tr_id_to_gene_name.get(x.split(':')[1])))

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


def dtu_section(dtu_file, txid_to_gene_name):
    """Plot dtu section."""
    dtu_results = pd.read_csv(dtu_file, sep='\t')
    dtu_results["gene_name"] = dtu_results["txID"].apply(
        lambda x: txid_to_gene_name.get(x))

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


def dge_section(df, pval_thresh):
    """Create DGE table and MA plot."""
    h5("Differential gene expression")
    df[['logFC', 'logCPM', 'F']] = df[
        ['logFC', 'logCPM', 'F']].round(2)

    p("""Table showing the genes from the edgeR analysis.
    Information shown includes the log2 fold change between
    experimental conditions, the log-scaled counts per million measure of abundance
    and the FDR-corrected p-value (False discovery rate - Benjamini-Hochberg).
    This table has not been
    filtered for genes that satisfy statistical or magnitudinal thresholds""")

    df = df.sort_values('FDR', ascending=True)
    df.index.name = 'gene_id'
    DataTable.from_pandas(df, use_index=True)

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
    df['sig'] = None
    df.loc[(df["logFC"] > 0) & (df['PValue'] < pval_thresh), 'sig'] = 'up'
    df.loc[(df["logFC"] <= 0) & (df['PValue'] < pval_thresh), 'sig'] = 'down'
    df.loc[(df["PValue"] >= pval_thresh), 'sig'] = 'not_sig'

    plot = scatterplot(
        data=df, x='logCPM', y='logFC', hue='sig',
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
    """Create gene_and transcript id mappings.

    Annotation can be stringtie-generated (GTF) or from the input
    reference annotation (GTF or GFF3) and the various attributes can differ
    """
    with open(gtf) as fh:
        txid_to_gene_name = {}
        gid_to_gene_name = {}
        tx_id_to_gene_id = {}

        def get_feature(row, feature):
            return row.split(feature)[1].split(
                ";")[0].replace('=', '').replace("\"", "").strip()

        for gff_entry in fh:
            # Process transcripts features only
            if gff_entry.startswith("#") or gff_entry.split('\t')[2] != 'transcript':
                continue
            # Different gtf/gff formats contain different attributes
            # and different formating (eg. gene_name="xyz" or gene_name "xyz")
            gene_name = gene_id = transcript_id = 'unknown'

            if 'ref_gene_id' in gff_entry:
                # Favour ref_gene_id over gene_id. The latter can be multi-locus merged
                # genes from stringtie
                gene_id = get_feature(gff_entry, 'ref_gene_id')
            elif 'gene_id' in gff_entry:
                gene_id = get_feature(gff_entry, 'gene_id')
            else:
                gene_id = get_feature(gff_entry, 'gene')

            if 'transcript_id' in gff_entry:
                transcript_id = get_feature(gff_entry, 'transcript_id')

            if 'gene_name' in gff_entry:
                gene_name = get_feature(gff_entry, 'gene_name')
            else:
                # Fallback to gene_id if gene_name is not present
                gene_name = gene_id

            txid_to_gene_name[transcript_id] = gene_name
            tx_id_to_gene_id[transcript_id] = gene_id
            gid_to_gene_name[gene_id] = gene_name
    return txid_to_gene_name, tx_id_to_gene_id, gid_to_gene_name


def de_section(
        annotation, dge, dexseq, dtu,
        tpm, report, filtered, unfiltered,
        gene_counts, aln_stats_dir, pval_threshold=0.01):
    """Differential expression sections."""
    with (report.add_section("Differential expression", "DE")):

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

        # Get translations for adding gene names to tables
        (
            txid_to_gene_name,  txid_to_gene_id, gid_to_gene_name
         ) = get_translations(annotation)

        # Add gene names columns to counts files and write out
        # for publishing to user dir.
        df_dge = pd.read_csv(dge, sep='\t')
        df_dge.insert(0, 'gene_name', df_dge.index.map(
            lambda x: gid_to_gene_name.get(x)))
        df_dge.to_csv('results_dge.tsv', index=True, index_label="gene_id", sep="\t")

        # write_dge(gene_counts, gid_to_gene_name, "all_gene_counts.tsv")
        df_gene_counts = pd.read_csv(gene_counts, sep='\t')
        df_gene_counts.insert(
            0, 'gene_name', df_gene_counts.index.map(
                lambda x: gid_to_gene_name.get(x)))
        df_gene_counts.to_csv(
            'results_dge.tsv', index=True, index_label="gene_id", sep="\t")

        df_filtered = pd.read_csv(filtered, sep='\t')
        df_filtered.insert(1, "gene_name", df_filtered.gene_id.map(
            lambda x: gid_to_gene_name.get(x)))
        df_filtered.to_csv(
            'filtered_transcript_counts_with_genes.tsv', index=False, sep='\t')

        df_unfiltered = pd.read_csv(unfiltered, sep='\t')
        df_unfiltered.insert(1, "gene_name", df_unfiltered.gene_id.map(
            lambda x: gid_to_gene_name.get(x)))
        df_unfiltered.to_csv(
            'unfiltered_transcript_counts_with_genes.tsv', index=False, sep='\t')

        df_tpm = pd.read_csv(tpm, sep='\t')
        df_tpm.insert(1, "gene_name", df_tpm.Reference.map(
            lambda x: txid_to_gene_name.get(x)))
        df_tpm.to_csv("unfiltered_tpm_transcript_counts.tsv", index=False, sep='\t')

        # Add tables to report
        dge_section(df_dge, pval_threshold)
        dexseq_section(dexseq, txid_to_gene_name, txid_to_gene_id, pval_threshold)
        dtu_section(dtu, txid_to_gene_name)
