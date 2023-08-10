#!/usr/bin/env python
"""Create de report section."""

from glob import glob
import os

from aplanat import hist, points
from aplanat.bars import boxplot_series
from aplanat.util import Colors
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


def dtu_table(gene_id, dtu_file, alignment_stats, condition_sheet):
    """Create DTU table and plot."""
    dtu_results = pd.read_csv(dtu_file, sep='\t')
    f = open(condition_sheet)
    df = pd.read_csv(f, sep='\t')
    treated_df = df.loc[df['condition'] == "treated"]
    untreated_df = df.loc[df['condition'] == "untreated"]
    treated_samples = treated_df['sample'].tolist()
    control_samples = untreated_df['sample'].tolist()
    # parmaterise these
    control_name = "condition1"
    treated_name = "condition2"
    table = dtu_results.loc[(dtu_results["geneID"] == gene_id)]
    msg = "Gene ID \"{}\" does not exist in the dataset, please select another"
    assert not table.empty, msg.format(gene_id)
    alignment_stats[alignment_stats["Ref"].isin([gene_id])]
    alignment_stats["Transcript"] = alignment_stats["Ref"].apply(
        lambda x: x.split(".")[0])
    gene_alignments = alignment_stats[alignment_stats["Transcript"].isin(
        table["txID"])]
    gene_alignments = gene_alignments.loc[(
        gene_alignments["Type"] == "Primary")]
    gene_alignments["condition"] = gene_alignments.apply(
        lambda x: control_name if x["fname"] in
        control_samples else treated_name, axis=1)
    groups = gene_alignments.groupby(
        ["Transcript", "fname", "condition"]).agg(
            **{"Transcript_count": ("Read", "size")})
    df = None
    for gene_id, group in gene_alignments.groupby("Transcript"):
        temp = group.groupby("fname").agg(**{
            gene_id: ("Read", "size")
        }).transpose()
    if df is None:
        df = temp
    else:
        df = pd.concat([df, temp])
    df = df.reset_index().rename(columns={"index": "Transcript_ID"})
    table = table.rename(columns={
        "txID": "Transcript_ID",
        "gene": "p_gene",
        "transcript": "p_transcript"
        })
    gene_table = pd.merge(df, table)
    gene_table = gene_table.set_index(["geneID", "Transcript_ID"])
    file_names = set(control_samples.keys())
    file_names.update(treated_samples.keys())
    gene_table = gene_table.fillna(0)
    table = groups.reset_index()[
        ["condition", "Transcript", "Transcript_count"]]
    table['transcript, condition'] = \
        table['Transcript'].astype(str) + ', ' + table['condition'].astype(str)
    repeats = groups.groupby(level=['Transcript', 'condition']).size()
    min_rep, max_rep = min(repeats), max(repeats)
    plot = boxplot_series(
        table['transcript, condition'], table['Transcript_count'],
        x_axis_label='Transcript, condition',
        y_axis_label='Transcript count',
        height=200, width=200,
        title="Transcript counts (from {}-{} replicates)".format(
            min_rep, max_rep))
    plot.xaxis.major_label_orientation = 3.1452/2
    if min_rep < 7:
        for renderer in plot.renderers:
            renderer.glyph.line_alpha = 0.2
            try:
                renderer.glyph.fill_alpha = 0.2
            except Exception:
                pass
        plot.circle(
            table['transcript, condition'], table['Transcript_count'],
            fill_color='black', line_color='black')
    return (table, plot)


def pool_csvs(folder):
    """Concat seqkit stats."""
    files = glob(folder + "/*.seqkit.stats")
    dfs = [parse_seqkit(f) for f in files]
    return pd.concat(dfs)


def abundance_histogram(filtered_counts, gene_counts, section):
    """Create plot for abundance of transcripts across all samples."""
    section.markdown("""
Histogram showing the abundance of transcript
counts for genes identified in the analysis.
    """)
    filtered_count_file = filtered_counts
    gene_count_file = gene_counts
    transcripts_per_gene = pd.read_csv(
        filtered_count_file, sep='\t',
        usecols=['gene_id', 'feature_id']).groupby(['gene_id']).agg(['count'])
    transcripts_per_gene.columns = transcripts_per_gene.columns.droplevel()
    gene_ids = pd.read_csv(
        gene_count_file, sep='\t',
        usecols=[0]).index.values.tolist()
    singletons = [
        gene_id for gene_id in gene_ids if gene_id not in transcripts_per_gene.index.values.tolist()] # noqa
    singletons = pd.DataFrame(index=singletons, columns=['count']).fillna(1)
    transcripts_per_gene = pd.concat([singletons, transcripts_per_gene])
    transcript_plot = hist.histogram(
        [transcripts_per_gene['count'].tolist()],
        binwidth=1, colors=[Colors.cerulean])
    transcript_plot.xaxis.axis_label = "Number of isoforms per gene (n)"
    transcript_plot.yaxis.axis_label = "Number of occurences"
    section.markdown("### Transcripts per gene")
    section.plot(transcript_plot)


def dexseq_section(dexseq_file, section, id_dic):
    """Add gene isoforms table and plot."""
    section.markdown("### Differential Isoform usage")
    dexseq_caption = '''Table showing gene isoforms, ranked by adjusted
    p-value, from the DEXSeq analysis. Information shown includes the log2 fold
    change between experimental conditions, the log-scaled transcript
    abundance and the false discovery corrected p-value (FDR).
    This table has not been filtered
    for genes that satisfy statistical or magnitudinal thresholds'''
    section.markdown(dexseq_caption)
    dexseq_results = pd.read_csv(dexseq_file, sep='\t')
    dexseq_results.index.name = "gene_id:trancript_id"
    # Replace gene id with more useful gene name where possible
    dexseq_results.index = dexseq_results.index.map(
        lambda x: str(id_dic.get(x.split(':')[0])) + ':' + str(x.split(':')[1]))
    dexseq_pvals = dexseq_results.sort_values(by='pvalue', ascending=True)
    section.table(dexseq_results.loc[dexseq_pvals.index], index=True)
    section.markdown("""
The figure below presents the MA plot from the DEXSeq analysis.
M is the log2 ratio of isoform transcript abundance between conditions.
A is the log2 transformed mean abundance value.
Transcripts that satisfy the logFC and FDR corrected p-value
thresholds defined are shaded as 'Up-' or 'Down-' regulated.""")
    pval_limit = 0.01
    up = dexseq_results.loc[
        (dexseq_results["Log2FC"] > 0) & (
            dexseq_results['pvalue'] < pval_limit)]
    down = dexseq_results.loc[
        (dexseq_results["Log2FC"] <= 0) & (
            dexseq_results['pvalue'] < pval_limit)]
    not_sig = dexseq_results.loc[(dexseq_results["pvalue"] >= pval_limit)]

    dexseq_plot = points.points(
        x_datas=[
                up["Log2MeanExon"],
                down["Log2MeanExon"],
                not_sig["Log2MeanExon"],
                ],
        y_datas=[
                up["Log2FC"],
                down["Log2FC"],
                not_sig["Log2FC"],
                ],
        title="Average copy per million (CPM) vs Log-fold change (LFC)",
        colors=["red", "blue", "black"],
        names=["Up", "Down", "NotSig"]
        )

    dexseq_plot.xaxis.axis_label = "A (log2 transformed mean exon read counts)"
    dexseq_plot.yaxis.axis_label = """
    M (log2 transformed differential abundance)
    """
    dexseq_results_caption = "### Dexseq results"
    section.markdown(dexseq_results_caption)
    section.plot(dexseq_plot)


def dtu_section(dtu_file, section, gt_dic, ge_dic):
    """Plot dtu section."""
    dtu_results = pd.read_csv(dtu_file, sep='\t')
    dtu_results["gene_name"] = dtu_results["txID"].apply(
        lambda x: gt_dic.get(x))
    dtu_results["geneID"] = dtu_results["geneID"].apply(
        lambda x: ge_dic.get(x))
    dtu_pvals = dtu_results.sort_values(by='gene', ascending=True)
    dtu_caption = '''Table showing gene and transcript identifiers
    and their FDR corrected probabilities
    for the genes and their isoforms that have been
    identified as showing DTU using the R packages DEXSeq and StageR.
    This list has been shortened requiring that both gene and transcript
    must satisfy the p-value
    threshold'''
    section.markdown(dtu_caption)
    section.table(dtu_results.loc[dtu_pvals.index])


def dge_names(dge_file, geid_gname):
    """Add gene name column to DGE tsv."""
    dge_results = pd.read_csv(dge_file, sep='\t')
    dge_results["gene_name"] = dge_results.index.map(lambda x: geid_gname.get(x))
    dge_results.to_csv('results_dge.tsv', index=True, index_label="gene_id")


def dge_section(dge_file, section, ids_dic):
    """Create DGE table and plot."""
    section.markdown('### Differential gene expression')
    dge_results = pd.read_csv(dge_file, sep='\t')
    dge_pvals = dge_results.sort_values(by='FDR', ascending=True)
    dge_results[['logFC', 'logCPM', 'F']] = dge_results[
        ['logFC', 'logCPM', 'F']].round(2)
    dge_caption = """
Table showing the genes from the edgeR analysis.
Information shown includes the log2 fold change between
experimental conditions, the log-scaled counts per million measure of abundance
and the false discovery corrected p-value (FDR). This table has not been
filtered for genes that satisfy statistical or magnitudinal thresholds"""
    section.markdown(dge_caption)
    dge_results.index = dge_results.index.map(lambda x: ids_dic.get(x))
    dge_pvals.index = dge_pvals.index.map(lambda x: ids_dic.get(x))
    section.table(dge_results.loc[dge_pvals.index], index=True)
    dge = pd.read_csv(dge_file, sep="\t")
    section.markdown("""
This plot visualises differences in measurements between the
two experimental conditions. M is the log2 ratio of gene expression
calculated between the conditions.
A is a log2 transformed mean expression value.
The figure below presents the MA figure from this edgeR analysis.
Genes that satisfy the logFC and FDR corrected p-value thresholds
defined are shaded as 'Up-' or 'Down-' regulated.
    """)
    pval_limit = 0.01
    up = dge.loc[(dge["logFC"] > 0) & (dge['PValue'] < pval_limit)]
    down = dge.loc[(dge["logFC"] <= 0) & (dge['PValue'] < pval_limit)]
    not_sig = dge.loc[(dge["PValue"] >= pval_limit)]
    logcpm_vs_logfc = points.points(
        x_datas=[
                up["logCPM"],
                down["logCPM"],
                not_sig["logCPM"],
                ],
        y_datas=[
                up["logFC"],
                down["logFC"],
                not_sig["logFC"],
                ],
        title="Average copy per million (CPM) vs Log-fold change (LFC)",
        colors=["red", "blue", "black"],
        names=["Up", "Down", "NotSig"]
        )
    logcpm_vs_logfc.xaxis.axis_label = "Average log CPM"
    logcpm_vs_logfc.yaxis.axis_label = "Log-fold change"
    logcpm_caption = """### Results of the edgeR Analysis."""
    section.markdown(logcpm_caption)
    section.plot(logcpm_vs_logfc)


def salmon_table(salmon_counts, section):
    """Create salmon counts summary table."""
    salmon_counts = pd.read_csv(salmon_counts, sep='\t')
    salmon_counts.set_index("Reference", drop=True, append=False, inplace=True)
    salmon_size_top = salmon_counts.sum(axis=1).sort_values(ascending=False)
    salmon_counts = salmon_counts.applymap(np.int64)
    salmon_count_caption = """
    Table showing the annotated Transcripts Per Million
    identified by Minimap2 mapping and Salmon transcript
    detection with the highest
    number of mapped reads"""
    section.markdown("### Transcripts Per Million ")
    section.markdown(salmon_count_caption, "salmon-head-caption")
    section.table(
        salmon_counts.loc[salmon_size_top.index].head(n=100), index=True)


def get_translations(gtf):
    """Create dict with gene_name and gene_references."""
    fn = open(gtf).readlines()
    gene_txid = {}
    gene_geid = {}
    geid_gname = {}

    def get_feature(row, feature):
        return row.split(feature)[1].split(
            ";")[0].replace('=', '').replace("\"", "").strip()

    for i in fn:
        if i.startswith("#"):
            continue
        # Different gtf/gff formats contain different attributes
        # and different formating (eg. gene_name="xyz" or gene_name "xyz")
        if 'gene_name' in i:
            gene_name = get_feature(i, "gene_name")
        elif 'gene_id' in i:
            gene_name = get_feature(i, 'gene_id')
        elif 'gene' in i:
            gene_name = get_feature(i, "gene")
        else:
            continue

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
        tpm, report):
    """Differential expression sections."""
    section = report.add_section()
    section.markdown("# Differential expression.")
    section.markdown("""
This section shows differential gene expression
and differential isoform usage. Salmon was used to
assign reads to individual annotated isoforms defined by
the GTF-format annotation.
These counts were used to perform a statistical analysis to identify
the genes and isoforms that show differences in abundance between
the experimental conditions.
Any novel genes or transcripts that do not have relevant gene or transcript IDs
are prefixed with MSTRG for use in differential expression analysis.
Find the full sequences of any transcripts in the
`final_non_redundant_transcriptome.fasta` file.
    """)
    section.markdown("### Alignment summary stats")
    alignment_stats = pool_csvs("seqkit")
    alignment_summary_df = create_summary_table(alignment_stats)
    alignment_summary_df = alignment_summary_df.fillna(0).applymap(np.int64)
    section.table(alignment_summary_df, key='alignment-stats', index=True)
    salmon_table(tpm, section)
    gene_txid, gene_name, geid_gname = get_translations(stringtie)
    dge_section(dge, section, gene_name)
    dge_names(dge, geid_gname)
    dexseq_section(dexseq, section, gene_name)
    dtu_section(dtu, section, gene_txid, gene_name)
    # missing dtu plots at the moment as too many
    section.markdown("""
### View dtu_plots.pdf file to see plots of differential isoform usage
""")
