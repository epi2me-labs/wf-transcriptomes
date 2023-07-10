process checkSampleSheetCondition {
    label "isoforms"
    input:
        path "sample_sheet.csv"
    """
    workflow-glue check_sample_sheet_condition "sample_sheet.csv"
    """
}



process count_transcripts {
    // Count transcripts using Salmon.
    // library type is specified as forward stranded (-l SF) as it should have either been through pychopper or come from direct RNA reads.
    label "isoforms"
    cpus params.threads
    input:
        tuple val(meta), path(bam), path(ref_transcriptome)
    output:
        path "*transcript_counts.tsv", emit: counts
        path "*seqkit.stats", emit: seqkit_stats
    """
    salmon quant --noErrorModel -p "${task.cpus}" -t "${ref_transcriptome}" -l SF -a "${bam}" -o counts
    mv counts/quant.sf "${meta.alias}.transcript_counts.tsv"
    seqkit bam  "${bam}" 2>  "${meta.alias}.seqkit.stats"
    """
}


process mergeCounts {
    label "isoforms"
    input:
        path counts
    output:
        path "all_counts.tsv"
    """
    workflow-glue merge_count_tsvs -z -o all_counts.tsv -tsvs ${counts}
    """
}

process mergeTPM {
    label "isoforms"
    input:
        path counts
    output:
        path "tpm_counts.tsv"
    """
    workflow-glue merge_count_tsvs -o tpm_counts.tsv -z -tpm True -tsvs $counts
    """
}


process deAnalysis {
    label "isoforms"
    errorStrategy "retry"
    maxRetries 1
    input:
        path sample_sheet
        path merged_tsv 
        path "annotation.gtf"
    output:
        path "de_analysis/results_dtu_stageR.tsv", emit: stageR
        path "merged/all_counts_filtered.tsv", emit: flt_counts
        path "merged/all_gene_counts.tsv", emit: gene_counts
        path "de_analysis/results_dge.tsv", emit: dge
        path "de_analysis/results_dexseq.tsv", emit: dexseq
        path "de_analysis", emit: de_analysis
    script:
    // Just try both annotation file type because a .gff extension may be gff2(gtf) or gff3
    String annotation_type = "gtf"
    if (task.attempt == 2){
        annotation_type = "gff3"
    }
    """
    mkdir merged
    mkdir de_analysis
    mv $merged_tsv merged/all_counts.tsv
    mv $sample_sheet de_analysis/coldata.tsv
    de_analysis.R annotation.gtf $params.min_samps_gene_expr $params.min_samps_feature_expr $params.min_gene_expr $params.min_feature_expr $annotation_type
   
    """
}


process plotResults {
    label "isoforms"
    input:
        path flt_count
        path res_dtu
        path sample_sheet 
        path de_analysis
    output:
        path "de_analysis/dtu_plots.pdf", emit: dtu_plots
        path "sample_sheet.tsv", emit: sample_sheet_csv
        path "de_analysis", emit: stageR
    """
    mkdir merged
    mv $sample_sheet de_analysis/coldata.tsv
    mv $flt_count merged/all_counts_filtered.tsv
    plot_dtu_results.R
    mv de_analysis/coldata.tsv sample_sheet.tsv
    """
}

process build_minimap_index_transcriptome{
    /*
    Build minimap index from reference genome
    */
    label "isoforms"
    cpus params.threads
    input:
        path reference
    output:
        tuple path("genome_index.mmi"), path(reference), emit: index
    script:
    """
    minimap2 -t "${task.cpus}" ${params.minimap_index_opts}  -I 1000G -d "genome_index.mmi" "${reference}"
  
    """
}


process map_transcriptome{
    /*
    Map reads to reference using minimap2.
    Filter reads by mapping quality.
    Filter internally-primed reads.
    */
    label "isoforms"
    cpus params.threads

    input:
       tuple val(meta), path (fastq_reads), path(index)
    output:
       tuple val(meta), path("${meta.alias}_reads_aln_sorted.bam"), emit: bam
    """
    minimap2 -t ${task.cpus} -ax splice -uf -p 1.0 "${index}" "${fastq_reads}" \
    | samtools view -Sb > "output.bam"
    samtools sort -@ ${task.cpus} "output.bam" -o "${meta.alias}_reads_aln_sorted.bam"
    """
}


workflow differential_expression {
    take:
       ref_transcriptome
       full_len_reads
       sample_sheet
       ref_annotation
    main:
        checkSampleSheetCondition(sample_sheet)
        t_index = build_minimap_index_transcriptome(ref_transcriptome)
        mapped = map_transcriptome(full_len_reads.combine(t_index))
        count_transcripts(mapped.bam.combine(t_index.map{ mmi, reference -> reference}))
        merged = mergeCounts(count_transcripts.out.counts.collect())
        merged_TPM = mergeTPM(count_transcripts.out.counts.collect())
        analysis = deAnalysis(sample_sheet, merged, ref_annotation)
        plotResults(analysis.flt_counts, analysis.stageR, sample_sheet, analysis.de_analysis)
        de_report = analysis.flt_counts.combine(analysis.gene_counts).combine(analysis.dge).combine(analysis.dexseq).combine(
                    analysis.stageR).combine(plotResults.out.sample_sheet_csv).combine(merged).combine(
                    ref_annotation).combine(merged_TPM)
        count_transcripts_file = count_transcripts.out.seqkit_stats.collect()
emit:
       all_de = de_report
       count_transcripts = count_transcripts_file
       dtu_plots = plotResults.out.dtu_plots
       de_outputs = plotResults.out.stageR
}
