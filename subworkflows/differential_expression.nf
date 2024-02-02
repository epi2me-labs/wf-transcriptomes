process checkSampleSheetCondition {
    label "isoforms"
    cpus 1
    memory "2 GB"
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
    memory "16 GB"
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
    cpus 1
    memory "2 GB"
    input:
        path counts
    output:
        path "unfiltered_transcript_counts.tsv"
    """
    workflow-glue merge_count_tsvs -z -o unfiltered_transcript_counts.tsv -tsvs ${counts}
    """
}

process mergeTPM {
    label "isoforms"
    cpus 1
    memory "2 GB"
    input:
        path counts
    output:
        path "unfiltered_tpm_transcript_counts.tsv"
    // Use tpm parameter with merge_counts_tsvs.py to out transcript per million file
    """
    workflow-glue merge_count_tsvs -o unfiltered_tpm_transcript_counts.tsv -z -tpm True -tsvs $counts
    """
}


process deAnalysis {
    label "isoforms"
    errorStrategy "retry"
    maxRetries 3
    cpus 4
    memory "16 GB"
    input:
        path "sample_sheet.csv"
        path "all_counts.tsv" 
        path "annotation.gtf"
    output:
        path "de_analysis/results_dtu_stageR.tsv", emit: stageR
        path "merged/filtered_transcript_counts_with_genes.tsv", emit: flt_counts
        path "de_analysis/unfiltered_transcript_counts_with_genes.tsv", emit: unflt_counts
        path "merged/all_gene_counts.tsv", emit: gene_counts
        path "de_analysis/results_dge.tsv", emit: dge
        path "de_analysis/results_dexseq.tsv", emit: dexseq
        path "de_analysis", emit: de_analysis
    script:
    // Just try both annotation file type because a .gff extension may be gff2(gtf) or gff3
    String annotation_type = "gtf"
    String strip_version  = "false"
    if (task.attempt == 2){
        annotation_type = "gff3"
        strip_version  = "false"
        log.info("Retry deAnalysis with gff format setting.")
    }
    else if (task.attempt == 3){
        annotation_type = "gff3"
        strip_version  = "true"
        log.info("Retry deAnalysis with gff format setting and version removal.")
    }
    else if (task.attempt == 4){
        strip_version  = "true"
        log.info("Retry deAnalysis with gtf format setting and version removal.")
    }

    """
    mkdir merged
    mkdir de_analysis
    de_analysis.R annotation.gtf $params.min_samps_gene_expr $params.min_samps_feature_expr $params.min_gene_expr $params.min_feature_expr $annotation_type $strip_version
    """
}


process plotResults {
    label "isoforms"
    cpus 2
    memory "2 GB"
    input:
        path "filtered_transcript_counts_with_genes.tsv"
        path "results_dtu_stageR.tsv"
        path "sample_sheet.tsv"
        path de_analysis
    output:
        path "de_analysis/dtu_plots.pdf", emit: dtu_plots
        path "sample_sheet.tsv", emit: sample_sheet_csv
        path "de_analysis/*", emit: stageR
    """
    plot_dtu_results.R
    """
}

process build_minimap_index_transcriptome{
    /*
    Build minimap index from reference genome
    */
    label "isoforms"
    cpus params.threads
    memory "16 GB"
    input:
        path reference
    output:
        tuple path("genome_index.mmi"), path(reference), emit: index
    script:
    """
    minimap2 -t "${task.cpus}" ${params.minimap2_index_opts}  -I 1000G -d "genome_index.mmi" "${reference}"
  
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
    memory "16 GB"

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
        mapped = map_transcriptome(full_len_reads.combine(t_index)
        .map{meta, fastq, reference, transcriptome -> tuple(meta, fastq, reference) })
        count_transcripts(mapped.bam.combine(t_index.map{ mmi, reference -> reference}))
        merged = mergeCounts(count_transcripts.out.counts.collect())
        merged_TPM = mergeTPM(count_transcripts.out.counts.collect())
        analysis = deAnalysis(sample_sheet, merged, ref_annotation)
        plotResults(analysis.flt_counts, analysis.stageR, sample_sheet, analysis.de_analysis)
        de_report = analysis.flt_counts.concat(analysis.gene_counts, analysis.dge, analysis.dexseq,
        analysis.stageR, plotResults.out.sample_sheet_csv, merged, ref_annotation, merged_TPM, analysis.unflt_counts).collect()
        count_transcripts_file = count_transcripts.out.seqkit_stats.collect()
        all_counts = merged_TPM.concat(analysis.flt_counts, analysis.gene_counts)
emit:
       all_de = de_report
       count_transcripts = count_transcripts_file
       dtu_plots = plotResults.out.dtu_plots
       de_outputs = plotResults.out.stageR
       counts = all_counts
}
