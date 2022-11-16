process count_transcripts {
    // Count transcripts using Salmon.
    // library type is specified as forward stranded (-l SF) as it should have either been through pychopper or come from direct RNA reads.
    label "isoforms"
    cpus params.threads
    input:
        tuple val(sample_id), path(bam)
        path ref_transcriptome
    output:
        path "*transcript_counts.tsv", emit: counts
        path "*seqkit.stats", emit: seqkit_stats
    """
    salmon quant --noErrorModel -p "${task.cpus}" -t "${ref_transcriptome}" -l SF -a "${bam}" -o counts
    mv counts/quant.sf "${sample_id}.transcript_counts.tsv"
    seqkit bam  "${bam}" 2>  "${sample_id}.seqkit.stats"
    """
}


process mergeCounts {
    label "isoforms"
    input:
        path counts
    output:
        path "all_counts.tsv"
    """
    merge_count_tsvs.py -z -o all_counts.tsv -tsvs ${counts}
    """
}

process mergeTPM {
    label "isoforms"
    input:
        path counts
    output:
        path "tpm_counts.tsv"
    """
    merge_count_tsvs.py -o tpm_counts.tsv -z -tpm True -tsvs $counts 
    """
}


process deAnalysis {
    label "isoforms"
    input:
        path condition_sheet
        path merged_tsv 
        path annotation
    output:
        path "de_analysis/results_dtu_stageR.tsv", emit: stageR
        path "merged/all_counts_filtered.tsv", emit: flt_counts
        path "merged/all_gene_counts.tsv", emit: gene_counts
        path "de_analysis/results_dge.tsv", emit: dge
        path  "de_analysis/results_dexseq.tsv", emit: dexseq
        path "de_analysis", emit: de_analysis
    
    """
    cp $annotation annotation.gtf
    echo \$(realpath annotation.gtf)
    echo Annotation\$'\t'min_samps_gene_expr\$'\t'min_samps_feature_expr\$'\t'min_gene_expr\$'\t'min_feature_expr > params.tsv
    echo \$(realpath $params.ref_annotation)\$'\t'$params.min_samps_gene_expr\$'\t'\
    $params.min_samps_feature_expr\$'\t'$params.min_gene_expr\$'\t'$params.min_feature_expr >> params.tsv
    mkdir merged
    mkdir de_analysis
    mv $merged_tsv merged/all_counts.tsv
    mv params.tsv de_analysis/de_params.tsv
    mv $condition_sheet de_analysis/coldata.tsv
    de_analysis.R
    """
}


process plotResults {
    label "isoforms"
    input:
        path flt_count
        path res_dtu
        path condition_sheet 
        path de_analysis
    output:
        path "de_analysis/dtu_plots.pdf", emit: dtu_plots
        path "condition_sheet.tsv", emit: condition_sheet_tsv
        path "de_analysis", emit: stageR
    """
    mkdir merged
    mv $condition_sheet de_analysis/coldata.tsv
    mv $flt_count merged/all_counts_filtered.tsv
    plot_dtu_results.R
    mv de_analysis/coldata.tsv condition_sheet.tsv
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
        path "genome_index.mmi", emit: index
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
       tuple val(sample_id), path (fastq_reads)
       file index
       file transcript_reference
    output:
       tuple val(sample_id), path("${sample_id}_reads_aln_sorted.bam"), emit: bam
    """
    minimap2 -t ${task.cpus} -ax splice -uf -p 1.0 "${index}" "${fastq_reads}" \
    | samtools view -Sb > "output.bam"
    samtools sort -@ ${task.cpus} "output.bam" -o "${sample_id}_reads_aln_sorted.bam"
    samtools index "${sample_id}_reads_aln_sorted.bam"
    """
}


workflow differential_expression {
    take:
       ref_transcriptome
       full_len_reads
       condition_sheet
       ref_annotation
    main:
        t_index = build_minimap_index_transcriptome(ref_transcriptome)
        mapped = map_transcriptome(full_len_reads, t_index, ref_transcriptome)
        count_transcripts(mapped.bam, ref_transcriptome)
        merged = mergeCounts(count_transcripts.out.counts.collect())
        merged_TPM = mergeTPM(count_transcripts.out.counts.collect())
        analysis = deAnalysis(condition_sheet, merged, ref_annotation)
        plotResults(analysis.flt_counts, analysis.stageR, condition_sheet, analysis.de_analysis)
        de_report = analysis.flt_counts.combine(analysis.gene_counts).combine(analysis.dge).combine(analysis.dexseq).combine(
                    analysis.stageR).combine(plotResults.out.condition_sheet_tsv).combine(merged).combine(
                    ref_annotation).combine(merged_TPM)
        count_transcripts_file = count_transcripts.out.seqkit_stats.collect()
emit:
       all_de = de_report
       count_transcripts = count_transcripts_file
       dtu_plots = plotResults.out.dtu_plots
       de_outputs = plotResults.out.stageR
}
