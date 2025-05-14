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
    memory "31 GB"
    input:
        tuple val(meta), path(bam), path(ref_transcriptome)
    output:
        path "*transcript_counts.tsv", emit: counts
    """
    salmon quant --noErrorModel -p "${task.cpus}" -t "${ref_transcriptome}" -l SF -a "${bam}" -o counts
    mv counts/quant.sf "${meta.alias}.transcript_counts.tsv"
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
    cpus 4
    memory "16 GB"
    input:
        path "sample_sheet.csv"
        path "all_counts.tsv" 
        path "annotation.gtf"
    output:
        path "de_analysis/results_dtu_stageR.tsv", emit: stageR
        path "merged/filtered_transcript_counts_with_genes.tsv", emit: flt_counts
        path "merged/all_gene_counts.tsv", emit: gene_counts
        path "de_analysis/unfiltered_transcript_counts_with_genes.tsv", emit: unflt_counts
        path "de_analysis/results_dge.tsv", emit: dge
        path "de_analysis/results_dexseq.tsv", emit: dexseq
        path "de_analysis/results_dge.pdf", emit: dge_pdf
        path "de_analysis/results_dge.tsv", emit: dge_tsv
        path "de_analysis/results_dtu_gene.tsv", emit: dtu_gene
        path "de_analysis/results_dtu_transcript.tsv", emit: dtu_transcript
        path "de_analysis/results_dtu_stageR.tsv", emit: dtu_stageR
        path "de_analysis/results_dtu.pdf", emit: dtu_pdf
        path "de_analysis/cpm_gene_counts.tsv", emit: cpm
    """
    de_analysis.R \
        --annotation annotation.gtf \
        --min_samps_gene_expr $params.min_samps_gene_expr \
        --min_samps_feature_expr $params.min_samps_feature_expr \
        --min_gene_expr $params.min_gene_expr \
        --min_feature_expr $params.min_feature_expr \
        --sample_sheet sample_sheet.csv \
        --all_counts all_counts.tsv \
        --de_out_dir de_analysis \
        --merged_out_dir merged

    # Check that the original aliases in the input TSV have not been mangled by R's read.csv or other functions
    head -1 all_counts.tsv | cut -f2- | tr '\t' '\n'   > expected_colnames
    
    head -1 de_analysis/cpm_gene_counts.tsv | cut -f2- | tr '\t' '\n'  > cpm_gene_counts_colnames
    head -1 merged/all_gene_counts.tsv | tr '\t' '\n'  > merged_counts_colnames
    head -1 merged/filtered_transcript_counts_with_genes.tsv | cut -f3- | tr '\t' '\n' > merged_filtered_colnames
    
    # Check for mismatches in sample column names
    for file in cpm_gene_counts_colnames merged_counts_colnames merged_filtered_colnames; do
    if ! diff -q \$file expected_colnames > /dev/null; then
        echo "Column names in \$file do not match expected aliases."
        exit 70
    fi
    done
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
    output:
        path "dtu_plots.pdf", emit: dtu_plots
    """
    plot_dtu_results.R \
    --counts filtered_transcript_counts_with_genes.tsv \
    --results_dtu results_dtu_stageR.tsv \
    --sample_sheet sample_sheet.tsv \
    --pdf_out dtu_plots.pdf
    """
}

process build_minimap_index_transcriptome{
    /*
    Build minimap index from reference genome
    */
    label "isoforms"
    cpus params.threads
    memory "31 GB"
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
       path("${meta.alias}.flagstat.stats"), emit: align_stats
    """
    minimap2 -t ${task.cpus} -ax splice -uf -p 1.0 "${index}" "${fastq_reads}" \
    | samtools view -Sb > "output.bam"
    samtools sort -@ ${task.cpus} "output.bam" -o "${meta.alias}_reads_aln_sorted.bam"
    samtools flagstat -O json "${meta.alias}_reads_aln_sorted.bam" > "${meta.alias}.flagstat.stats"
    """
}


workflow differential_expression {
    take:
       ref_transcriptome
       full_len_reads
       sample_sheet
       ref_annotation
    main:
        sample_sheet = Channel.fromPath(sample_sheet)
        checkSampleSheetCondition(sample_sheet)
        t_index = build_minimap_index_transcriptome(ref_transcriptome)
        mapped = map_transcriptome(full_len_reads.combine(t_index)
        .map{meta, fastq, reference, transcriptome -> tuple(meta, fastq, reference) })
        count_transcripts(mapped.bam.combine(t_index.map{ mmi, reference -> reference}))
        merged = mergeCounts(count_transcripts.out.counts.collect())
        merged_TPM = mergeTPM(count_transcripts.out.counts.collect())
        analysis = deAnalysis(sample_sheet, merged, ref_annotation)
        plotResults(analysis.flt_counts, analysis.stageR, sample_sheet)
        // Concat files required for making the report
        de_report = analysis.flt_counts.concat(
            analysis.gene_counts, analysis.dge, analysis.dexseq,
            analysis.stageR, sample_sheet, merged, ref_annotation, merged_TPM, analysis.unflt_counts).collect()
        // Concat files required to be output to user without any changes
        de_outputs_concat = analysis.cpm.concat(plotResults.out.dtu_plots, analysis.dge_pdf, analysis.dge_tsv,
        analysis.dtu_gene, analysis.dtu_transcript, analysis.dtu_stageR, analysis.dtu_pdf, merged_TPM).collect()
        collected_de_alignment_stats = mapped.align_stats.collect()
emit:
       all_de = de_report
       de_alignment_stats = collected_de_alignment_stats
       de_outputs = de_outputs_concat
}
