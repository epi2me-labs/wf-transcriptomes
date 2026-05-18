nextflow.enable.dsl = 2

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")


process bambuDiscover {
    label "wf_transcriptomes"
    cpus {
        int requested = (params.threads ?: 4) as int
        int sampleCount = aliases instanceof Collection ? aliases.size() : 1
        sampleCount > 1 ? requested : 1
    }
    memory "60 GB"
    input:
        tuple val(meta), val(aliases), path(bams, stageAs: "bams/??.bam"), path(bais, stageAs: "bams/??.bam.bai"), path(sample_sheet)
        path annotation, stageAs: "annotation/*"
        path reference, stageAs: "reference/*"
    output:
        tuple val(meta), path("discover"), emit: dir
    script:
        def bam_list = bams instanceof Collection ? bams : [bams]
        def alias_list = aliases instanceof Collection ? aliases : [aliases]
        String bams_arg = "--bams '${bam_list.join(",")}'"
        String aliases_arg = "--aliases '${alias_list.join(",")}'"
        String sample_sheet_arg = sample_sheet.name == OPTIONAL_FILE.name ? "" : "--sample_sheet '${sample_sheet}'"
        String ndr_arg = params.ndr != null ? "--ndr ${params.ndr}" : ""
    """
    supeRglue bambu discover \
        ${bams_arg} \
        ${aliases_arg} \
        ${sample_sheet_arg} \
        --annotation "${annotation}" \
        --genome "${reference}" \
        --transcriptome_mode "${params.transcriptome_mode}" \
        --threads ${task.cpus} \
        ${ndr_arg} \
        --out_dir discover
    """
}


process bambuQuant {
    label "wf_transcriptomes"
    cpus {
        int requested = (params.threads ?: 4) as int
        boolean isJoint = meta instanceof Map && meta.alias == 'cohort'
        isJoint ? requested : 1
    }
    memory { ["8.GB", "16.GB", "48.GB"][task.attempt - 1] }
    maxRetries 2
    errorStrategy 'retry'
    input:
        tuple val(meta), val(chunk_id), val(annotation_tx_count), path(chunk_rds), path(discovered_annotation)
        path reference, stageAs: "reference/*"
    output:
        tuple val(meta), val(chunk_id), path("${chunk_id}"), emit: dir
    script:
    """
    supeRglue bambu quant \
        --chunk_rds "${chunk_rds}" \
        --discovered_annotation_rds "${discovered_annotation}" \
        --genome "${reference}" \
        --threads ${task.cpus} \
        --out_dir "${chunk_id}"
    """
}


process bambuEmpty {
    label "wf_transcriptomes"
    cpus 1
    memory "4 GB"
    input:
        tuple val(meta), val(aliases)
    output:
        tuple val(meta), path("${meta.alias}"), emit: dir
        tuple val(meta), path("${meta.alias}/transcripts.gtf"), emit: gtf
        tuple val(meta), path("${meta.alias}/transcript_counts.tsv"), emit: transcript_counts
        tuple val(meta), path("${meta.alias}/gene_counts.tsv"), emit: gene_counts
        tuple val(meta), path("${meta.alias}/bambu_transcripts.rds"), emit: transcript_rds
        tuple val(meta), path("${meta.alias}/bambu_genes.rds"), emit: gene_rds
        tuple val(meta), path("${meta.alias}/transcript_metadata.tsv"), emit: transcript_metadata
    script:
        def alias_list = aliases instanceof Collection ? aliases : [aliases]
        String aliases_arg = "--aliases '${alias_list.join(",")}'"
    """
    supeRglue bambu empty \
        ${aliases_arg} \
        --transcriptome_mode "${params.transcriptome_mode}" \
        --out_dir "${meta.alias}"
    """
}


process collateBambuQuant {
    label "wf_transcriptomes"
    cpus 1
    memory "16 GB"
    input:
        tuple val(meta), path(chunk_dirs, stageAs: "chunks/*")
    output:
        tuple val(meta), path("${meta.alias}"), emit: dir
        tuple val(meta), path("${meta.alias}/transcripts.gtf"), emit: gtf
        tuple val(meta), path("${meta.alias}/transcript_counts.tsv"), emit: transcript_counts
        tuple val(meta), path("${meta.alias}/gene_counts.tsv"), emit: gene_counts
        tuple val(meta), path("${meta.alias}/bambu_transcripts.rds"), emit: transcript_rds
        tuple val(meta), path("${meta.alias}/bambu_genes.rds"), emit: gene_rds
        tuple val(meta), path("${meta.alias}/transcript_metadata.tsv"), emit: transcript_metadata
    script:
        def chunk_dir_list = chunk_dirs instanceof Collection ? chunk_dirs : [chunk_dirs]
        String chunk_dirs_arg = "--chunk_dirs '${chunk_dir_list.join(",")}'"
        String ndr_arg = params.ndr != null ? "--ndr ${params.ndr}" : ""
    """
    supeRglue bambu collate \
        ${chunk_dirs_arg} \
        --transcriptome_mode "${params.transcriptome_mode}" \
        ${ndr_arg} \
        --out_dir "${meta.alias}"
    """
}
