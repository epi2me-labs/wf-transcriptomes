nextflow.enable.dsl = 2

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")


process prepareAnnotationReference {
    label "wf_transcriptomes"
    cpus 1
    memory "6 GB"
    input:
        path ref_annotation
        tuple path(ref), path(ref_idx)
    output:
        stdout emit: warnings
        path "annotation.gtf", emit: annotation
        path "reference.fasta", emit: reference
        path "annotation_reference_summary.json", emit: summary
        path "unstranded_annotation.gtf", optional: true, emit: unstranded
    script:
    """
    workflow-glue prepare_annotation_reference \
        --annotation "${ref_annotation}" \
        --reference "${ref}" \
        --out_dir prepared
    mv prepared/* .
    """
}


process runJointBambu {
    label "wf_transcriptomes"
    cpus { params.threads ?: 4 }
    memory "32 GB"
    input:
        tuple val(aliases), path(bams, stageAs: "bams/??.bam"), path(bais, stageAs: "bams/??.bam.bai")
        path sample_sheet
        path annotation, stageAs: "annotation/*"
        path reference, stageAs: "reference/*"
    output:
        path "cohort", emit: dir
        path "cohort/transcripts.gtf", emit: gtf
        path "cohort/transcript_counts.tsv", emit: transcript_counts
        path "cohort/gene_counts.tsv", emit: gene_counts
        path "cohort/bambu_transcripts.rds", emit: transcript_rds
        path "cohort/bambu_genes.rds", emit: gene_rds
        path "cohort/transcript_metadata.tsv", emit: transcript_metadata
    script:
        def bam_list = bams instanceof Collection ? bams : [bams]  // todo dont run joint on single sample anyway
        def alias_list = aliases instanceof Collection ? aliases : [aliases]
        String bams_arg = "--bams '${bam_list.join(",")}'"
        String aliases_arg = "--aliases '${alias_list.join(",")}'"
        String sample_sheet_arg = sample_sheet.name == OPTIONAL_FILE.name ? "" : "--sample_sheet ${sample_sheet}"
        String ndr_arg = params.ndr != null ? "--ndr ${params.ndr}" : ""
    """
    supeRglue bambu \
        ${bams_arg} \
        ${aliases_arg} \
        ${sample_sheet_arg} \
        --annotation "${annotation}" \
        --genome "${reference}" \
        --transcriptome_mode "${params.transcriptome_mode}" \
        --threads ${task.cpus} \
        ${ndr_arg} \
        --out_dir cohort \
    """
}


process runPerSampleBambu {
    label "wf_transcriptomes"
    cpus { params.threads ?: 4 }
    memory "24 GB"
    input:
        tuple val(meta), path(bam), path(bai), path(stats)
        path annotation, stageAs: "annotation/*"
        path reference, stageAs: "reference/*"
    output:
        tuple val(meta), path("${meta.alias}"), emit: dir
        tuple val(meta), path("${meta.alias}/transcripts.gtf"), emit: gtf
        tuple val(meta), path("${meta.alias}/transcript_counts.tsv"), emit: transcript_counts
        tuple val(meta), path("${meta.alias}/gene_counts.tsv"), emit: gene_counts
        tuple val(meta), path("${meta.alias}/bambu_transcripts.rds"), emit: transcript_rds
        tuple val(meta), path("${meta.alias}/bambu_genes.rds"), emit: gene_rds
        tuple val(meta), path("${meta.alias}/transcript_metadata.tsv"), emit: transcript_metadata
    script:
        String bams_arg = "--bams '${bam.toString()}'"
        String aliases_arg = "--aliases '${meta.alias}'"
        String ndr_arg = params.ndr != null ? "--ndr ${params.ndr}" : ""
    """
    supeRglue bambu \
        ${bams_arg} \
        ${aliases_arg} \
        --annotation "${annotation}" \
        --genome "${reference}" \
        --transcriptome_mode "${params.transcriptome_mode}" \
        --threads ${task.cpus} \
        ${ndr_arg} \
        --out_dir "${meta.alias}"
    """
}


process buildCohortTranscriptomeFasta {
    label "wf_transcriptomes"
    cpus 1
    memory "4 GB"
    input:
        path "transcripts.gtf"
        path reference
    output:
        path "cohort.transcriptome.fa", emit: fasta
    script:
    """
    gffread -g "${reference}" -w cohort.transcriptome.fa transcripts.gtf
    """
}


process buildSampleTranscriptomeFasta {
    label "wf_transcriptomes"
    cpus 1
    memory "4 GB"
    input:
        tuple val(meta), path("transcripts.gtf")
        path reference
    output:
        tuple val(meta), path("${meta.alias}.transcriptome.fa"), emit: fasta
    script:
    """
    gffread -g "${reference}" -w "${meta.alias}.transcriptome.fa" transcripts.gtf
    """
}


process runJointSqanti {
    label "wf_transcriptomes_sqanti"
    cpus { params.threads ?: 4 }
    memory "24 GB"
    input:
        path gtf
        path annotation, stageAs: "annotation/*"
        path reference, stageAs: "reference/*"
    output:
        path "sqanti_cohort", emit: dir
        path "sqanti_cohort/classification_summary.tsv", emit: summary
    script:
        String extra = params.sqanti_extra_args ?: ""
        String skip_orf = params.sqanti_skip_orf ? "--skipORF" : ""
    """
    mkdir sqanti_cohort
    sqanti3_qc.py \
        --isoforms "${gtf}" \
        --refGTF "${annotation}" \
        --refFasta "${reference}" \
        ${skip_orf} \
        --force_id_ignore \
        --report skip \
        -t ${task.cpus} \
        -d sqanti_cohort \
        -o cohort \
        ${extra}
    workflow-glue summarise_sqanti --sqanti_dir sqanti_cohort \
        --output sqanti_cohort/classification_summary.tsv
    """
}


process runPerSampleSqanti {
    label "wf_transcriptomes_sqanti"
    cpus { params.threads ?: 4 }
    memory "24 GB"
    input:
        tuple val(meta), path(gtf)
        path annotation, stageAs: "annotation/*"
        path reference, stageAs: "reference/*"
    output:
        tuple val(meta), path("${meta.alias}_sqanti"), emit: dir
        tuple val(meta), path("${meta.alias}_sqanti/classification_summary.tsv"), emit: summary
    script:
        String extra = params.sqanti_extra_args ?: ""
        String skip_orf = params.sqanti_skip_orf ? "--skipORF" : ""
    """
    mkdir "${meta.alias}_sqanti"
    sqanti3_qc.py \
        --isoforms "${gtf}" \
        --refGTF "${annotation}" \
        --refFasta "${reference}" \
        ${skip_orf} \
        --force_id_ignore \
        --report skip \
        -t ${task.cpus} \
        -d "${meta.alias}_sqanti" \
        -o "${meta.alias}" \
        ${extra}
    workflow-glue summarise_sqanti --sqanti_dir "${meta.alias}_sqanti" \
        --output "${meta.alias}_sqanti/classification_summary.tsv"
    """
}


workflow transcriptome_analysis {
    take:
        alignments
        ref_genome
        ref_annotation
        sample_sheet
    main:
        prepared_reference_annotation = prepareAnnotationReference(ref_annotation, ref_genome)
        prepared_reference_annotation.warnings.map { stdoutput ->
            if (stdoutput) {
                log.warn(stdoutput.trim())
            }
        }
        analysis_annotation = prepared_reference_annotation.annotation.first()
        analysis_reference = prepared_reference_annotation.reference.first()

        joint_bambu = runJointBambu(
            alignments
            | collect(flat: false)
            | map { rows ->
                // transform [meta, bam, bai] to [[alias1...aliasN], [bam1...bamN], [bai1...baiN]]
                tuple(
                    rows.collect { it[0].alias },
                    rows.collect { it[1] },
                    rows.collect { it[2] }
                )
            },
            sample_sheet,
            analysis_annotation,
            analysis_reference
        )

        sample_bambu = runPerSampleBambu(alignments, analysis_annotation, analysis_reference)

        joint_fasta = buildCohortTranscriptomeFasta(joint_bambu.gtf, analysis_reference)
        sample_fastas = buildSampleTranscriptomeFasta(sample_bambu.gtf, analysis_reference)

        if (params.skip_sqanti) {
            joint_sqanti_dir = Channel.empty()
            sample_sqanti_dirs = Channel.empty()
        } else {
            joint_sqanti = runJointSqanti(joint_bambu.gtf, analysis_annotation, analysis_reference)
            sample_sqanti = runPerSampleSqanti(sample_bambu.gtf, analysis_annotation, analysis_reference)
            joint_sqanti_dir = joint_sqanti.dir
            sample_sqanti_dirs = sample_sqanti.dir
        }

    emit:
        annotation = ref_annotation
        annotation_reference_summary = prepared_reference_annotation.summary
        unstranded_annotation = prepared_reference_annotation.unstranded
        joint_dir = joint_bambu.dir
        joint_gtf = joint_bambu.gtf
        joint_fasta = joint_fasta.fasta
        joint_transcript_counts = joint_bambu.transcript_counts
        joint_gene_counts = joint_bambu.gene_counts
        joint_transcript_rds = joint_bambu.transcript_rds
        joint_gene_rds = joint_bambu.gene_rds
        joint_metadata = joint_bambu.transcript_metadata
        sample_dirs = sample_bambu.dir
        sample_gtf = sample_bambu.gtf
        sample_fastas = sample_fastas.fasta
        sample_transcript_counts = sample_bambu.transcript_counts
        sample_gene_counts = sample_bambu.gene_counts
        sample_transcript_rds = sample_bambu.transcript_rds
        sample_gene_rds = sample_bambu.gene_rds
        sample_metadata = sample_bambu.transcript_metadata
        joint_sqanti_dir = joint_sqanti_dir
        sample_sqanti_dirs = sample_sqanti_dirs
}
