nextflow.enable.dsl = 2

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")


process decompressReference {
    label "wf_transcriptomes"
    cpus 1
    memory "2 GB"
    input:
        path compressed_ref
    output:
        path "${compressed_ref.baseName}", emit: decompressed_ref
    script:
    """
    gzip -dc "${compressed_ref}" > "${compressed_ref.baseName}"
    """
}


process decompressAnnotation {
    label "wf_transcriptomes"
    cpus 1
    memory "2 GB"
    input:
        path compressed_annotation
    output:
        path "${compressed_annotation.baseName}", emit: decompressed_annotation
    script:
    """
    gzip -dc "${compressed_annotation}" > "${compressed_annotation.baseName}"
    """
}


process normaliseAnnotationToGtf {
    label "wf_transcriptomes"
    cpus 1
    memory "4 GB"
    input:
        path annotation
    output:
        path "annotation.gtf", emit: gtf
    script:
    """
    gffread -T "${annotation}" -o annotation.gtf
    """
}


process filterUnstrandedAnnotation {
    label "wf_transcriptomes"
    cpus 1
    memory "2 GB"
    input:
        path "annotation.gtf"
    output:
        tuple stdout, path("annotation_stranded.gtf"), emit: filtered
    script:
    """
    awk '
        BEGIN { OFS = "\\t" }
        /^#/ { print; next }
        (\$7 == "+" || \$7 == "-") { print; next }
        { print \$0 >> "unstranded.gtf" }
    ' annotation.gtf > annotation_stranded.gtf

    if [ -s unstranded.gtf ]; then
        echo "Warning: Unstranded entries found and excluded from differential expression analysis."
        echo "If running with reference-guided transcriptome source, consider increasing read depth to reduce unstranded annotations."
        echo "If running with precomputed transcriptome source, ensure your ref_annotation gtf file contains only '+' or '-' strand entries."
        echo "A sample of unstranded entries:"
        head -n 20 unstranded.gtf
    fi
    """
}


process validateReferenceAnnotation {
    label "wf_transcriptomes"
    cpus 1
    memory "2 GB"
    input:
        path "annotation.gtf"
        path "reference.fasta"
    output:
        stdout emit: warnings
    script:
    """
    grep -v '^#' annotation.gtf | cut -f1 | sort -u > annotation_ids.txt
    awk '/^>/ {print substr(\$1,2)}' reference.fasta | sort -u > reference_ids.txt
    matches=\$(comm -12 annotation_ids.txt reference_ids.txt || true)
    only_in_annotation=\$(comm -23 annotation_ids.txt reference_ids.txt || true)
    only_in_reference=\$(comm -13 annotation_ids.txt reference_ids.txt || true)

    if [[ -z "\$matches" ]]; then
        echo "ERROR: No overlapping seqnames were found between the reference annotation and the reference genome." >&2
        echo "Annotation ID examples:" >&2
        head -n 5 annotation_ids.txt >&2
        echo "Reference ID examples:" >&2
        head -n 5 reference_ids.txt >&2
        exit 78
    fi

    if [[ -n "\$only_in_annotation" ]]; then
        echo "Warning: Some seqnames are present in the annotation but not the genome:"
        echo "\$only_in_annotation" | head -n 5
    fi

    if [[ -n "\$only_in_reference" ]]; then
        echo "Warning: Some seqnames are present in the genome but not the annotation:"
        echo "\$only_in_reference" | head -n 5
    fi
    """
}


process buildMinimapIndex {
    label "wf_transcriptomes"
    cpus 4
    memory "32 GB"
    input:
        path reference
    output:
        tuple path("genome.mmi"), path(reference), emit: index
    script:
        String index_opts = params.minimap2_index_opts ?: ""
    """
    minimap2 -t ${task.cpus} ${index_opts} -d genome.mmi "${reference}"
    """
}


process alignReads {
    label "wf_transcriptomes"
    cpus { params.threads ?: 4 }
    memory "16 GB"
    input:
        tuple val(meta), path(reads), path(index), path(reference)
    output:
        tuple val(meta),
            path("${meta.alias}.aligned.sorted.bam"),
            path("${meta.alias}.aligned.sorted.bam.bai"),
            path("${meta.alias}.flagstat.txt"),
            emit: bam
    script:
        String preset = params.direct_rna ? "-ax splice -uf -k14" : "-ax splice -uf"
        String extra = params.minimap2_opts ?: ""
        def read_args = reads instanceof java.util.ArrayList ? reads.join(" ") : reads
        int sort_threads = Math.max(task.cpus as int - 1, 1)
    """
    minimap2 -t ${task.cpus} ${preset} ${extra} "${index}" ${read_args} \
        | samtools sort --write-index -@ ${sort_threads} \
            -o "${meta.alias}.aligned.sorted.bam##idx##${meta.alias}.aligned.sorted.bam.bai" -
    samtools flagstat "${meta.alias}.aligned.sorted.bam" > "${meta.alias}.flagstat.txt"
    """
}


process runJointBambu {
    label "wf_transcriptomes"
    cpus { params.threads ?: 4 }
    memory "32 GB"
    input:
        path bam_files, stageAs: "bams/*"
        path bam_indexes, stageAs: "bams/*"
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
        String sample_sheet_arg = sample_sheet.name == OPTIONAL_FILE.name ? "" : "--sample_sheet ${sample_sheet}"
        String ndr_arg = params.ndr != null ? "--ndr ${params.ndr}" : ""
    """
    run_bambu.R \
        --bam_dir bams \
        ${sample_sheet_arg} \
        --annotation "${annotation}" \
        --genome "${reference}" \
        --transcriptome_mode "${params.transcriptome_mode}" \
        --threads ${task.cpus} \
        ${ndr_arg} \
        --out_dir cohort
    """
}


process runPerSampleBambu {
    label "wf_transcriptomes"
    cpus { params.threads ?: 4 }
    memory "24 GB"
    input:
        tuple val(meta), path(bam), path(bai), path(flagstat)
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
    """
    run_bambu.R \
        --bam_path "${bam}" \
        --sample_alias "${meta.alias}" \
        --annotation "${annotation}" \
        --genome "${reference}" \
        --transcriptome_mode "${params.transcriptome_mode}" \
        --threads ${task.cpus} \
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
    awk '
        BEGIN { OFS = "\\t" }
        /^#/ { print; next }
        (\$7 == "+" || \$7 == "-") { print; next }
    ' "${annotation}" > reference_annotation_stranded.gtf
    sqanti3_qc.py \
        --isoforms "${gtf}" \
        --refGTF "reference_annotation_stranded.gtf" \
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
    awk '
        BEGIN { OFS = "\\t" }
        /^#/ { print; next }
        (\$7 == "+" || \$7 == "-") { print; next }
    ' "${annotation}" > reference_annotation_stranded.gtf
    sqanti3_qc.py \
        --isoforms "${gtf}" \
        --refGTF "reference_annotation_stranded.gtf" \
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


process faidx {
    label "wf_transcriptomes"
    cpus 1
    memory "4 GB"
    input:
        path ref
    output:
        path("${ref}.fai"), emit: fai
    script:
    """
    samtools faidx "${ref}"
    """
}


process gzFaidx {
    label "wf_transcriptomes"
    cpus 1
    memory "4 GB"
    errorStrategy "ignore"
    input:
        path ref
    output:
        tuple path("${ref}.fai"), path("${ref}.gzi"), emit: indexes
    script:
    """
    samtools faidx "${ref}"
    """
}


workflow transcriptome_analysis {
    take:
        reads
        ref_genome
        ref_annotation
        sample_sheet
    main:
        analysis_reference = ref_genome
        if (params.ref_genome.toLowerCase().endsWith("gz")) {
            decompressed_reference = decompressReference(ref_genome)
            analysis_reference = decompressed_reference.decompressed_ref
        }

        analysis_annotation = ref_annotation
        if (params.ref_annotation.toLowerCase().endsWith("gz")) {
            decompressed_annotation = decompressAnnotation(ref_annotation)
            analysis_annotation = decompressed_annotation.decompressed_annotation
        }
        if (params.ref_annotation.toLowerCase() ==~ /.*\.gff3?(\.gz)?$/) {
            normalised_annotation = normaliseAnnotationToGtf(analysis_annotation)
            analysis_annotation = normalised_annotation.gtf
        }
        filtered_annotation = filterUnstrandedAnnotation(analysis_annotation)
        analysis_annotation = filtered_annotation.filtered.map { warning_text, annotation ->
            if (warning_text?.trim()) {
                log.warn(warning_text.trim())
            }
            annotation
        }

        validateReferenceAnnotation(analysis_annotation, analysis_reference).map { stdoutput ->
            if (stdoutput) {
                log.warn(stdoutput)
            }
        }

        genome_index = buildMinimapIndex(analysis_reference)

        aligned = alignReads(
            reads.map { meta, sample_reads, stats -> [meta, sample_reads] }
                .combine(genome_index.index)
        )

        joint_bambu = runJointBambu(
            aligned.bam.map { meta, bam, bai, flagstat -> bam }.collect(),
            aligned.bam.map { meta, bam, bai, flagstat -> bai }.collect(),
            sample_sheet,
            analysis_annotation,
            analysis_reference
        )

        sample_bambu = runPerSampleBambu(aligned.bam, analysis_annotation, analysis_reference)

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

        compressed_reference = params.ref_genome.toLowerCase().endsWith("gz")
        if (compressed_reference) {
            ref_indexes = gzFaidx(ref_genome)
            reference_fai = ref_indexes.indexes.map { fai, gzi -> fai }
            reference_gzi = ref_indexes.indexes.map { fai, gzi -> gzi }
        } else {
            ref_indexes = faidx(ref_genome)
            reference_fai = ref_indexes.fai
            reference_gzi = Channel.empty()
        }
    emit:
        reference = ref_genome
        annotation = ref_annotation
        alignments = aligned.bam
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
        reference_fai = reference_fai
        reference_gzi = reference_gzi
}
