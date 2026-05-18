nextflow.enable.dsl = 2

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

// use nextflows classic rename trick
include {
    // joint
    bambuDiscover as runJointBambuDiscover
    bambuQuant as runJointBambuQuant
    bambuEmpty as runJointBambuEmpty
    collateBambuQuant as collateJointBambuQuant
    // persample
    bambuDiscover as runPerSampleBambuDiscover
    bambuQuant as runPerSampleBambuQuant
    bambuEmpty as runPerSampleBambuEmpty
    collateBambuQuant as collatePerSampleBambuQuant
} from '../modules/local/bambu_chunked'


def bambu_discover_to_quant_inputs(discover_channel) {
    discover_channel
        .map { meta, discover_dir ->
            def sampleAliases = discover_dir.resolve("samples.csv")
                .readLines()
                .drop(1)
                .findAll { it?.trim() }
                .collect { it.split(",", 2)[0] }
            tuple(
                meta,
                sampleAliases,
                discover_dir.resolve("bambu_discovered_annotations.rds"),
                discover_dir.resolve("chunks"),
                discover_dir.resolve("chunk_manifest.tsv")
            )
        }
        .splitCsv(header: true, sep: '\t', elem: 4)
        .map { meta, sample_aliases, discovered_annotation, chunks_dir, row ->
            def txCountRaw = row.annotation_tx_count?.toString()?.trim()
            Integer annotationTxCount = (!txCountRaw || txCountRaw == 'NA') ? null : txCountRaw as Integer
            tuple(
                meta,
                sample_aliases,
                row.chunk_id,
                annotationTxCount,
                chunks_dir.resolve("${row.chunk_id}.rds"),
                discovered_annotation
            )
        }
}


def bambu_filter_quant_inputs_with_warning(quant_inputs) {
    quant_inputs.filter { meta, sample_aliases, chunk_id, annotation_tx_count, chunk_rds, discovered_annotation ->
        boolean keep = annotation_tx_count == null || annotation_tx_count > 0
        if (!keep) {
            log.warn("Dropping bambu quant chunk '${chunk_id}' for '${meta.alias}' because annotation_tx_count=0")
        }
        keep
    }
}


def bambu_empty_inputs(quant_inputs) {
    quant_inputs
        .map { meta, sample_aliases, chunk_id, annotation_tx_count, chunk_rds, discovered_annotation ->
            tuple(meta.alias, meta, sample_aliases, annotation_tx_count)
        }
        .groupTuple()
        .filter { alias, metas, sample_aliases_sets, annotation_tx_counts ->
            !annotation_tx_counts.any { it == null || it > 0 }
        }
        .map { alias, metas, sample_aliases_sets, annotation_tx_counts ->
            tuple(metas[0], sample_aliases_sets[0])
        }
}


def bambu_quant_process_inputs(quant_inputs) {
    quant_inputs.map { meta, sample_aliases, chunk_id, annotation_tx_count, chunk_rds, discovered_annotation ->
        tuple(meta, chunk_id, annotation_tx_count, chunk_rds, discovered_annotation)
    }
}

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

        joint_meta = [alias: "cohort"]

        joint_discover = runJointBambuDiscover(
            alignments
                .collect(flat: false)
                .map { rows ->
                    // transform [meta, bam, bai] rows to 
                    // [meta, [alias1...aliasN], [bam1...bamN], [bai1...baiN], sample_sheet]
                    tuple(
                        joint_meta,
                        rows.collect { it[0].alias },
                        rows.collect { it[1] },
                        rows.collect { it[2] },
                        sample_sheet
                    )
                },
            analysis_annotation,
            analysis_reference
        )
        joint_quant_inputs_all = bambu_discover_to_quant_inputs(joint_discover.dir)
        joint_quant = runJointBambuQuant(
            bambu_quant_process_inputs(
                bambu_filter_quant_inputs_with_warning(joint_quant_inputs_all)
            ),
            analysis_reference
        )
        joint_bambu_real = collateJointBambuQuant(
            joint_quant.dir
                .map { meta, chunk_id, chunk_dir ->
                    tuple(meta.alias, meta, chunk_dir)
                }
                .groupTuple()
                .map { alias, metas, chunk_dirs ->
                    tuple(metas[0], chunk_dirs)
                }
        )
        joint_bambu_empty = runJointBambuEmpty(
            bambu_empty_inputs(joint_quant_inputs_all)
        )

        joint_bambu_dir = joint_bambu_real.dir.mix(joint_bambu_empty.dir)
        joint_bambu_gtf = joint_bambu_real.gtf.mix(joint_bambu_empty.gtf)
        joint_bambu_transcript_counts = joint_bambu_real.transcript_counts.mix(joint_bambu_empty.transcript_counts)
        joint_bambu_gene_counts = joint_bambu_real.gene_counts.mix(joint_bambu_empty.gene_counts)
        joint_bambu_transcript_rds = joint_bambu_real.transcript_rds.mix(joint_bambu_empty.transcript_rds)
        joint_bambu_gene_rds = joint_bambu_real.gene_rds.mix(joint_bambu_empty.gene_rds)
        joint_bambu_metadata = joint_bambu_real.transcript_metadata.mix(joint_bambu_empty.transcript_metadata)

        sample_discover = runPerSampleBambuDiscover(
            alignments.map { meta, bam, bai, stats ->
                tuple(
                    meta,
                    [meta.alias],
                    bam,
                    bai,
                    sample_sheet
                )
            },
            analysis_annotation,
            analysis_reference
        )
        sample_quant_inputs_all = bambu_discover_to_quant_inputs(sample_discover.dir)
        sample_quant = runPerSampleBambuQuant(
            bambu_quant_process_inputs(
                bambu_filter_quant_inputs_with_warning(sample_quant_inputs_all)
            ),
            analysis_reference
        )
        sample_bambu_real = collatePerSampleBambuQuant(
            sample_quant.dir
                .map { meta, chunk_id, chunk_dir ->
                    tuple(meta.alias, meta, chunk_dir)
                }
                .groupTuple()
                .map { alias, metas, chunk_dirs ->
                    tuple(metas[0], chunk_dirs)
                }
        )
        sample_bambu_empty = runPerSampleBambuEmpty(
            bambu_empty_inputs(sample_quant_inputs_all)
        )

        sample_bambu_dirs = sample_bambu_real.dir.mix(sample_bambu_empty.dir)
        sample_bambu_gtf = sample_bambu_real.gtf.mix(sample_bambu_empty.gtf)
        sample_bambu_transcript_counts = sample_bambu_real.transcript_counts.mix(sample_bambu_empty.transcript_counts)
        sample_bambu_gene_counts = sample_bambu_real.gene_counts.mix(sample_bambu_empty.gene_counts)
        sample_bambu_transcript_rds = sample_bambu_real.transcript_rds.mix(sample_bambu_empty.transcript_rds)
        sample_bambu_gene_rds = sample_bambu_real.gene_rds.mix(sample_bambu_empty.gene_rds)
        sample_bambu_metadata = sample_bambu_real.transcript_metadata.mix(sample_bambu_empty.transcript_metadata)

        joint_fasta = buildCohortTranscriptomeFasta(
            joint_bambu_real.gtf.map { meta, gtf -> gtf },
            analysis_reference
        )
        sample_fastas = buildSampleTranscriptomeFasta(sample_bambu_real.gtf, analysis_reference)

        if (params.skip_sqanti) {
            joint_sqanti_dir = Channel.empty()
            sample_sqanti_dirs = Channel.empty()
        } else {
            joint_sqanti = runJointSqanti(
                joint_bambu_real.gtf.map { meta, gtf -> gtf },
                analysis_annotation,
                analysis_reference
            )
            sample_sqanti = runPerSampleSqanti(sample_bambu_real.gtf, analysis_annotation, analysis_reference)
            joint_sqanti_dir = joint_sqanti.dir
            sample_sqanti_dirs = sample_sqanti.dir
        }

    emit:
        annotation = analysis_annotation
        annotation_reference_summary = prepared_reference_annotation.summary
        unstranded_annotation = prepared_reference_annotation.unstranded
        joint_dir = joint_bambu_dir.map { meta, dir -> dir }
        joint_gtf = joint_bambu_gtf.map { meta, gtf -> gtf }
        joint_fasta = joint_fasta.fasta
        joint_transcript_counts = joint_bambu_transcript_counts.map { meta, counts -> counts }
        joint_gene_counts = joint_bambu_gene_counts.map { meta, counts -> counts }
        joint_transcript_rds = joint_bambu_transcript_rds.map { meta, rds -> rds }
        joint_gene_rds = joint_bambu_gene_rds.map { meta, rds -> rds }
        joint_metadata = joint_bambu_metadata.map { meta, metadata -> metadata }
        sample_dirs = sample_bambu_dirs
        sample_gtf = sample_bambu_gtf
        sample_fastas = sample_fastas.fasta
        sample_transcript_counts = sample_bambu_transcript_counts
        sample_gene_counts = sample_bambu_gene_counts
        sample_transcript_rds = sample_bambu_transcript_rds
        sample_gene_rds = sample_bambu_gene_rds
        sample_metadata = sample_bambu_metadata
        joint_sqanti_dir = joint_sqanti_dir
        sample_sqanti_dirs = sample_sqanti_dirs
}
