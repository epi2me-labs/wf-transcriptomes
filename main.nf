#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { getParams; configure_igv } from './lib/common'
include { prepare_reference } from './lib/reference'
include { transcriptome_analysis } from './subworkflows/transcriptome'
include { differential_expression } from './subworkflows/differential_expression'
include { mod_analysis } from './subworkflows/mods'



OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")


process getVersions {
    label "wf_transcriptomes"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "versions.txt"
    cpus 1
    memory "2 GB"
    output:
        path "versions.txt"
    script:
    """
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    gffread --version | sed 's/^/gffread,/' >> versions.txt || true
    Rscript -e 'pkgs <- c("bambu", "DESeq2", "DEXSeq"); for (pkg in pkgs) {cat(pkg, as.character(packageVersion(pkg)), sep = ","); cat("\\n")}' >> versions.txt
    """
}


process makeReport {
    label "wf_common"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "wf-transcriptomes-report.html"
    cpus 1
    memory 8.GB
    input:
        tuple val(metadata), path(stats, stageAs: "stats_*")
        path "versions/*"
        path "params.json"
        path cohort_dir, stageAs: "cohort"
        path sample_dirs, stageAs: "samples/*"
        path sqanti_dirs, stageAs: "sqanti/*"
        path de_files
        path "annotation_reference_summary.tsv"
        val wf_version
    output:
        path "wf-transcriptomes-report.html", emit: report
    script:
        String metadata_json = new JsonBuilder(metadata).toPrettyString().replaceAll("'", "'\\\\''")
        def report_stats = (stats instanceof java.util.Collection ? stats : (stats ? [stats] : []))
            .findAll { it.name != OPTIONAL_FILE.name }
        String stats_args = report_stats ? "--stats ${report_stats.join(' ')}" : ""
        String de_args = de_files.name == OPTIONAL_FILE.name ? "" : "--de_dir de_analysis"
    """
    echo '${metadata_json}' > metadata.json
    workflow-glue report wf-transcriptomes-report.html \
        --metadata metadata.json \
        ${stats_args} \
        --cohort_dir cohort \
        --samples_dir samples \
        --sqanti_dir sqanti \
        ${de_args} \
        --versions versions \
        --params params.json \
        --ref_summary annotation_reference_summary.tsv \
        --wf_version ${wf_version}
    """
}


process publishResults {
    label "wf_common"
    memory 2.GB
    cpus 1
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    script:
    """
    """
}


workflow pipeline {
    take:
        reads
        sample_sheet
        ref_genome
        ref_annotation
    main:
        software_versions = getVersions()
        workflow_params = getParams()

        transcriptome = transcriptome_analysis(reads, ref_genome, ref_annotation, sample_sheet)
        mods = mod_analysis(reads, ref_genome)

        if (params.de_analysis) {
            de_results = differential_expression(
                transcriptome.joint_transcript_rds,
                transcriptome.joint_gene_rds,
                sample_sheet
            )
            de_dir = de_results.dir
        } else {
            de_dir = Channel.of(OPTIONAL_FILE)
        }

        report_input = reads
            .collect(flat: false)
            .map { rows ->
                def metadata = rows.collect { meta, xam, xai, stats ->
                    meta + [has_stats: stats as boolean]
                }
                def stats_rows = rows.findAll { meta, xam, xai, stats -> stats != null }
                [
                    metadata,
                    stats_rows ? stats_rows.collect { meta, xam, xai, stats -> stats } : [OPTIONAL_FILE]
                ]
            }

        sample_dirs_for_report = transcriptome.sample_dirs
            .map { meta, sample_dir -> sample_dir }
            .collect()

        sqanti_dirs_for_report = transcriptome.joint_sqanti_dir
            .concat(transcriptome.sample_sqanti_dirs.map { meta, sqanti_dir -> sqanti_dir })
            .ifEmpty(OPTIONAL_FILE)
            .collect()

        // meta.src_xam is non-null if BAMs are "passed through"
        generated_alignment_outputs = reads
            .filter { meta, bam, bai, stats -> meta.src_xam == null }
            .flatMap { meta, bam, bai, stats ->
                def outdir = "samples/${meta.alias}/alignment"
                [
                    [bam, outdir],
                    [bai, outdir],
                    [stats.resolve("bamstats.flagstat.tsv"), outdir],
                ]
            }

        report = makeReport(
            report_input,
            software_versions,
            workflow_params,
            transcriptome.joint_dir.ifEmpty(OPTIONAL_FILE),
            sample_dirs_for_report,
            sqanti_dirs_for_report,
            de_dir,
            transcriptome.annotation_reference_summary,
            workflow.manifest.version
        )

        results = Channel.empty()
            .concat(report.report.map { [it, null] })
            .concat(workflow_params.map { [it, null] })
            .concat(transcriptome.annotation_reference_summary.map { [it, "cohort/reference"] })
            .concat(transcriptome.unstranded_annotation.map { [it, "cohort/reference"] })
            .concat(transcriptome.joint_gtf.map { [it, "cohort"] })
            .concat(transcriptome.joint_fasta.map { [it, "cohort"] })
            .concat(transcriptome.joint_transcript_counts.map { [it, "cohort"] })
            .concat(transcriptome.joint_gene_counts.map { [it, "cohort"] })
            .concat(transcriptome.joint_transcript_rds.map { [it, "cohort"] })
            .concat(transcriptome.joint_gene_rds.map { [it, "cohort"] })
            .concat(transcriptome.joint_metadata.map { [it, "cohort"] })
            .concat(transcriptome.sample_gtf.map { meta, gtf -> [gtf, "samples/${meta.alias}"] })
            .concat(transcriptome.sample_fastas.map { meta, fasta -> [fasta, "samples/${meta.alias}"] })
            .concat(transcriptome.sample_transcript_counts.map { meta, counts -> [counts, "samples/${meta.alias}"] })
            .concat(transcriptome.sample_gene_counts.map { meta, counts -> [counts, "samples/${meta.alias}"] })
            .concat(transcriptome.sample_transcript_rds.map { meta, rds -> [rds, "samples/${meta.alias}"] })
            .concat(transcriptome.sample_gene_rds.map { meta, rds -> [rds, "samples/${meta.alias}"] })
            .concat(transcriptome.sample_metadata.map { meta, metadata -> [metadata, "samples/${meta.alias}"] })
            .concat(generated_alignment_outputs)

        if (params.de_analysis) {
            results = results.concat(de_dir.map { [it, null] })
        }
    emit:
        results = results
        bigwigs = mods.bigwig
}


WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)


    sample_sheet = params.sample_sheet ? file(params.sample_sheet, type: "file") : OPTIONAL_FILE
    ref_annotation = file(params.ref_annotation, type: "file")

    prepared_reference = prepare_reference(
        params.ref_genome, [
        "output_cache": false,
        "output_mmi": false,
    ])
    ref_genome = prepared_reference.ref_tuple.first()

    if (!ref_annotation.exists()) {
        throw new Exception("--ref_annotation does not exist.")
    }
    if (sample_sheet != OPTIONAL_FILE && !sample_sheet.exists()) {
        throw new Exception("--sample_sheet does not exist.")
    }

    def ingress_args = [
        "minimap2_memory": ["31GB", "62GB"],
        "minimap2_opts": params.direct_rna ? "-ax splice -uf -k14" : "-ax splice -uf",
        "alignment_threads": 12,
        "output_xam_fmt": "bam",
        "sample": params.sample,
        "sample_sheet": params.sample_sheet,
        "analyse_unclassified": params.analyse_unclassified,
        "analyse_fail": params.analyse_fail,
        "fastcat_extra_args": "",
        "required_sample_types": [],
    ]

    if (params.fastq) {
        samples = fastq_ingress([
            "input": params.fastq,
        ] + ingress_args, ref_genome)
    } else {
        samples = xam_ingress([
            "input": params.bam,
            "force_alignment": params.force_alignment,
        ] + ingress_args, ref_genome)
    }
    analysis_samples = samples
        .filter { meta, xam, xai, stats ->
            boolean is_excluded = false
            String excluded_reason = null
            if (meta.n_primary == 0) {
                excluded_reason = "has no reads"
                is_excluded = true
            } else if (meta.n_primary == null) {
                excluded_reason = "was not found during ingress"
                is_excluded = true
            }
            if (is_excluded) {
                log.warn("Sample ${meta.alias} ${excluded_reason} and will be excluded from transcriptome analysis.")
                if (params.de_analysis) {
                    throw new Exception("""\
                    Differential gene expression and differential transcript analyses
                    requires all the samples present in the sample sheet to have reads.
                    """.stripIndent())
                }
            }
            return !is_excluded
        }
        .ifEmpty {
            throw new Exception("No samples with reads were available for transcriptome analysis.")
        }

    processed_samples = analysis_samples

    pipeline_run = pipeline(processed_samples, sample_sheet, ref_genome, ref_annotation)
    results = pipeline_run.results

    reference_basename = file(params.ref_genome).getName()
    if (params.igv) {
        results = results
            .concat(ref_genome.map { fasta, faidx -> [fasta, "igv_reference"] })
            .concat(ref_genome.map { fasta, faidx -> [faidx, "igv_reference"] })

        is_compressed = params.ref_genome.toLowerCase().endsWith("gz")

        if (is_compressed) {
            // ref files are directly publish into output
            igv_files = Channel.of("${reference_basename}")
            igv_index_paths = prepared_reference.ref_gzidx.map {
                    fasta, faidx, gzidx -> "${faidx.getName()}"
                }
                .concat(prepared_reference.ref_gzidx.map {
                    fasta, faidx, gzidx -> "${gzidx.getName()}"
                })
        } else {
            igv_files = Channel.of("igv_reference/${reference_basename}")
            igv_index_paths = ref_genome.map { fasta, faidx -> "igv_reference/${faidx.getName()}"}
        }

        igv_alignment_paths = processed_samples
            .map { meta, bam, bai, stat -> [
                meta.src_xam ?: "${meta.alias},samples/${meta.alias}/alignment/reads.bam",
                meta.src_xai ?: "${meta.alias},samples/${meta.alias}/alignment/reads.bam.bai"
            ] }
            .flatten()

        // convert [alias0, [bw00...bw0N]] to [alias0, bw00] ... [aliasN, bwNN]
        // allowing for [aliasM, bwM0] if only one bw is output because ... nextflow
        // and use the anticipated output location
        igv_bigwigs = pipeline_run.bigwigs
            .flatMap { alias, paths ->
                (paths instanceof List ? paths : [paths]).collect { path -> "${alias},samples/${alias}/mods/${path.name}" }
            }

        igv_files = igv_files
            .concat(igv_index_paths)
            .concat(igv_alignment_paths)
            .concat(igv_bigwigs)
            .collectFile(name: "igv-files.txt", newLine: true, sort: false)

        igv_conf = configure_igv(
            igv_files,
            "",
            [displayMode: "SQUISHED", colorBy: "strand"],
            [:],
            false
        )
        results = results.concat(igv_conf.map { [it, null] })
    }
    publishResults(results)
}


workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}

workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
