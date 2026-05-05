#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { getParams; configure_igv } from './lib/common'
include { transcriptome_analysis } from './subworkflows/transcriptome'
include { differential_expression } from './subworkflows/differential_expression'


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


process collectIngressResultsInDir {
    label "wf_common"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path(reads, stageAs: "reads/*"), path(stats, stageAs: "stats/*")
    output:
        path "out/*", emit: out
    script:
        String outdir = "out/${meta["alias"].replaceAll("'", "'\\\\''")}"
        String meta_json = new JsonBuilder(meta).toPrettyString().replaceAll("'", "'\\\\''")
        String stats_arg = stats.fileName.name == OPTIONAL_FILE.name ? "" : stats
        def read_args = reads instanceof java.util.ArrayList ? reads.join(" ") : reads
    """
    mkdir -p '${outdir}'
    echo '${meta_json}' > metamap.json
    mv metamap.json ${read_args} ${stats_arg} '${outdir}'
    """
}


process preprocess_reads {
    label "wf_transcriptomes_pychopper"
    cpus { params.threads ?: 4 }
    memory "8 GB"
    input:
        tuple val(meta), path(reads, stageAs: "reads/*")
    output:
        tuple val(meta.alias), path("${meta.alias}_pychopper_output/${meta.alias}_full_length_reads.fastq"), emit: full_len_reads
        tuple val(meta.alias), path("${meta.alias}_pychopper_output"), emit: dir
    script:
        String backend = params.pychopper_backend ?: "edlib"
        String extra = params.pychopper_opts ?: ""
        String cdna_kit = params.cdna_kit ? params.cdna_kit.tokenize("-")[-1] : ""
        def read_args = reads instanceof java.util.ArrayList ? reads.join(" ") : reads
    """
    cat ${read_args} > seqs.fastq.gz
    pychopper -t ${task.cpus} -k ${cdna_kit} -m ${backend} ${extra} \
        seqs.fastq.gz "${meta.alias}_full_length_reads.fastq"

    awk '
        BEGIN { FS = OFS = "\\t" }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                idx[\$i] = i
            }
            print "Classification", "Value"
            next
        }
        \$idx["Category"] == "Classification" {
            print \$idx["Name"], \$idx["Value"]
        }
    ' pychopper.tsv > pychopper_summary.tsv

    mkdir -p "${meta.alias}_pychopper_output"
    find . -maxdepth 1 -mindepth 1 \
        ! -name "seqs.fastq.gz" \
        ! -name "${meta.alias}_pychopper_output" \
        -exec mv -t "${meta.alias}_pychopper_output" {} +
    """
}


process makeReport {
    label "wf_common"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "wf-transcriptomes-report.html"
    input:
        tuple val(metadata), path(stats, stageAs: "stats_*")
        path "versions/*"
        path "params.json"
        path alignment_stats, stageAs: "alignment_stats/*"
        path cohort_dir, stageAs: "cohort"
        path sample_dirs, stageAs: "samples/*"
        path pychopper_dirs, stageAs: "pychopper/*"
        path sqanti_dirs, stageAs: "sqanti/*"
        path de_files
        val wf_version
    output:
        path "wf-transcriptomes-report.html", emit: report
    script:
        String metadata_json = new JsonBuilder(metadata).toPrettyString().replaceAll("'", "'\\\\''")
        def report_stats = stats instanceof java.util.Collection ? stats : (stats ? [stats] : [])
        def report_pychopper = pychopper_dirs instanceof java.util.Collection ? pychopper_dirs : (pychopper_dirs ? [pychopper_dirs] : [])
        String stats_args = report_stats ? "--stats ${report_stats.join(' ')}" : ""
        String pychopper_args = report_pychopper.find { it.name != OPTIONAL_FILE.name } ? "--pychopper_dir pychopper" : ""
        String de_args = de_files.name == OPTIONAL_FILE.name ? "" : "--de_dir de_analysis"
    """
    echo '${metadata_json}' > metadata.json
    workflow-glue report wf-transcriptomes-report.html \
        --metadata metadata.json \
        ${stats_args} \
        --alignment_stats_dir alignment_stats \
        --cohort_dir cohort \
        --samples_dir samples \
        ${pychopper_args} \
        --sqanti_dir sqanti \
        ${de_args} \
        --versions versions \
        --params params.json \
        --wf_version ${wf_version}
    """
}


process publishResults {
    label "wf_common"
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    """
}


def coerceBooleanParam(value) {
    if (value == null || value instanceof Boolean) {
        return value
    }
    if (value instanceof CharSequence) {
        switch (value.toString().trim().toLowerCase()) {
            case "true":
            case "1":
            case "yes":
                return true
            case "false":
            case "0":
            case "no":
                return false
        }
    }
    return value
}


[
    "help",
    "version",
    "igv",
    "direct_rna",
    "cdna_preprocess",
    "de_analysis",
    "analyse_unclassified",
    "analyse_fail",
    "skip_sqanti",
    "sqanti_skip_orf",
    "disable_ping",
    "monochrome_logs",
    "validate_params",
    "show_hidden_params",
].each { name ->
    params[name] = coerceBooleanParam(params[name])
}

[
    "keep_unaligned",
    "return_fastq",
    "per_read_stats",
    "allow_multiple_basecall_models",
].each { name ->
    if (params.wf?.containsKey(name)) {
        params.wf[name] = coerceBooleanParam(params.wf[name])
    }
}


workflow pipeline {
    take:
        reads
        sample_sheet
        ref_genome
        ref_annotation
        pychopper_dirs
    main:
        software_versions = getVersions()
        workflow_params = getParams()

        ingress_results = collectIngressResultsInDir(
            reads.map { meta, sample_reads, stats ->
                [meta, sample_reads, stats ?: OPTIONAL_FILE]
            }
        )

        transcriptome = transcriptome_analysis(reads, ref_genome, ref_annotation, sample_sheet)

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
            .map { meta, sample_reads, stats ->
                ["all_samples", meta + [has_stats: stats as boolean], stats]
            }
            .groupTuple()
            .map { report_group, metas, stats ->
                [metas, stats.findAll { it != null }]
            }

        alignment_stats = transcriptome.alignments
            .map { meta, bam, bai, flagstat -> flagstat }
            .collect()

        sample_dirs_for_report = transcriptome.sample_dirs
            .map { meta, sample_dir -> sample_dir }
            .collect()

        pychopper_dirs_for_report = pychopper_dirs
            .map { meta, pychopper_dir -> pychopper_dir }
            .ifEmpty(OPTIONAL_FILE)
            .collect()

        sqanti_dirs_for_report = transcriptome.joint_sqanti_dir
            .concat(transcriptome.sample_sqanti_dirs.map { meta, sqanti_dir -> sqanti_dir })
            .ifEmpty(OPTIONAL_FILE)
            .collect()

        report = makeReport(
            report_input,
            software_versions,
            workflow_params,
            alignment_stats,
            transcriptome.joint_dir,
            sample_dirs_for_report,
            pychopper_dirs_for_report,
            sqanti_dirs_for_report,
            de_dir,
            workflow.manifest.version
        )

        results = Channel.empty()
            .concat(ingress_results.out.map { [it, "ingress_results"] })
            .concat(report.report.map { [it, null] })
            .concat(workflow_params.map { [it, null] })
            .concat(pychopper_dirs.map { meta, pychopper_dir -> [pychopper_dir, "ingress_results/${meta.alias}"] })
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
            .concat(transcriptome.alignments.map { meta, bam, bai, flagstat -> [bam, "cohort/alignments"] })
            .concat(transcriptome.alignments.map { meta, bam, bai, flagstat -> [bai, "cohort/alignments"] })
            .concat(transcriptome.alignments.map { meta, bam, bai, flagstat -> [flagstat, "cohort/alignments"] })
            .concat(transcriptome.joint_sqanti_dir.map { [it, "cohort"] })
            .concat(transcriptome.sample_sqanti_dirs.map { meta, sqanti_dir -> [sqanti_dir, "samples/${meta.alias}"] })

        reference_basename = file(params.ref_genome).getName()
        if (params.igv) {
            results = results
                .concat(transcriptome.reference.map { [it, "igv_reference"] })
                .concat(transcriptome.reference_fai.map { [it, "igv_reference"] })
                .concat(transcriptome.reference_gzi.map { [it, "igv_reference"] })

            igv_index_paths = transcriptome.reference_fai
                .map { "igv_reference/${it.getName()}" }
                .concat(transcriptome.reference_gzi.map { "igv_reference/${it.getName()}" })

            igv_alignment_paths = transcriptome.alignments
                .map { meta, bam, bai, flagstat -> [
                    "cohort/alignments/${bam.getName()}",
                    "cohort/alignments/${bai.getName()}"
                ] }
                .flatten()

            igv_files = Channel.of("igv_reference/${reference_basename}")
                .concat(igv_index_paths)
                .concat(igv_alignment_paths)
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

        if (params.de_analysis) {
            results = results.concat(de_dir.map { [it, null] })
        }
    emit:
        results = results
}


WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)

    if (params.containsKey("ref_transcriptome")) {
        throw new Exception("--ref_transcriptome has been removed. Use --transcriptome_mode fixed_annotation with --ref_genome and --ref_annotation.")
    }
    if (params.containsKey("transcriptome_source")) {
        throw new Exception("--transcriptome_source has been removed. Use --transcriptome_mode with either discover or fixed_annotation.")
    }

    if (!!params.fastq == !!params.bam) {
        throw new Exception("Provide exactly one of --fastq or --bam.")
    }
    if (!params.ref_genome) {
        throw new Exception("Provide --ref_genome.")
    }
    if (!params.ref_annotation) {
        throw new Exception("Provide --ref_annotation.")
    }
    if (!(params.transcriptome_mode in ["discover", "fixed_annotation"])) {
        throw new Exception("--transcriptome_mode must be one of: discover, fixed_annotation.")
    }
    if (params.direct_rna && params.cdna_preprocess) {
        throw new Exception("--cdna_preprocess cannot be used together with --direct_rna.")
    }
    if (params.de_analysis && !params.sample_sheet) {
        throw new Exception("Provide --sample_sheet when running with --de_analysis.")
    }

    sample_sheet = params.sample_sheet ? file(params.sample_sheet, type: "file") : OPTIONAL_FILE
    ref_genome = file(params.ref_genome, type: "file")
    ref_annotation = file(params.ref_annotation, type: "file")

    if (!ref_genome.exists()) {
        throw new Exception("--ref_genome does not exist.")
    }
    if (!ref_annotation.exists()) {
        throw new Exception("--ref_annotation does not exist.")
    }
    if (sample_sheet != OPTIONAL_FILE && !sample_sheet.exists()) {
        throw new Exception("--sample_sheet does not exist.")
    }

    def samples
    if (params.fastq) {
        samples = fastq_ingress([
            "input": params.fastq,
            "sample": params.sample,
            "sample_sheet": params.sample_sheet,
            "analyse_unclassified": params.analyse_unclassified,
            "analyse_fail": params.analyse_fail,
            "fastcat_extra_args": "",
            "required_sample_types": [],
            "fastq_chunk": params.fastq_chunk,
            "per_read_stats": params.wf.per_read_stats,
            "allow_multiple_basecall_models": params.wf.allow_multiple_basecall_models,
        ])
    } else {
        samples = xam_ingress([
            "input": params.bam,
            "sample": params.sample,
            "sample_sheet": params.sample_sheet,
            "analyse_unclassified": params.analyse_unclassified,
            "analyse_fail": params.analyse_fail,
            "keep_unaligned": params.wf.keep_unaligned,
            "return_fastq": params.wf.return_fastq,
            "fastq_chunk": params.fastq_chunk,
            "per_read_stats": params.wf.per_read_stats,
            "allow_multiple_basecall_models": params.wf.allow_multiple_basecall_models,
        ])
    }

    decorated_samples = samples
        .map { meta, fname, stats -> [meta["group_key"], meta, fname, stats] }
        .groupTuple()
        .map { key, metas, fnames, statss ->
            if (fnames[0] == null) {
                fnames = null
            }
            [
                metas[0] + ["group_index": metas.collect { it["group_index"] }],
                fnames,
                statss[0]
            ]
        }

    analysis_samples = decorated_samples
        .filter { meta, sample_reads, stats ->
            if (meta.n_seqs == 0) {
                log.warn("Sample ${meta.alias} has no reads - excluded from transcriptome analysis.")
                return false
            }
            true
        }
        .ifEmpty {
            throw new Exception("No samples with reads were available for transcriptome analysis.")
        }

    pychopper_results = Channel.empty()
    processed_samples = analysis_samples

    if (params.cdna_preprocess) {
        grouped_samples = analysis_samples.branch { meta, sample_reads, stats ->
            to_process: sample_reads != null
            passthrough: sample_reads == null
        }

        preprocessed = preprocess_reads(
            grouped_samples.to_process
                .map { meta, sample_reads, stats -> [meta, sample_reads] }
        )

        processed_samples = grouped_samples.passthrough
            .mix(
                grouped_samples.to_process
                    .map { meta, sample_reads, stats -> [meta.alias, meta, stats] }
                    .join(preprocessed.full_len_reads)
                    .map { alias, meta, stats, full_length_reads ->
                        [meta, full_length_reads, stats]
                    }
            )

        pychopper_results = grouped_samples.to_process
            .map { meta, sample_reads, stats -> [meta.alias, meta] }
            .join(preprocessed.dir)
            .map { alias, meta, pychopper_dir ->
                [meta, pychopper_dir]
            }
    }

    pipeline_run = pipeline(processed_samples, sample_sheet, ref_genome, ref_annotation, pychopper_results)
    publishResults(pipeline_run.results)
}


workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}

workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
