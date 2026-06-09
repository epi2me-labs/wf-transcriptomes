#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { getParams; configure_igv } from './lib/common'
include { prepare_reference } from './lib/reference'
include { transcriptome } from './subworkflows/transcriptome'
include { differential_expression } from './subworkflows/differential_expression'
include { mods } from './subworkflows/mods'



OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")


process getVersions {
    label "wf_transcriptomes"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "versions.txt"
    cpus 1
    memory "2 GB"
    input:
        path "additional_versions.txt"
    output:
        path "versions.txt"
    script:
    """
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    gffread --version | sed 's/^/gffread,/' >> versions.txt || true
    Rscript -e 'pkgs <- c("bambu", "DESeq2", "DEXSeq"); for (pkg in pkgs) {cat(pkg, as.character(packageVersion(pkg)), sep = ","); cat("\\n")}' >> versions.txt
    cat additional_versions.txt >> versions.txt
    """
}


process getSqantiVersion {
    label "wf_transcriptomes_sqanti"
    cpus 1
    memory "2 GB"
    output:
        path "versions.txt"
    script:
    """
    sqanti3_qc.py --version | sed 's/ /,/' >> versions.txt
    """
}


process getModkitVersion {
    label "modkit"
    cpus 1
    memory "2 GB"
    input:
        path "old_versions.txt"
    output:
        path "versions.txt"
    script:
    """
    cp old_versions.txt versions.txt
    modkit --version | sed 's/ /,/' >> versions.txt
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
        path mod_summary_files, stageAs: "mod_summaries/*"
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
        --mod_summary_dir mod_summaries \
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


workflow wf {
    take:
        reads
        sample_sheet
        ref_genome
        ref_annotation
    main:
        software_versions = getVersions(getModkitVersion(getSqantiVersion()))
        workflow_params = getParams()

        transcriptome_results = transcriptome(reads, ref_genome, ref_annotation, sample_sheet)
        mod_results = mods(reads, ref_genome)

        if (params.de_analysis) {
            de_results = differential_expression(
                transcriptome_results.joint_transcript_rds,
                transcriptome_results.joint_gene_rds,
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

        sample_dirs_for_report = transcriptome_results.sample_dirs
            .map { meta, sample_dir -> sample_dir }
            .collect()

        mod_summaries_for_report = mod_results.summary
            .map { alias, summary -> summary }
            .ifEmpty(OPTIONAL_FILE)
            .collect()

        sqanti_dirs_for_report = transcriptome_results.joint_sqanti_dir
            .concat(transcriptome_results.sample_sqanti_dirs.map { meta, sqanti_dir -> sqanti_dir })
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
            transcriptome_results.joint_dir.ifEmpty(OPTIONAL_FILE),
            sample_dirs_for_report,
            mod_summaries_for_report,
            sqanti_dirs_for_report,
            de_dir,
            transcriptome_results.annotation_reference_summary,
            workflow.manifest.version
        )

        results = Channel.empty()
            .concat(report.report.map { [it, null] })
            .concat(workflow_params.map { [it, null] })
            .concat(transcriptome_results.annotation_reference_summary.map { [it, "cohort/reference"] })
            .concat(transcriptome_results.unstranded_annotation.map { [it, "cohort/reference"] })
            .concat(transcriptome_results.joint_gtf.map { [it, "cohort"] })
            .concat(transcriptome_results.joint_fasta.map { [it, "cohort"] })
            .concat(transcriptome_results.joint_transcript_counts.map { [it, "cohort"] })
            .concat(transcriptome_results.joint_gene_counts.map { [it, "cohort"] })
            .concat(transcriptome_results.joint_transcript_rds.map { [it, "cohort"] })
            .concat(transcriptome_results.joint_gene_rds.map { [it, "cohort"] })
            .concat(transcriptome_results.joint_metadata.map { [it, "cohort"] })
            .concat(transcriptome_results.sample_gtf.map { meta, gtf -> [gtf, "samples/${meta.alias}"] })
            .concat(transcriptome_results.sample_fastas.map { meta, fasta -> [fasta, "samples/${meta.alias}"] })
            .concat(transcriptome_results.sample_transcript_counts.map { meta, counts -> [counts, "samples/${meta.alias}"] })
            .concat(transcriptome_results.sample_gene_counts.map { meta, counts -> [counts, "samples/${meta.alias}"] })
            .concat(transcriptome_results.sample_transcript_rds.map { meta, rds -> [rds, "samples/${meta.alias}"] })
            .concat(transcriptome_results.sample_gene_rds.map { meta, rds -> [rds, "samples/${meta.alias}"] })
            .concat(transcriptome_results.sample_metadata.map { meta, metadata -> [metadata, "samples/${meta.alias}"] })
            .concat(generated_alignment_outputs)

        if (params.de_analysis) {
            results = results.concat(de_dir.map { [it, null] })
        }
    emit:
        results = results
        bigwigs = mod_results.bigwig
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
        "sort_threads": 3,
    ]

    def checked_bam_hook = { checked_bams ->
        if (params.force_alignment) {
            // short-circuit this hook to allow for forced realignment
            return checked_bams
        }

        // determine BAMs with and without splice-aware CIGAR evidence
        ArrayList aliases_with_splice_cigars = []
        ArrayList aliases_without_splice_cigars = []
        checked_bams.each {
            def meta = it[0]
            if (meta.has_reads) {
                if (meta.has_splice_cigars) {
                    aliases_with_splice_cigars.add(meta.alias)
                } else {
                    aliases_without_splice_cigars.add(meta.alias)
                }
            }
        }
        boolean any_with_splice_cigars = aliases_with_splice_cigars
        boolean any_without_splice_cigars = aliases_without_splice_cigars

        // explode on mixed splice-aware status
        if (any_with_splice_cigars && any_without_splice_cigars) {
            log.error """
                Cannot proceed with mixed splice-aware CIGAR evidence in input BAMs.

                This workflow requires splice-aware alignment. It cannot safely realign
                only the BAMs that appear to lack splice-aware CIGARs, because the
                original alignment parameters for the other BAMs are unknown.

                It is unlikely for a splice-aware aligned transcriptome BAM to have no
                reads with N CIGAR operations. Inspect your BAMs to understand why their
                splice-aware alignment evidence differs. Either remove the BAMs without
                splice-aware CIGARs from the analysis, or use --force_alignment to realign
                all inputs.

                Samples with splice CIGARs: ${aliases_with_splice_cigars.join(', ')}
                Samples without splice CIGARs: ${aliases_without_splice_cigars.join(', ')}
                """.stripIndent().trim()
            error "Mixed splice-aware CIGAR evidence in input BAMs."
        }

        // all read-containing BAMs have splice-aware CIGAR evidence, or there are no reads
        if (!any_without_splice_cigars) {
            return checked_bams
        }

        // if we're here, all read-containing BAMs lack splice-aware CIGAR evidence
        // fiddle the metamap to set requires_alignment (if there are reads to align)
        log.warn "No input BAMs appear to contain splice-aware CIGAR evidence. " +
            "The workflow will realign all inputs."
        checked_bams.collect {
            def meta = it[0]
            def paths = it[1]
            [
                meta + [
                    requires_alignment: meta.has_reads,
                ],
                paths
            ]
        }
    }

    if (params.fastq) {
        samples = fastq_ingress([
            "input": params.fastq,
        ] + ingress_args, ref_genome)
    } else {
        samples = xam_ingress([
            "input": params.bam,
            "force_alignment": params.force_alignment,
            "checked_bam_hook": checked_bam_hook,
        ] + ingress_args, ref_genome)
    }


    sample_sheet_aliases = sample_sheet == OPTIONAL_FILE ?
        null :
        sample_sheet
            .splitCsv(header: true, quote: '"')
            .collect { it.alias }
            .findAll { it != null }
            .toSet()
    samples.subscribe { meta, xam, xai, stats ->
        if (sample_sheet_aliases != null && !sample_sheet_aliases.contains(meta.alias)) {
            throw new Exception(
                "Sample alias '${meta.alias}' was not found in the sample_sheet alias column."
            )
        }
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

    pipeline_run = wf(processed_samples, sample_sheet, ref_genome, ref_annotation)
    results = pipeline_run.results

    if (params.igv) {
        // TODO lib/ref should be responsible for writing NEW outputs to a location of our choosing
        // until then, we'll handle emission here. we'll emit (path:str, to_publish:bool) tuples for ref-related files
        // and pass those to both igv_ref_paths and results (for publishing)
        is_compressed = params.ref_genome.toLowerCase().endsWith("gz")
        if (is_compressed) {
            ref_files = prepared_reference.ref_gzidx | flatten | map {
                boolean to_publish = it.toString().startsWith("${workflow.workDir}")
                [it, to_publish]
            }
        } else {
            ref_files = ref_genome | flatten | map {
                boolean to_publish = it.toString().startsWith("${workflow.workDir}")
                [it, to_publish]
            }
        }

        // convert files set to_publish to their IGV location
        igv_ref_paths = ref_files.map {
            path, to_publish -> to_publish ? "reference/${path.getName()}" : path.toString()
        }
        publish_ref_paths = ref_files
            .filter { it[1] }  // select files set to_publish
            .map { [ it[0], "reference" ] }

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

        igv_files = igv_ref_paths
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
        results = results
            .concat(publish_ref_paths)
            .concat(igv_conf.map { [it, null] })
    }
    publishResults(results)
}


workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}

workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
