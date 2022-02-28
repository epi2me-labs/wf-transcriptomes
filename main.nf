#!/usr/bin/env nextflow

/* This workflow is a adapted from two previous pipeline written in Snakemake:
- https://github.com/nanoporetech/pipeline-nanopore-ref-isoforms
- https://github.com/nanoporetech/pipeline-nanopore-denovo-isoforms
*/

import groovy.json.JsonBuilder;
import nextflow.util.BlankSeparatedList;
import java.util.ArrayList;
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress' 
include { start_ping; end_ping } from './lib/ping'
include { reference_assembly } from './reference_assembly'
include { denovo_assembly } from './denovo_assembly'


process summariseConcatReads {
    // concatenate fastq and fastq.gz in a dir write stats

    label "isoforms"
    cpus 1
    input:
        tuple path(directory), val(sample_id), val(type)
    output:
        tuple val(sample_id), path("${sample_id}.fastq"), emit: input_reads
        tuple val(sample_id), path('*.stats'), emit: summary
    script:
    """
    fastcat -s ${sample_id} -r ${sample_id}.stats -x ${directory} >  ${sample_id}.fastq
    """
}

process getVersions {
    label "isoforms"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import aplanat; print(f'aplanat,{aplanat.__version__}')" >> versions.txt
    python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    bedtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    python -c "import pychopper; print(f'pychopper,{pychopper.__version__}')" >> versions.txt
    gffread --version | sed 's/^/gffread,/' >> versions.txt
    seqkit version | head -n 1 | sed 's/ /,/' >> versions.txt
    stringtie --version | sed 's/^/stringtie,/' >> versions.txt
    gffcompare --version | head -n 1 | sed 's/ /,/' >> versions.txt
    spoa --version | sed 's/^/spoa,/' >> versions.txt
    isONclust2 version | sed 's/ version: /,/' >> versions.txt
    """
}


process getParams {
    label "isoforms"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process preprocess_reads {
    /* 
    Concatenate reads from a sample directory.
    Optionally classify, trim, and orient cDNA reads using cdna_classifier from pychopper
    */
 
    label "isoforms"
    cpus 4

    input:
        tuple val(sample_id), path(input_reads)
    output:
         tuple val(sample_id), path("${sample_id}_full_length_reads.fq"), emit: full_len_reads
         tuple val(sample_id), path('*.tsv'),  emit: report
    script:
        """
        cdna_classifier.py -t ${params.threads} ${params.pychopper_opts} ${input_reads} ${sample_id}_full_length_reads.fq
        mv cdna_classifier_report.tsv ${sample_id}_cdna_classifier_report.tsv
        generate_pychopper_stats.py --data ${sample_id}_cdna_classifier_report.tsv --output .
        """
}

process build_minimap_index{
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
        minimap2 -t ${params.threads} ${params.minimap_index_opts} -I 1000G -d "genome_index.mmi" ${reference}
    """
}

process split_bam{
    /*
    Partition BAM file into loci or bundles with `params.bundle_min_reads` minimum size
    If no splitting required, just create single symbolic link to a single bundle.

    Output tuples containing `sample_id` so bundles can be combined later in th pipeline.
    */

    label 'isoforms'
    cpus params.threads

    input:
        tuple val(sample_id), path(bam)
    output:
        tuple val(sample_id), path('*.bam'), emit: bundles
    script:
    if (params["bundle_min_reads"] != false)
        """
        seqkit bam -j ${params.threads} -N ${params.bundle_min_reads} ${bam} -o  bam_bundles/
        mv bam_bundles/* .
        for f in *:*; do mv -v "\$f" \$(echo "\$f" | tr ':' '-'); done
        """
    else
        """
        mkdir -p ./${sample_id}_bam_bundles
        ln -s ${bam} ${sample_id}_bam_bundles-000000000-ALL-0-1_bundle.bam
        """
}


process assemble_transcripts{
    /*
    Assemble transcripts using stringtie.
    Take aligned reads in bam format that may be a chunk of a larger alignment file.
    Optionally use reference annotation to guide assembly.

    Output gff annotation files in a tuple with `sample_id` for combining into samples late rin the pipeline.
    */
    label 'isoforms'
    cpus params.threads

    input:
        tuple val(sample_id), path(bam)
        path ref_annotation
    output:
        tuple val(sample_id), path('*.gff'), emit: gff_bundles
    script:
        def out_filename = bam.name.replaceFirst(~/\.[^\.]+$/, '') + "_${sample_id}.gff"
        def G_FLAG = ref_annotation.name.startsWith('OPTIONAL_FILE') ? '' : "-G ${ref_annotation}"
    """
    stringtie --rf ${G_FLAG} -L -v -A gene_abund.tab -p ${params.threads} ${params.stringtie_opts} -o  ${out_filename} \
    ${bam} 2>/dev/null
    """
}


process merge_gff_bundles{
    /*
    Merge gff bundles into a single gff file.
    */
    label 'isoforms'

    input:
        tuple val(sample_id), path (gff_bundle)
    output:
        tuple val(sample_id), path('*.gff'), emit: gff
    script:
    def merged_gff = "transcripts_${sample_id}.gff"
    """
    echo '##gff-version 2' >> $merged_gff;
    echo '#pipeline-nanopore-isoforms: stringtie' >> $merged_gff;

    for fn in ${gff_bundle};
    do
        grep -v '#' \$fn >> $merged_gff

    done
    """
}

process run_gffcompare{
    /*
    Compare query and reference annotations.
    If ref_annotation is an optional file, just make an empty directory to satisfy
    the requirements of the downstream processes.
    */

    label 'isoforms'

    input:
       tuple val(sample_id), path(query_annotation)
       path ref_annotation
    output:
        tuple val(sample_id), path("${sample_id}_gffcompare"), emit: gffcmp_dir
    script:
    def out_dir = "${sample_id}_gffcompare"

    if ( ref_annotation.name.startsWith('OPTIONAL_FILE') ){
        """
        mkdir $out_dir
        """
    } else {
        """
        mkdir $out_dir
        echo "Doing comparison of reference annotation: ${ref_annotation} and the query annotation"

       gffcompare -o ${out_dir}/str_merged -r ${ref_annotation} \
        ${params.gffcompare_opts} ${query_annotation}

        generate_tracking_summary.py --tracking $out_dir/str_merged.tracking \
        --output_dir ${out_dir} --annotation ${ref_annotation}

       mv *.tmap $out_dir
       mv *.refmap $out_dir
        """
    }
}


process get_transcriptome{
        /*
        Write out a transcriptome file based on the query gff annotations.
        TODO: Do we need to touch merged_transcriptome or can we just pass it?
        */
        label 'isoforms'

        input:
            tuple val(sample_id), path(transcripts_gff), path(gffcmp_dir), path(reference_seq)
        output:
            tuple val(sample_id), path("*.fas"), emit: transcriptome

        script:
        def transcriptome = "${sample_id}_transcriptome.fas"
        def merged_transcriptome = "${sample_id}_merged_transcriptome.fas"
        """

        gffread -g ${reference_seq} -w ${transcriptome} ${transcripts_gff}
        if  [ "\$(ls -A $gffcmp_dir)" ];
        then
            echo "Yes"
            gffread -F -g ${reference_seq} -w ${merged_transcriptome} \
            $gffcmp_dir/str_merged.annotated.gtf
        fi
        """
}

process makeReport {

    label "isoforms"

    input:
        path versions
        path "params.json"
        val denovo
        tuple val(sample_ids),
                path(seq_summaries),
                path(aln_stats),
                path(gffcmp_dir),
                path(cdna_class_report),
                path(gff_annotation)
    output:
        path("wf-isoforms-*.html"), emit: report
    script:
        // Convert the sample_id arrayList.
        sids = new BlankSeparatedList(sample_ids)

        def report_name = "wf-isoforms-report.html"
        def OPT_ALN = denovo ?  '' : "--alignment_stats ${aln_stats}"
        def OPT_DENOVO = denovo ? "--denovo" : ''
    """
    report.py --report $report_name \
    --versions $versions \
    --params params.json \
    $OPT_ALN \
    --pychop_report $cdna_class_report \
    --sample_ids $sids \
    --summaries $seq_summaries \
    --gffcompare_dir $gffcmp_dir \
    --gff_annotation $gff_annotation \
    --transcript_table_cov_thresh $params.transcript_table_cov_thresh \
    $OPT_DENOVO
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}

// workflow module
workflow pipeline {
    take:
        reads
        ref_genome
        ref_annotation
    main:
        map_sample_ids_cls = {it ->
        /* Harmonize tuples
        output:
            tuple val(sample_id), path('*.gff')
        When there are multiple paths, will emit:
            [sample_id, [path, path ..]]
        when there's a single path, this:
            [sample_id, path]
        This closure makes both cases:
            [[sample_id, path][sample_id, path]].
        */
            if (it[1].getClass() != java.util.ArrayList){
                // If only one path, `it` will be [sample_id, path]
                return [it]
            }
            l = [];
            for (x in it[1]){
                l.add(tuple(it[0], x))
            }
            return l
        }

        summariseConcatReads(reads)
        sample_ids = summariseConcatReads.out.summary.flatMap({it -> it[0]})

        software_versions = getVersions()
        workflow_params = getParams()
        preprocess_reads(summariseConcatReads.out.input_reads)

        if (params.denovo){
            println("Doing de novo assembly")
             m = denovo_assembly(preprocess_reads.out.full_len_reads, ref_genome)

        } else {
            build_minimap_index(ref_genome)
            println("Doing reference based transcript analysis")
            m = reference_assembly(build_minimap_index.out.index, ref_genome, preprocess_reads.out.full_len_reads)
        }

        split_bam(m.bam)

        assemble_transcripts(split_bam.out.bundles.flatMap(map_sample_ids_cls), ref_annotation)

        merge_gff_bundles(assemble_transcripts.out.gff_bundles.groupTuple())

        use_ref_ann = !ref_annotation.name.startsWith('OPTIONAL_FILE')

        run_gffcompare(merge_gff_bundles.out.gff, ref_annotation)

        if (params.denovo){
            // Use the perd-sample, de novo-assembled CDS
            seq_for_transcriptome_build = m.cds
        }else {
            // If doing reference based assembly, there is only one reference
            // So map this reference to all sample_ids
            seq_for_transcriptome_build = sample_ids.flatten().combine(Channel.fromPath(params.ref_genome))
        }

        makeReport(software_versions,
                    workflow_params,
                    params.denovo,
                    summariseConcatReads.out.summary
                    .join(m.stats)
                    .join(run_gffcompare.out.gffcmp_dir)
                    .join(preprocess_reads.out.report)
                    .join(merge_gff_bundles.out.gff)
                    .toList().transpose().toList()
                    )

        report = makeReport.out.report

       get_transcriptome(merge_gff_bundles.out.gff
            .join(run_gffcompare.out.gffcmp_dir)
            .join(seq_for_transcriptome_build))

       if (use_ref_ann){
            results = preprocess_reads.out.report
                        .concat(run_gffcompare.output.gffcmp_dir,
                        m.stats,
                        get_transcriptome.out.flatMap(map_sample_ids_cls))
                        .map {it -> it[1]}
                        .concat(makeReport.out.report)

       }
        if (!use_ref_ann && !params.denovo){
            results = preprocess_reads.out.report
                        .concat(m.stats,
                        get_transcriptome.out.flatMap(map_sample_ids_cls))
                        .map {it -> it[1]}
                        .concat(makeReport.out.report)

        }
        if (params.denovo){
            results = m.cds
                      .concat(m.stats,
                      seq_for_transcriptome_build,
                      get_transcriptome.out.flatMap(map_sample_ids_cls),
                      merge_gff_bundles.out.gff,
                      m.opt_qual_ch.flatMap {it ->
                      l = []
                       for (x in it[1..-1]){
                            l.add(tuple(it[0], x))
                            }
                        return l
                      })
                      .map {it -> it[1]}
                      .concat(makeReport.out.report)
        }

    emit:
        results
        telemetry = workflow_params
}

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    start_ping()

    fastq = file(params.fastq, type: "file")

     if (!fastq.exists()) {
        println("--fastq: File doesn't exist, check path.")
        exit 1
    }

    if (!params.denovo && !params.ref_genome){
        println("--ref_genome must be supplied unless doing de novo assembly (--denovo)")
        exit 1
    }


    if (params.ref_genome){
        ref_genome = file(params.ref_genome, type: "file")
        if (!ref_genome.exists()) {
            println("--ref_genome: File doesn't exist, check path.")
            exit 1
        }
    }else {
        ref_genome = file("$projectDir/data/OPTIONAL_FILE")
    }

    if (params.denovo && params.ref_annotation) {
        println("Reference annotation with de denovo assembly is not supported")
        exit 1
    }

    if (params.ref_annotation){
        ref_annotation = file(params.ref_annotation, type: "file")
        if (!ref_annotation.exists()) {
            println("--annotation: File doesn't exist, check path.")
            exit 1
        }
    }else{
        ref_annotation = file("$projectDir/data/OPTIONAL_FILE")
    }

    reads = fastq_ingress(
        params.fastq, params.out_dir, params.sample, params.sample_sheet, params.sanitize_fastq
    )

    pipeline(reads, ref_genome, ref_annotation)

    output(pipeline.out.results)

    end_ping(pipeline.out.telemetry)
}
