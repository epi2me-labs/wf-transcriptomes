#!/usr/bin/env nextflow

// Developer notes
// 
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion.
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
import java.util.ArrayList;
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress' 
include { start_ping; end_ping } from './lib/ping'


process summariseReads {
    // concatenate fastq and fastq.gz in a dir

    label "isoforms"
    cpus 1
    input:
        tuple path(directory), val(sample_name), val(type)
    output:
        tuple val(sample_name), path('*'), emit: summary
    script:
    """
    fastcat -s ${sample_name} -r ${sample_name}.stats -x ${directory} > /dev/null
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
    python -c "import seaborn; print(f'seaborn,{seaborn.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    bedtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    python -c "import pychopper; print(f'pychopper,{pychopper.__version__}')" >> versions.txt
    gffread --version | sed 's/^/gffread,/' >> versions.txt
    seqkit version | head -n 1 | sed 's/ /,/' >> versions.txt
    csvtk version | head -n 1 | sed 's/ /,/' >> versions.txt
    stringtie --version | sed 's/^/stringtie,/' >> versions.txt
    gffcompare --version | head -n 1 | sed 's/ /,/' >> versions.txt
    spoa --version | sed 's/^/spoa,/' >> versions.txt
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
    cpus params.threads

    input:
        tuple path(directory), val(sample_id), val(type) 
    output:
         tuple val(sample_id), path("full_length_reads.fq"), emit: full_len_reads
         tuple val(sample_id), path('cdna_classifier_report.tsv'),  emit: report
//          val "${sample_id}", emit: sample_id
//          path 'cdna_classifier_report.tsv', optional: true, emit: cdna_class_report
         // Not sure if this is an antipattern, but it's function is to publish all the files to the publishDir dir
    script:
        """
        fastcat -s ${sample_id} -r ${sample_id}.stats -x ${directory} > input_reads.fq
    
        if [[ ${params.use_pychopper} == true ]];
        then
            cdna_classifier.py -t ${params.threads} ${params.pychopper_opts} input_reads.fq full_length_reads.fq
            generate_pychopper_stats.py --data cdna_classifier_report.tsv --output .
        else
            ln -s `realpath input_reads.fq` full_length_reads.fq
        fi
        """
}

process generate_fq_stats {
    label "isoforms"

    input:
        tuple(sample_id), path(fastq), val(unused)

    output:
        path "*"
    script:
    """
        run_fastq_qc.py --fastq ${fastq} --output .
    """
}

process build_minimap_index{
    /*
    Build minimap index from reference genome
    */
    label "isoforms"
    cpus params.threads

    input:
        file reference
    output:
        path "genome_index.mmi", emit: index
    script:
    """
        minimap2 -t ${params.threads} ${params.minimap_index_opts} -I 1000G -d "genome_index.mmi" ${reference}
    """
}


process map_reads{
    /*
    Map reads to reference using minimap2.
    Filter reads by mapping quality.
    Filter reads where length of poly(A) > max_poly_run at either ends of the read (defined by poly_context)
    */
    label "isoforms"
    cpus params.threads

    input:
       file index
       file reference
       tuple val(sample_id), file(fastq_reads)

    output:
       tuple val(sample_id), path("reads_aln_sorted.bam"), emit: bam
       tuple val(sample_id), path("read_aln_stats.tsv"), emit: stats
    script:
        def ab = "reads_aln_sorted.bam"
        def af = "internal_priming_fail.tsv"
        def fs = "context_internal_priming_fail_start.fasta"
        def fe = "context_internal_priming_fail_end.fasta"
        def fasta_reads = "reads.fa"
        def ContextFilter = """AlnContext: { Ref: "${reference}", LeftShift: -${params.poly_context}, RightShift: ${params.poly_context},
           RegexEnd: "[Aa]{${params.max_poly_run},}",
           Stranded: True,Invert: True, Tsv: "internal_priming_fail.tsv"} """
       """
       seqkit fq2fa ${fastq_reads} -o ${fasta_reads};
       minimap2 -t ${params.threads} -ax splice ${params.minimap2_opts} ${index} ${fasta_reads}\
       | samtools view -q ${params.minimum_mapping_quality} -F 2304 -Sb -\
       | seqkit bam -j ${params.threads} -x -T '${ContextFilter}' -\
       | samtools sort -@ ${params.threads} -o ${ab} -;

       ((seqkit bam -s -j ${params.threads} ${ab} 2>&1)  | tee read_aln_stats.tsv ) || true

       if [[ -s ${af} ]];
       then
           tail -n +2 ${af} | awk '{{print ">" \$1 "\\n" \$4 }}' - > ${fs}
           tail -n +2 ${af} | awk '{{print ">" \$1 "\\n" \$6 }}' - > ${fe}
       fi
       """
}

process plot_aln_stats{
    /*
    Create a pdf of alignemnt statistis. 
    */
    label 'isoforms'

    input:
        tuple val(sample_id), path (aln_stats)
    output:
        path "*"
    script:
    """
    plot_aln_stats.py ${aln_stats} -r read_aln_stats.pdf
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
        """
    else
        """
        mkdir -p ./${sample_id}_bam_bundles
        ln -s ${bam} ${sample_id}_bam_bundles/000000000_ALL:0:1_bundle.bam
        """
}


process stringtie{
    /*
    Takes in aligned reads in bam format that may be a chunk of a larger alignment file.
    $G_FLAG specifies whether or not to use reference annotation as a guide in transcript assembly. 

    Output gff annotation files in a tuple with `sample_id` for combining into samples late rin the pipeline.
    */
    label 'isoforms'
    cpus params.threads

    input:
        tuple val(sample_id), path(bam)
        file ref_annotation
    output:
        tuple val(sample_id), path('*.gff'), emit: gff_bundles
    script:
        def out_filename = bam.name.replaceFirst(~/\.[^\.]+$/, '') + '.gff'
        def label = "STR.${bam.name.split('_')[0].toInteger()}."
        def G_FLAG = ref_annotation.name.startsWith('OPTIONAL_FILE') ? '' : "-G ${ref_annotation}"
    """
    stringtie --rf ${G_FLAG} -L -v -p ${params.threads} ${params.stringtie_opts} -o  ${out_filename} \
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
    def merged_gff = "str_merged_${sample_id}.gff"
    """
    echo '#gff-version 2' >> $merged_gff;
    echo '#pipeline-nanopore-isoforms: stringtie' >> $merged_gff;

    for fn in ${gff_bundle};
    do
        grep -v '#' \$fn >> $merged_gff

    done
    """
}

process run_gff_compare{
    /*
    Compare query and reference annotations.
    */

    label 'isoforms'

    input:
       tuple val(sample_id), path(query_annotation)
       path ref_annotation
    output:
        tuple val(sample_id), path('str_merged.annotated.gtf'), emit: merged_annotated
        tuple val(sample_id), path('str_merged.stats'), emit: stats
        tuple val(sample_id), path('str_merged.tracking'), emit: tracking
    script:
        """
        echo "Doing comparison of reference annotation: ${ref_annotation} and the current annotation"
        gffcompare -o str_merged -r ${ref_annotation} ${params.gffcompare_opts} ${query_annotation}
        generate_tracking_summary.py --tracking str_merged.tracking --output_dir . --annotation ${ref_annotation}

        if [[ ${params.plot_gffcmp_stats} == true ]];
        then
            plot_gffcmp_stats.py -r str_gffcmp_report.pdf -t str_merged.tracking str_merged.stats;
        fi
        """
}


process run_gffread{
        /*
        Write out a transctiptome file based on the gff annotations.
        */
        label 'isoforms'

        input:
            tuple val(sample_id), path(gff_merged), path(merged_ann_gff)
            path reference_seq
        output:
            tuple val(sample_id), path('*.fas'), emit: transcriptome
        script:
        def str_transcriptome = "${sample_id}_str_transcriptome.fas"
        def merged_transcriptome = "${sample_id}_merged_transcriptome.fas"
        """

        gffread -g ${reference_seq} -w ${str_transcriptome} ${gff_merged}
        if [ -f ${merged_ann_gff} ]
        then
            gffread -F -g ${reference_seq} -w ${merged_transcriptome} ${merged_ann_gff}
        else
            touch ${merged_transcriptome}
        fi
        """
}

process makeReport {

    label "isoforms"

    input:
        path versions
        path params
        tuple val(sample_id), path(seqs), path(aln_stats), path(gffcompare_tracking), path(gffcompare_stats),
                path(cdna_class_report)
    output:
        tuple val(sample_id), path("wf-isoforms-*.html"), emit: report
    script:
        def report_name = "wf-isoforms-${sample_id}_report.html"
        def opt_gff_track = gffcompare_tracking.name.startsWith('OPTIONAL_FILE') ? '' : "--gffcompare_tracking ${gffcompare_tracking}"
        def opt_gff_stats = gffcompare_stats.name.startsWith('OPTIONAL_FILE') ? '' : "--gffcompare_stats ${gffcompare_stats}"
    """
    report.py ${report_name} --versions ${versions} ${seqs} --params params.json \
    --alignment_stats ${aln_stats} \
    ${opt_gff_track} \
    ${opt_gff_stats} \
    --pychop_report ${cdna_class_report}
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "isoforms"
    publishDir "${params.results_dir}/${sample_id}", mode: 'copy', pattern: "*"

    input:
        tuple val(sample_id), file(fname)
    output:
        file fname
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

        summariseReads(reads)
        sample_ids = summariseReads.out.summary.collect({it -> it[0]})
        software_versions = getVersions()
        workflow_params = getParams()
        preprocess_reads(reads)
//         generate_fq_stats(preprocess_reads.out.sample)  skip for now

        build_minimap_index(ref_genome)
        map_reads(build_minimap_index.out.index, ref_genome, preprocess_reads.out.full_len_reads)
//         plot_aln_stats(map_reads.out.bam)
        split_bam(map_reads.out.bam)
        stringtie(split_bam.out.bundles.flatMap(map_sample_ids_cls), ref_annotation)
        merge_gff_bundles(stringtie.out.gff_bundles.groupTuple())

        use_ref_ann = !ref_annotation.name.startsWith('OPTIONAL_FILE')

        if (use_ref_ann){
            run_gff_compare(merge_gff_bundles.out.gff, ref_annotation)
            run_gffread(merge_gff_bundles.out.gff.join(run_gff_compare.out.merged_annotated),
                        ref_genome)
            gff_tracking = run_gff_compare.out.tracking
            gff_stats = run_gff_compare.out.stats
        }else{
            // Create dummy file paths to satisfy required path inputs of processes
            // Add to tuple with sample_id to guide to correct sample
            gff_tracking = sample_ids.combine(Channel.fromPath("$projectDir/data/OPTIONAL_FILE")).view()
            gff_stats = sample_ids.combine(Channel.fromPath("$projectDir/data/OPTIONAL_FILE_1"))
        }

        makeReport(
                    software_versions,
                    workflow_params,
                    summariseReads.out.summary
                    .join(map_reads.out.stats)
                    .join(gff_tracking)
                    .join(gff_stats)
                    .join(preprocess_reads.out.report)
                    )
        if (use_ref_ann){
            results = preprocess_reads.out.report
                        .concat(merge_gff_bundles.out.gff,
                        run_gff_compare.out.stats,
                        run_gff_compare.out.merged_annotated,
                        run_gffread.out.transcriptome,
                        makeReport.out.report
                        )
        }else{  // Write a minimal report if no reference annotation is given
            results = preprocess_reads.out.report
                        .concat(merge_gff_bundles.out.gff,
                                makeReport.out.report
                                )
            }
    emit:
        results
        telemetry = workflow_params
}

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    start_ping()
    params.results_dir = "${params.out_dir}/output"

    fastq = file(params.fastq, type: "file")

     if (!fastq.exists()) {
        println("--fastq: File doesn't exist, check path.")
        exit 1
    }

    if (params.ref_genome){
        ref_genome = file(params.ref_genome, type: "file")
        if (!ref_genome.exists()) {
            println("--reference: File doesn't exist, check path.")
            exit 1
        }
    }

    ref_annotation = null
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

    output(
        pipeline.out.results
        )

    end_ping(pipeline.out.telemetry)
}
