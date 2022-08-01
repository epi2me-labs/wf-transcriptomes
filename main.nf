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
include { reference_assembly } from './subworkflows/reference_assembly'
include { denovo_assembly } from './subworkflows/denovo_assembly'
include { gene_fusions } from './subworkflows/JAFFAL/gene_fusions'


process summariseConcatReads {
    // concatenate fastq and fastq.gz in a dir write stats

    label "isoforms"
    cpus 1
    input:
        tuple path(directory), val(meta)
    output:
        tuple val(meta.sample_id), path("${meta.sample_id}.fastq"), emit: input_reads
        tuple val(meta.sample_id), path('*.stats'), emit: summary
    script:
    """
    fastcat -s ${meta.sample_id} -r ${meta.sample_id}.stats -x ${directory} >  ${meta.sample_id}.fastq
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
    python -c "import sklearn; print(f'scikit-learn,{sklearn.__version__}')" >> versions.txt
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
    Optionally classify, trim, and orient cDNA reads using pychopper
    */

    label "isoforms"
    cpus 4

    input:
        tuple val(sample_id), path(input_reads)
    output:
         tuple val(sample_id), path("${sample_id}_full_length_reads.fastq"), emit: full_len_reads
         path '*.tsv',  emit: report
    script:
        """
        pychopper -t ${params.threads} ${params.pychopper_opts} ${input_reads} ${sample_id}_full_length_reads.fastq
        mv pychopper.tsv ${sample_id}_pychopper.tsv
        generate_pychopper_stats.py --data ${sample_id}_pychopper.tsv --output .

        # Add sample id column
        sed "1s/\$/\tsample_id/; 1 ! s/\$/\t${sample_id}/" ${sample_id}_pychopper.tsv > tmp
        mv tmp ${sample_id}_pychopper.tsv
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
    """
    n=`samtools view -c $bam`
    if [[ n -lt 1 ]]
    then
        echo 'There are no reads mapping for $sample_id. Exiting!'
        exit 1
    fi

    re='^[0-9]+\$'

    if [[ $params.bundle_min_reads =~ \$re ]]
    then
        echo "Bundling up the bams"
        seqkit bam -j ${params.threads} -N ${params.bundle_min_reads} ${bam} -o  bam_bundles/
        let i=1
        for b in bam_bundles/*.bam; do
            echo \$b
            newname="${sample_id}_batch_\${i}.bam"
            mv \$b \$newname
           ((i++))
        done
    else
        echo 'no bundling'
        ln -s ${bam} ${sample_id}_batch_1.bam
    fi
    """
}


process assemble_transcripts{
    /*
    Assemble transcripts using stringtie.
    Take aligned reads in bam format that may be a chunk of a larger alignment file.
    Optionally use reference annotation to guide assembly.

    Output gff annotation files in a tuple with `sample_id` for combining into samples later in the pipeline.
    */
    label 'isoforms'
    cpus params.threads

    input:
        tuple val(sample_id), path(bam)
        path ref_annotation
    output:
        tuple val(sample_id), path('*.gff'), emit: gff_bundles
    script:
        def G_FLAG = ref_annotation.name.startsWith('OPTIONAL_FILE') ? '' : "-G ${ref_annotation}"
        def prefix =  bam.name.split(/\./)[0]

    """
    stringtie --rf ${G_FLAG} -L -v -p ${params.threads} ${params.stringtie_opts} \
    -o  ${prefix}.gff -l ${prefix} ${bam} 2>/dev/null
     """
}


process merge_gff_bundles{
    /*
    Merge gff bundles into a single gff file per sample.
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
                gffread -F -g ${reference_seq} -w ${merged_transcriptome} $gffcmp_dir/str_merged.annotated.gtf
        fi
        """
}

process makeReport {

    label "isoforms"

    input:
        path versions
        path "params.json"
        val denovo
        path pychopper_report
        path jaffal_csv
        tuple val(sample_ids),
              path(seq_summaries),
              path(aln_stats),
              path(gffcmp_dir),
              path(gff_annotation)
    output:
        path("wf-transcriptomes-*.html"), emit: report
    script:
        // Convert the sample_id arrayList.
        sids = new BlankSeparatedList(sample_ids)
        def report_name = "wf-transcriptomes-report.html"
        def OPT_ALN = denovo ?  '' : "--alignment_stats ${aln_stats}"
        def OPT_DENOVO = denovo ? "--denovo" : ''
        def OPT_PC_REPORT = pychopper_report.name.startsWith('OPTIONAL_FILE') ? '' : "--pychop_report ${pychopper_report}"
        def OPT_JAFFAL_CSV = jaffal_csv.name.startsWith('OPTIONAL_FILE') ? '' : "--jaffal_csv ${jaffal_csv}"

    """
    report.py --report $report_name \
    --versions $versions \
    --params params.json \
    $OPT_ALN \
    $OPT_PC_REPORT \
    --sample_ids $sids \
    --summaries $seq_summaries \
    --gffcompare_dir $gffcmp_dir \
    --gff_annotation $gff_annotation \
    --transcript_table_cov_thresh $params.transcript_table_cov_thresh \
    $OPT_JAFFAL_CSV \
    $OPT_DENOVO
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    label "isoforms"
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
        jaffal_refBase
        jaffal_genome
        jaffal_annotation
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

        if (!params.direct_rna){
            preprocess_reads(summariseConcatReads.out.input_reads)
            full_len_reads = preprocess_reads.out.full_len_reads
            pychopper_report = preprocess_reads.out.report.collectFile(keepHeader: true)
        }
        else{
            full_len_reads = summariseConcatReads.out.input_reads
            pychopper_report = file("$projectDir/data/OPTIONAL_FILE")
        }

        if (params.denovo){
            println("Doing de novo assembly")
             m = denovo_assembly(full_len_reads, ref_genome)

        } else {
            build_minimap_index(ref_genome)
            println("Doing reference based transcript analysis")
            m = reference_assembly(build_minimap_index.out.index, ref_genome, full_len_reads)
        }

        split_bam(m.bam)

        assemble_transcripts(split_bam.out.bundles.flatMap(map_sample_ids_cls), ref_annotation)

        merge_gff_bundles(assemble_transcripts.out.gff_bundles.groupTuple())

        use_ref_ann = !ref_annotation.name.startsWith('OPTIONAL_FILE')

        run_gffcompare(merge_gff_bundles.out.gff, ref_annotation)

        if (params.denovo){
            // Use the per-sample, de novo-assembled CDS
            seq_for_transcriptome_build = m.cds
        }else {
            // For reference based assembly, there is only one reference
            // So map this reference to all sample_ids
            seq_for_transcriptome_build = sample_ids.flatten().combine(Channel.fromPath(params.ref_genome))
        }

        if (jaffal_refBase){
            gene_fusions(full_len_reads, jaffal_refBase, jaffal_genome, jaffal_annotation)
            jaffal_out = gene_fusions.out.results_csv.collectFile(keepHeader: true, name: 'jaffal.csv')
        }else{
            jaffal_out = file("$projectDir/data/OPTIONAL_FILE_1")
        }

        makeReport(
            software_versions,
            workflow_params,
            params.denovo,
            pychopper_report,
            jaffal_out,
            summariseConcatReads.out.summary
            .join(m.stats)
            .join(run_gffcompare.out.gffcmp_dir)
            .join(merge_gff_bundles.out.gff)
            .toList().transpose().toList())

        report = makeReport.out.report

        get_transcriptome(
            merge_gff_bundles.out.gff
            .join(run_gffcompare.out.gffcmp_dir)
            .join(seq_for_transcriptome_build))

       if (use_ref_ann){
            results = run_gffcompare.output.gffcmp_dir
                      .concat(
                      m.stats,
                      get_transcriptome.out.flatMap(map_sample_ids_cls))
                      .map {it -> it[1]}
                      .concat(makeReport.out.report)

       }
        if (!use_ref_ann && !params.denovo){
            results =  m.stats
                      .concat(
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
                      m.opt_qual_ch.flatMap {
                          it ->
                          l = []
                          for (x in it[1..-1]){
                              l.add(tuple(it[0], x))
                          }
                        return l
                      })
                      .map {it -> it[1]}
                      .concat(makeReport.out.report)
        }
        if (params.jaffal_refBase){
            results = results
                .concat(gene_fusions.out.results
                .map {it -> it[1]})
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

    error = null

    if (!fastq.exists()) {
        error = "--fastq: File doesn't exist, check path."
    }

    if (!params.denovo && !params.ref_genome){
        error = "--ref_genome must be supplied unless doing de novo assembly (--denovo)"
    }

    if (params.ref_genome){
        ref_genome = file(params.ref_genome, type: "file")
        if (!ref_genome.exists()) {
            error = "--ref_genome: File doesn't exist, check path."
        }
    }else {
        ref_genome = file("$projectDir/data/OPTIONAL_FILE")
    }

    if (params.denovo && params.ref_annotation) {
        error = "Reference annotation with de denovo assembly is not supported"
    }

    if (params.ref_annotation){
        ref_annotation = file(params.ref_annotation, type: "file")
        if (!ref_annotation.exists()) {
            error = "--annotation: File doesn't exist, check path."
        }
    }else{
        ref_annotation = file("$projectDir/data/OPTIONAL_FILE")
    }
     if (params.jaffal_refBase){
        jaffal_refBase = file(params.jaffal_refBase, type: "dir")
        if (!jaffal_refBase.exists()) {
            error = "--jaffa_refBase: Directory doesn't exist, check path."
        }
     }else{
        jaffal_refBase = null
     }

    if (error){
        println(error)
    }else{
        reads = fastq_ingress([
            "input":params.fastq,
            "sample":params.sample,
            "sample_sheet":params.sample_sheet,
            "sanitize": params.sanitize_fastq,
            "output":params.out_dir])

        pipeline(reads, ref_genome, ref_annotation, jaffal_refBase, params.jaffal_genome, params.jaffal_annotation)

        output(pipeline.out.results)

        end_ping(pipeline.out.telemetry)
    }
}
