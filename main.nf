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
include { reference_assembly } from './subworkflows/reference_assembly'
include { denovo_assembly } from './subworkflows/denovo_assembly'
include { gene_fusions } from './subworkflows/JAFFAL/gene_fusions'
include { differential_expression } from './subworkflows/differential_expression'



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
    stringtie --rf ${G_FLAG} -L -v -p ${task.cpus} ${params.stringtie_opts} \
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
        path ("${sample_id}_annotated.gtf"), emit: gtf, optional: true
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
        cp ${out_dir}/str_merged.annotated.gtf ${sample_id}_annotated.gtf
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

process merge_transcriptomes {
    // Merge the transcriptomes from all samples
    label 'isoforms'
    input:
        path "query_annotations/*"
        path ref_annotation
        path ref_genome
    output:
        path "non_redundant.fasta", emit: fasta
        path "stringtie.gtf", emit: gtf
    """
    stringtie --merge -G $ref_annotation -p ${task.cpus} -o stringtie.gtf query_annotations/*
    seqkit subseq --feature "transcript" --gtf-tag "transcript_id" --gtf stringtie.gtf $ref_genome > temp_transcriptome.fasta
    seqkit rmdup -s < temp_transcriptome.fasta > temp_del_repeats.fasta
    cat temp_del_repeats.fasta | sed 's/>.* />/'  | sed -e 's/_[0-9]* \\[/ \\[/' > temp_rm_empty_seq.fasta
    awk 'BEGIN {RS = ">" ; FS = "\\n" ; ORS = ""} \$2 {print ">"\$0}' temp_rm_empty_seq.fasta > non_redundant.fasta
    rm temp_transcriptome.fasta
    rm temp_del_repeats.fasta
    rm temp_rm_empty_seq.fasta
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
        path "de_report/*"
        path "seqkit/*"
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
    if [ -e "de_report/OPTIONAL_FILE" ]; then
        dereport=""
    else
        dereport="--de_report true --de_stats "seqkit/*""
        mv de_report/*.gtf de_report/stringtie_merged.gtf
    fi

    report.py --report $report_name \
    --versions $versions \
    --params params.json \
    $OPT_ALN \
    $OPT_PC_REPORT \
    --sample_ids $sids \
    --summaries $seq_summaries \
    --gffcompare_dir $gffcmp_dir \
    --gff_annotation $gff_annotation \
    --isoform_table_nrows $params.isoform_table_nrows \
    $OPT_JAFFAL_CSV \
    $OPT_DENOVO \
    \$dereport 

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
        condition_sheet
        ref_transcriptome
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
             assembly = denovo_assembly(full_len_reads, ref_genome)

        } else {
            build_minimap_index(ref_genome)
            println("Doing reference based transcript analysis")
            assembly = reference_assembly(build_minimap_index.out.index, ref_genome, full_len_reads)
        }

        split_bam(assembly.bam)

        assemble_transcripts(split_bam.out.bundles.flatMap(map_sample_ids_cls), ref_annotation)

        merge_gff_bundles(assemble_transcripts.out.gff_bundles.groupTuple())

        use_ref_ann = !ref_annotation.name.startsWith('OPTIONAL_FILE')

        run_gffcompare(merge_gff_bundles.out.gff, ref_annotation)

        if (params.denovo){
            // Use the per-sample, de novo-assembled CDS
            seq_for_transcriptome_build = assembly.cds
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

        get_transcriptome(
            merge_gff_bundles.out.gff
            .join(run_gffcompare.out.gffcmp_dir)
            .join(seq_for_transcriptome_build))

        if (params.de_analysis){

            if (!params.ref_transcriptome){
                merge_transcriptomes(run_gffcompare.output.gtf.collect(), ref_annotation, ref_genome)
                transcriptome = merge_transcriptomes.out.fasta
                gtf = merge_transcriptomes.out.gtf
            }
            else {
                transcriptome = ref_transcriptome
                gtf = Channel.fromPath(ref_annotation)
            }
            check_match = Channel.fromPath(params.condition_sheet)
            check_condition_sheet = check_match.splitCsv(header: true).map{ row -> tuple(
                row.sample_id)
            }
            check_condition_sheet.join(summariseConcatReads.out.input_reads, failOnMismatch: true)
            de = differential_expression(transcriptome, summariseConcatReads.out.input_reads, condition_sheet, gtf)
            de_report = de.all_de
            count_transcripts_file = de.count_transcripts
            dtu_plots = de.dtu_plots
            de_outputs = de.de_outputs
        } else{
            de_report = file("$projectDir/data/OPTIONAL_FILE")
            count_transcripts_file = file("$projectDir/data/OPTIONAL_FILE")
        }
        makeReport(
            software_versions,
            workflow_params,
            params.denovo,
            pychopper_report,
            jaffal_out,
            summariseConcatReads.out.summary
            .join(assembly.stats)
            .join(run_gffcompare.out.gffcmp_dir)
            .join(merge_gff_bundles.out.gff)
            .toList().transpose().toList(),
            de_report,
            count_transcripts_file)

        report = makeReport.out.report



       if (use_ref_ann){
            results = run_gffcompare.output.gffcmp_dir
                      .concat(
                      assembly.stats,
                      get_transcriptome.out.transcriptome.flatMap(map_sample_ids_cls))
                      .map {it -> it[1]}
                      .concat(makeReport.out.report)

       }
        if (!use_ref_ann && !params.denovo){
            results =  assembly.stats
                      .concat(
                       get_transcriptome.out.transcriptome.flatMap(map_sample_ids_cls))
                      .map {it -> it[1]}
                      .concat(makeReport.out.report)

        }
        if (params.denovo){
            results = assembly.cds
                      .concat(assembly.stats,
                      seq_for_transcriptome_build,
                      get_transcriptome.out.transcriptome.flatMap(map_sample_ids_cls),
                      merge_gff_bundles.out.gff,
                      assembly.opt_qual_ch.flatMap {
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

        if (params.de_analysis){
            results = results.concat(de.dtu_plots, de_outputs)
        }

    emit:
        results
        telemetry = workflow_params
}

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    if (params.disable_ping == false) {
            Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }

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
    ref_transcriptome = file("$projectDir/data/OPTIONAL_FILE")
    if (params.ref_transcriptome){
        ref_transcriptome = file(params.ref_transcriptome, type:"file")
    }
    if (params.de_analysis){
        if (!params.ref_annotation){
            error = "You must provide a reference annotation."
        }
        if (!params.condition_sheet){

            error = "You must provide a condition_sheet or set de_analysis to false."
        }
        condition_sheet = file(params.condition_sheet, type:"file")
    } else{
        condition_sheet = file("$projectDir/data/OPTIONAL_FILE")
    }
    if (error){
        println(error)
    }else{
        reads = fastq_ingress([
            "input":params.fastq,
            "sample":params.sample,
            "sample_sheet":params.sample_sheet])

        pipeline(reads, ref_genome, ref_annotation,
            jaffal_refBase, params.jaffal_genome, params.jaffal_annotation,
            condition_sheet, ref_transcriptome)

        output(pipeline.out.results)
    }
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }
    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
