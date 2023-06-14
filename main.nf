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

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

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
#     isONclust2 version | sed 's/ version: /,/' >> versions.txt
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


process decompress_ref {
    label "isoforms"
    cpus 1
    input:
        path compressed_ref
    output:
        path "${compressed_ref.baseName}", emit: decompressed_ref
    """
    gzip -df ${compressed_ref}
    """
}


process decompress_annotation {
    label "isoforms"
    cpus 1
    input:
        path compressed_annotation
    output:
        path "${compressed_annotation.baseName}"
    """
    gzip -df ${compressed_annotation}
    """
}

// Remove empty transcript ID fields
process preprocess_ref_annotation {
    label "isoforms"
    cpus 1
    input:
        path ref_annotation
    output:
        path "ammended.${ref_annotation}"
    """
    sed -i -e 's/transcript_id "";//g' ${ref_annotation}
    mv ${ref_annotation} "ammended.${ref_annotation}"
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
        tuple val(meta), path(input_reads)
    output:
         tuple val("${meta.alias}"), path("${meta.alias}_full_length_reads.fastq"), emit: full_len_reads
         path '*.tsv',  emit: report
    script:
        """
        pychopper -t ${params.threads} ${params.pychopper_opts} ${input_reads} ${meta.alias}_full_length_reads.fastq
        mv pychopper.tsv ${meta.alias}_pychopper.tsv
        workflow-glue generate_pychopper_stats --data ${meta.alias}_pychopper.tsv --output .

        # Add sample id column
        sed "1s/\$/\tsample_id/; 1 ! s/\$/\t${meta.alias}/" ${meta.alias}_pychopper.tsv > tmp
        mv tmp ${meta.alias}_pychopper.tsv
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
        tuple val(sample_id), path(bam), path(ref_annotation)
        val use_ref_ann
    output:
        tuple val(sample_id), path('*.gff'), emit: gff_bundles
    script:
        def G_FLAG = use_ref_ann == false ? '' : "-G ${ref_annotation}"
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

    if (params.transcriptome_source == "denovo"){
        """
        mkdir $out_dir
        """
    } else {
        """
        mkdir $out_dir
        echo "Doing comparison of reference annotation: ${ref_annotation} and the query annotation"

        gffcompare -o ${out_dir}/str_merged -r ${ref_annotation} \
            ${params.gffcompare_opts} ${query_annotation}

        workflow-glue generate_tracking_summary --tracking $out_dir/str_merged.tracking \
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
        path "pychopper_report/*"
        path"jaffal_csv/*"
        val sample_ids
        path per_read_stats
        path "aln_stats/*"
        path gffcmp_dir
        path "gff_annotation/*"
        path "de_report/*"
        path "seqkit/*"
    output:
        path("wf-transcriptomes-*.html"), emit: report
    script:
        // Convert the sample_id arrayList.
        sids = new BlankSeparatedList(sample_ids)
        def report_name = "wf-transcriptomes-report.html"
        def OPT_DENOVO = params.transcriptome_source == "denovo" ? "--denovo" : ''
    """
    if [ -f "de_report/OPTIONAL_FILE" ]; then
        dereport=""
    else
        dereport="--de_report true --de_stats "seqkit/*""
        mv de_report/*.g*f* de_report/stringtie_merged.gtf
    fi
    if [ -f "gff_annotation/OPTIONAL_FILE" ]; then
        OPT_GFF=""
    else
        OPT_GFF="--gffcompare_dir ${gffcmp_dir} --gff_annotation gff_annotation/*"
        
    fi
    if [ -f "jaffal_csv/OPTIONAL_FILE" ]; then
        OPT_JAFFAL_CSV=""
    else
        OPT_JAFFAL_CSV="--jaffal_csv jaffal_csv/*"
    fi
    if [ -f "aln_stats/OPTIONAL_FILE" ]; then
        OPT_ALN=""
    else
        OPT_ALN="--alignment_stats aln_stats/*"
    fi
    if [ -f "pychopper_report/OPTIONAL_FILE" ]; then
        OPT_PC_REPORT=""
    else
        OPT_PC_REPORT="--pychop_report pychopper_report/*"
    fi
    workflow-glue report --report $report_name \
    --versions $versions \
    --params params.json \
    \$OPT_ALN \
    \$OPT_PC_REPORT \
    --sample_ids $sids \
    --stats $per_read_stats \
    \$OPT_GFF \
    --isoform_table_nrows $params.isoform_table_nrows \
    \$OPT_JAFFAL_CSV \
    $OPT_DENOVO \
    \$dereport 

    """
}


// Creates a new directory named after the sample alias and moves the fastcat results
// into it.
process collectFastqIngressResultsInDir {
    label "isoforms"
    input:
        // both the fastcat seqs as well as stats might be `OPTIONAL_FILE` --> stage in
        // different sub-directories to avoid name collisions
        tuple val(meta), path(concat_seqs, stageAs: "seqs/*"), path(fastcat_stats,
            stageAs: "stats/*")
    output:
        // use sub-dir to avoid name clashes (in the unlikely event of a sample alias
        // being `seq` or `stats`)
        path "out/*"
    script:
    String outdir = "out/${meta["alias"]}"
    String metaJson = new JsonBuilder(meta).toPrettyString()
    String concat_seqs = \
        (concat_seqs.fileName.name == OPTIONAL_FILE.name) ? "" : concat_seqs
    String fastcat_stats = \
        (fastcat_stats.fileName.name == OPTIONAL_FILE.name) ? "" : fastcat_stats
    """
    mkdir -p $outdir
    echo '$metaJson' > metamap.json
    mv metamap.json $concat_seqs $fastcat_stats $outdir
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process output {
    // publish inputs to output directory
    label "isoforms"
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
        use_ref_ann
    main:
        if (params.ref_genome && file(params.ref_genome).extension == "gz") {
            // gzipped ref not supported by some downstream tools
            // easier to just decompress and pass it around.
            ref_genome = decompress_ref(ref_genome)
        }else {
            ref_genome = Channel.fromPath(ref_genome)
        }
        if (params.ref_annotation && file(params.ref_annotation).extension == "gz") {
            // gzipped ref not supported by some downstream tools
            // easier to just decompress and pass it around.
            decompress_annot= decompress_annotation(ref_annotation)
            ref_annotation = preprocess_ref_annotation(decompress_annot)
        }else {
            ref_annotation = preprocess_ref_annotation(ref_annotation)
        }
        
        fastq_ingress_results = reads
        // replace `null` with path to optional file
        | map { [ it[0], it[1] ?: OPTIONAL_FILE, it[2] ?: OPTIONAL_FILE ] }
        | collectFastqIngressResultsInDir
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

        

        software_versions = getVersions()
        workflow_params = getParams()
        input_reads = reads.map{ meta, samples, stats -> [meta, samples]}
        sample_ids = input_reads.flatMap({meta,samples -> meta.alias})
        stats = reads.map {
            it[2] ? it[2].resolve('per-read-stats.tsv') : null
        }

        if (!params.direct_rna){
            preprocess_reads(input_reads)
            full_len_reads = preprocess_reads.out.full_len_reads
            pychopper_report = preprocess_reads.out.report.collectFile(keepHeader: true)
        }
        else{
            full_len_reads = input_reads.map{ meta, reads -> [meta.alias, reads]}
            pychopper_report = file("$projectDir/data/OPTIONAL_FILE")
        }
        if (params.transcriptome_source != "precomputed"){
        
            if (params.transcriptome_source == "denovo"){
                log.info("Doing de novo assembly")
                log.info("WARNING: The `--transcriptome_source` denovo option may have unexpected results and errors. If possible it is preferable to use the reference-guided pipeline.")
                assembly = denovo_assembly(full_len_reads, ref_genome)

            } else {
                build_minimap_index(ref_genome)
                log.info("Doing reference based transcript analysis")
                assembly = reference_assembly(build_minimap_index.out.index, ref_genome, full_len_reads)
            }
            assembly_stats = assembly.stats.map{ it -> it[1]}.collect()
     
            split_bam(assembly.bam)

            assemble_transcripts(split_bam.out.bundles.flatMap(map_sample_ids_cls).combine(ref_annotation),use_ref_ann)

            merge_gff_bundles(assemble_transcripts.out.gff_bundles.groupTuple())
            run_gffcompare(merge_gff_bundles.out.gff, ref_annotation)

            if (params.transcriptome_source == "denovo"){
                // Use the per-sample, de novo-assembled CDS
                seq_for_transcriptome_build = assembly.cds
            }else {
                // For reference based assembly, there is only one reference
                // So map this reference to all sample_ids
                seq_for_transcriptome_build = sample_ids.flatten().combine(ref_genome)
            }
            get_transcriptome(
                merge_gff_bundles.out.gff
                .join(run_gffcompare.out.gffcmp_dir)
                .join(seq_for_transcriptome_build))

            gff_compare = run_gffcompare.out.gffcmp_dir.map{ it -> it[1]}.collect()
            merge_gff = merge_gff_bundles.out.gff.map{ it -> it[1]}.collect()
            results = Channel.empty()
        }else
        {
            gff_compare = file("$projectDir/data/OPTIONAL_FILE")
            merge_gff = file("$projectDir/data/OPTIONAL_FILE")
            assembly_stats = file("$projectDir/data/OPTIONAL_FILE")
            use_ref_ann = false
            results = Channel.empty()
        }
        if (jaffal_refBase){
                gene_fusions(full_len_reads, jaffal_refBase, jaffal_genome, jaffal_annotation)
                jaffal_out = gene_fusions.out.results_csv.collectFile(keepHeader: true, name: 'jaffal.csv')
            }else{
                jaffal_out = file("$projectDir/data/OPTIONAL_FILE")
        }


        if (params.de_analysis){

            if (!params.ref_transcriptome){
                merge_transcriptomes(run_gffcompare.output.gtf.collect(), ref_annotation, ref_genome)
                transcriptome = merge_transcriptomes.out.fasta
                gtf = merge_transcriptomes.out.gtf
            }
            else {
                transcriptome = ref_transcriptome
                gtf = ref_annotation
            }
            check_match = Channel.fromPath(params.condition_sheet)
            check_condition_sheet = check_match.splitCsv(header: true).map{ row -> tuple(
                row.sample_id)
            }
            join_reads = input_reads.map{ meta, reads -> [meta.alias, reads]}
            check_condition_sheet.join(join_reads, failOnMismatch: true)
            de = differential_expression(transcriptome, input_reads, condition_sheet, gtf)
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
            pychopper_report,
            jaffal_out,
            input_reads.map{ meta, fastq -> meta.alias}.collect(),
            stats,
            assembly_stats,
            gff_compare,
            merge_gff,
            de_report,
            count_transcripts_file)

        report = makeReport.out.report
       
        results = results.concat(makeReport.out.report)
       
       if (use_ref_ann){
            results = run_gffcompare.output.gffcmp_dir.concat(
                      assembly.stats,
                      get_transcriptome.out.transcriptome.flatMap(map_sample_ids_cls))
                      .map {it -> it[1]}
                      .concat(results)

       }
    
        if (!use_ref_ann && params.transcriptome_source == "reference-guided"){
            results =  assembly.stats.concat(
                       get_transcriptome.out.transcriptome.flatMap(map_sample_ids_cls))
                      .map {it -> it[1]}
                      .concat(results)

        }
        if (params.transcriptome_source == "denovo"){
            results = assembly.cds.concat(
                       assembly.stats,
                      seq_for_transcriptome_build,
                      get_transcriptome.out.transcriptome.flatMap(map_sample_ids_cls),
                      assembly.opt_qual_ch.flatMap {
                          it ->
                          l = []
                          for (x in it[1..-1]){
                              l.add(tuple(it[0], x))
                          }
                        return l
                      })
                      .map {it -> it[1]}
                      .concat(results)
        }
        if (params.jaffal_refBase){
            results = results
                .concat(gene_fusions.out.results
                .map {it -> it[1]})
        }

        if (params.de_analysis){
            results = results.concat(de.dtu_plots, de_outputs)
        }

        results  = fastq_ingress_results.map { [it, "fastq_ingress_results"] }.concat(results.map{ [it, null]})
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

    if (params.transcriptome_source == "precomputed" && !params.ref_transcriptome){
        error = "As transcriptome source parameter is precomputed you must include a ref_transcriptome parameter"
    }
    if (params.transcriptome_source == "reference-guided" && !params.ref_genome){
        error = "As transcriptome source is reference guided you must include a ref_genome parameter"
    }
    if (params.ref_genome){
        ref_genome = file(params.ref_genome, type: "file")
        if (!ref_genome.exists()) {
            error = "--ref_genome: File doesn't exist, check path."
        }
    }else {
        ref_genome = file("$projectDir/data/OPTIONAL_FILE")
    }
    if (params.transcriptome_source == "denovo" && params.ref_annotation) {
        error = "Reference annotation with de denovo assembly is not supported"
    }

    if (params.ref_annotation){
        ref_annotation = file(params.ref_annotation, type: "file")

        if (!ref_annotation.exists()) {
            error = "--ref_annotation: File doesn't exist, check path."
        }
        use_ref_ann = true
    }else{
        ref_annotation= file("$projectDir/data/OPTIONAL_FILE")
        use_ref_ann = false
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
        log.info("Reference Transcriptome provided will be used for differential expression.")
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
        throw new Exception(error)
    }else{
        reads = samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "analyse_unclassified":params.analyse_unclassified,
        "fastcat_stats": true,
        "fastcat_extra_args": ""])

        pipeline(reads, ref_genome, ref_annotation,
            jaffal_refBase, params.jaffal_genome, params.jaffal_annotation,
            condition_sheet, ref_transcriptome, use_ref_ann)

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
