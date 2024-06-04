#!/usr/bin/env nextflow

/* This workflow is a adapted from two previous pipeline written in Snakemake:
- https://github.com/nanoporetech/pipeline-nanopore-ref-isoforms
*/

import groovy.json.JsonBuilder;
import nextflow.util.BlankSeparatedList;
import java.util.ArrayList;
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { reference_assembly } from './subworkflows/reference_assembly'
include { gene_fusions } from './subworkflows/JAFFAL/gene_fusions'
include { differential_expression } from './subworkflows/differential_expression'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process getVersions {
    label "isoforms"
    cpus 1
    memory "2 GB"
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
    """
}


process getParams {
    label "isoforms"
    cpus 1
    memory "2 GB"
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
    memory "2 GB"
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
    memory "2 GB"
    input:
        path compressed_annotation
    output:
        path "${compressed_annotation.baseName}"
    """
    gzip -df ${compressed_annotation}
    """
}



process decompress_transcriptome {
    label "isoforms"
    cpus 1
    memory "2 GB"
    input:
        path "compressed_ref.gz"
    output:
        path "compressed_ref", emit: decompressed_ref
    """
    gzip -df "compressed_ref.gz"
    """
}


// Remove empty transcript ID fields
process preprocess_ref_annotation {
    label "isoforms"
    cpus 1
    memory "2 GB"
    input:
        path ref_annotation
    output:
        path "amended.${ref_annotation}"
    """
    sed -i -e 's/transcript_id "";//g' ${ref_annotation}
    mv ${ref_annotation} "amended.${ref_annotation}"
    """
}

// Just keep transcript ID for each transcriptome fasta
process preprocess_ref_transcriptome {
    label "isoforms"
    cpus 1
    memory "2 GB"
    input:
        path "ref_transcriptome"
    output:
        path "amended.${ref_transcriptome}"
    """
    sed -i -e 's/|.*//' ${ref_transcriptome}
    mv ${ref_transcriptome} "amended.${ref_transcriptome}"
    """
}



process preprocess_reads {
    /*
    Concatenate reads from a sample directory.
    Optionally classify, trim, and orient cDNA reads using pychopper
    */

    label "isoforms"
    cpus 4
    memory "2 GB"
    input:
        tuple val(meta), path('seqs.fastq.gz')
    output:
        tuple val("${meta.alias}"),
              path("${meta.alias}_pychopper_output/${meta.alias}_full_length_reads.fastq"),
              emit: full_len_reads
        tuple val("${meta.alias}"),
              path("${meta.alias}_pychopper_output/"),
              emit: pychopper_output
        path("${meta.alias}_pychopper_output/pychopper.tsv"),
              emit: report
    script:
        def cdna_kit = params.cdna_kit.split("-")[-1]
	    def extra_params = params.pychopper_opts ?: ''
        """
        pychopper -t ${params.threads} -k ${cdna_kit} -m ${params.pychopper_backend} ${extra_params} 'seqs.fastq.gz' ${meta.alias}_full_length_reads.fastq
        workflow-glue generate_pychopper_stats --data pychopper.tsv --output .

        # Add sample id column
        sed "1s/\$/\tsample_id/; 1 ! s/\$/\t${meta.alias}/" pychopper.tsv > tmp
        mv tmp pychopper.tsv


        mkdir "${meta.alias}_pychopper_output/"
        shopt -s extglob  # Allow extended pattern matching so we can exclude files from the mv
        mv !("${meta.alias}_pychopper_output"|seqs.fastq.gz) "${meta.alias}_pychopper_output/"
        """
}

process build_minimap_index{
    /*
    Build minimap index from reference genome
    */
    label "isoforms"
    cpus params.threads
    memory "31 GB"

    input:
        path reference
    output:
        path "genome_index.mmi", emit: index
    script:
    """
    minimap2 -t ${params.threads} ${params.minimap2_index_opts} -I 1000G -d "genome_index.mmi" ${reference}
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
    memory "2 GB"

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
    memory "2 GB"

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
    cpus params.threads
    memory "2 GB"

    input:
        tuple val(sample_id), path (gff_bundle)
    output:
        tuple val(sample_id), path("${sample_id}.gff"), emit: gff
    script:
    def merged_gff = "${sample_id}.gff"
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
    cpus 1
    memory "2 GB"
    input:
       tuple val(sample_id), path(query_annotation)
       path ref_annotation
    output:
        tuple val(sample_id), path("${sample_id}"), emit: gffcmp_dir
        path ("${sample_id}_annotated.gtf"), emit: gtf
        tuple val(sample_id), path("${sample_id}_transcripts_table.tsv"),
            emit: isoforms_table
    script:
    def out_dir = "${sample_id}"
    """
        mkdir $out_dir
        echo "Doing comparison of reference annotation: ${ref_annotation} and the query annotation"

        gffcompare -o ${out_dir}/str_merged -r ${ref_annotation} \
            ${params.gffcompare_opts} ${query_annotation}

        workflow-glue generate_tracking_summary --tracking $out_dir/str_merged.tracking \
            --output_dir ${out_dir} --annotation ${ref_annotation}

        mv *.tmap "${out_dir}"
        mv *.refmap "${out_dir}"
        cp "${out_dir}/str_merged.annotated.gtf" "${sample_id}_annotated.gtf"

        # Make an isoform table for report and user output.
        workflow-glue make_isoform_table \
            --sample_id "${sample_id}" \
            --gffcompare_dir "${out_dir}"
        """
}

process get_transcriptome{
        /*
        Write out a transcriptome file based on the query gff annotations.
        */
        label 'isoforms'
        cpus 1
        memory "2 GB"
        input:
            tuple val(sample_id), path(transcripts_gff), path(gffcmp_dir), path(reference_seq)
        output:
            tuple val(sample_id), path("*.fas"), emit: transcriptome

        script:
        def transcriptome = "${sample_id}_transcriptome.fas"
        def merged_transcriptome = "${sample_id}_merged_transcriptome.fas"
        // if no ref_annotation gffcmp_dir will be optional file
        // so skip getting transcriptome FASTA from the annotated files.
        if (params.ref_annotation){
        """
        gffread -g ${reference_seq} -w ${transcriptome} ${transcripts_gff}
        if  [ "\$(ls -A $gffcmp_dir)" ];
            then
                gffread -F -g ${reference_seq} -w ${merged_transcriptome} $gffcmp_dir/str_merged.annotated.gtf
        fi
        """
        } else {
        """
        gffread -g ${reference_seq} -w ${transcriptome} ${transcripts_gff}
        """
        }
}

process merge_transcriptomes {
    // Merge the transcriptomes from all samples
    label 'isoforms'
    cpus 2
    memory "2 GB"
    input:
        path "query_annotations/*"
        path ref_annotation
        path ref_genome
    output:
        path "final_non_redundant_transcriptome.fasta", emit: fasta
        path "stringtie.gtf", emit: gtf
    """
    stringtie --merge -G $ref_annotation -p ${task.cpus} -o stringtie.gtf query_annotations/*
    seqkit subseq --feature "transcript" --gtf-tag "transcript_id" --gtf stringtie.gtf $ref_genome > temp_transcriptome.fasta
    seqkit rmdup -s < temp_transcriptome.fasta > temp_del_repeats.fasta
    cat temp_del_repeats.fasta | sed 's/>.* />/'  | sed -e 's/_[0-9]* \\[/ \\[/' > temp_rm_empty_seq.fasta
    awk 'BEGIN {RS = ">" ; FS = "\\n" ; ORS = ""} \$2 {print ">"\$0}' temp_rm_empty_seq.fasta > "final_non_redundant_transcriptome.fasta"
    rm temp_transcriptome.fasta
    rm temp_del_repeats.fasta
    rm temp_rm_empty_seq.fasta
    """
}

process makeReport {

    label "isoforms"
    cpus 2
    memory "4 GB"

    input:
        path versions
        path "params.json"
        path "pychopper_report/*"
        path"jaffal_csv/*"
        path "per_read_stats/?.gz"
        path "aln_stats/*"
        path "gffcmp_dir/*"
        path "gff_annotation/*"
        path "de_report/*"
        path "seqkit/*"
        path "isoforms_table/*"
    output:
        path ("wf-transcriptomes-*.html"), emit: report
        // If de analysis has been run output the counts files with gene name added.
        path ("results_dge.tsv"), emit: results_dge, optional: true
        path ("unfiltered_tpm_transcript_counts.tsv"), emit: tpm, optional: true
        path ("unfiltered_transcript_counts_with_genes.tsv"), emit: unfiltered, optional: true
        path ("filtered_transcript_counts_with_genes.tsv"), emit: filtered, optional: true
        path ("all_gene_counts.tsv"), emit: gene_counts, optional: true
    shell:
        report_name = "wf-transcriptomes-report.html"
    '''
    if [ -f "de_report/OPTIONAL_FILE" ]; then
        dereport=""
    else
        dereport="--de_report true --de_stats "seqkit/*""
        mv de_report/*.g*f* de_report/stringtie_merged.gtf
    fi
    if [ -f "gff_annotation/OPTIONAL_FILE" ]; then
        OPT_GFF_ANNOTATION=""
    else
        OPT_GFF_ANNOTATION="--gff_annotation gff_annotation/*"
    fi
    if [ -f "gffcmp_dir/OPTIONAL_FILE" ]; then
        OPT_GFFCMP_DIR=""
    else
        OPT_GFFCMP_DIR="--gffcompare_dir gffcmp_dir/"
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
    if [ -f "isoforms_table/OPTIONAL_FILE" ]; then
        OPT_ISO_TABLE=""
    else
        OPT_ISO_TABLE="--isoform_table isoforms_table"
    fi
    workflow-glue report --report !{report_name} \
    --versions !{versions} \
    --params params.json \
    ${OPT_ALN} \
    ${OPT_PC_REPORT} \
    --stats per_read_stats/* \
    ${OPT_GFF_ANNOTATION} \
    ${OPT_ISO_TABLE} \
    ${OPT_GFFCMP_DIR} \
    --isoform_table_nrows !{params.isoform_table_nrows} \
    ${OPT_JAFFAL_CSV} \
    ${dereport}
    '''
}


// Creates a new directory named after the sample alias and moves the fastcat results
// into it.
process collectFastqIngressResultsInDir {
    label "isoforms"
    cpus 1
    memory "2 GB"
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
    cpus 1
    memory "2 GB"
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

     
        results = Channel.empty()
        software_versions = getVersions()
        workflow_params = getParams()
        input_reads = reads.map{ meta, samples, stats -> [meta, samples]}
        sample_ids = input_reads.flatMap({meta,samples -> meta.alias})
        per_read_stats = reads.map{ meta, samples, stats -> stats.resolve("per-read-stats.tsv.gz") }.toList()

        if (!params.direct_rna){
            preprocess_reads(input_reads)
            full_len_reads = preprocess_reads.out.full_len_reads
            pychopper_report = preprocess_reads.out.report.collectFile(keepHeader: true)
            pychopper_results_dir = preprocess_reads.out.pychopper_output.map{ it -> it[1]}
            results = results.concat(pychopper_results_dir)
        }
        else{
            full_len_reads = input_reads.map{ meta, reads -> [meta.alias, reads]}
            pychopper_report = OPTIONAL_FILE
        }
        if (params.transcriptome_source != "precomputed"){
            build_minimap_index(ref_genome)
            log.info("Doing reference based transcript analysis")
            assembly = reference_assembly(build_minimap_index.out.index, ref_genome, full_len_reads)
        
            assembly_stats = assembly.stats.map{ it -> it[1]}.collect()
     
            split_bam(assembly.bam.map {sample_id, bam, bai -> [sample_id, bam]})

            assemble_transcripts(split_bam.out.bundles.flatMap(map_sample_ids_cls).combine(ref_annotation),use_ref_ann)

            merge_gff_bundles(assemble_transcripts.out.gff_bundles.groupTuple())
            // only run gffcompare if ref annotation provided. Otherwise create optional files and channels
            if (params.ref_annotation){
                run_gffcompare(merge_gff_bundles.out.gff, ref_annotation)
                gff_compare_dir = run_gffcompare.out.gffcmp_dir
                gff_compare = run_gffcompare.out.gffcmp_dir.map{ it -> it[1]}.collect()
                isoforms_table = run_gffcompare.out.isoforms_table.map{ it -> it[1]}.collect()
                // create per sample gff tuples with gff compare directories
                gff_tuple = merge_gff_bundles.out.gff
                .join(gff_compare_dir)
            } else {
                // create per sample gff tuples with optional files as no ref_annotation
                optional_channel = Channel.fromPath("$projectDir/data/OPTIONAL_FILE")
                gff_tuple = merge_gff_bundles.out.gff.combine(optional_channel)
                gff_compare = OPTIONAL_FILE
                isoforms_table = OPTIONAL_FILE
            }
            // For reference based assembly, there is only one reference
            // So map this reference to all sample_ids
            seq_for_transcriptome_build = sample_ids.flatten().combine(ref_genome)
            get_transcriptome(
                gff_tuple
                .join(seq_for_transcriptome_build))


            merge_gff = merge_gff_bundles.out.gff.map{ it -> it[1]}.collect()
            results = results.concat(assembly.bam.map {sample_id, bam, bai -> [bam, bai]}.flatten())
        }
        else{
            gff_compare = OPTIONAL_FILE
            isoforms_table = OPTIONAL_FILE
            merge_gff = OPTIONAL_FILE
            assembly_stats = OPTIONAL_FILE
            use_ref_ann = false
 
        }
        if (jaffal_refBase){
                gene_fusions(full_len_reads, jaffal_refBase, jaffal_genome, jaffal_annotation)
                jaffal_out = gene_fusions.out.results_csv.map{ alias, csv -> csv}.collectFile(keepHeader: true, name: 'jaffal.csv')
            } else{
                jaffal_out = OPTIONAL_FILE
        }

        if (params.de_analysis){
            sample_sheet = file(params.sample_sheet, type:"file")
            if (!params.ref_transcriptome){
                merge_transcriptomes(run_gffcompare.output.gtf.collect(), ref_annotation, ref_genome)
                transcriptome = merge_transcriptomes.out.fasta
                gtf = merge_transcriptomes.out.gtf
            }
            else {
                transcriptome =  Channel.fromPath(ref_transcriptome)
                if (file(params.ref_transcriptome).extension == "gz") {
                    transcriptome = decompress_transcriptome(ref_transcriptome)
                }
                transcriptome = preprocess_ref_transcriptome(transcriptome)
                gtf = ref_annotation
            }
            de = differential_expression(transcriptome, input_reads, sample_sheet, gtf)
            de_report = de.all_de
            count_transcripts_file = de.count_transcripts
            dtu_plots = de.dtu_plots
            de_outputs = de.de_outputs
            counts = de.counts
        } else{
            de_report = OPTIONAL_FILE
            count_transcripts_file = OPTIONAL_FILE
        }
        
        makeReport(
            software_versions,
            workflow_params,
            pychopper_report,
            jaffal_out,
            per_read_stats,
            assembly_stats,
            gff_compare,
            merge_gff,
            de_report,
            count_transcripts_file,
            isoforms_table)

       report = makeReport.out.report

       results = results.concat(report)
       
       if (use_ref_ann){
            results = run_gffcompare.output.gffcmp_dir.concat(
                      assembly.stats,
                      run_gffcompare.out.isoforms_table,
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
        if (params.jaffal_refBase){
            results = results
                .concat(gene_fusions.out.results
                .map {it -> it[1]})
        }
        
       results = results.map{ [it, null] }.concat(fastq_ingress_results.map { [it, "fastq_ingress_results"] })
        
        if (params.de_analysis){
           de_results = report.concat(
            transcriptome, de_outputs.flatten(), counts.flatten(),
            makeReport.out.results_dge,  makeReport.out.tpm,
            makeReport.out.filtered,  makeReport.out.unfiltered,
            makeReport.out.gene_counts)
            // Output de_analysis results in the dedicated directory.
            results = results.concat(de_results.map{ [it, "de_analysis"] })
        }

        results.concat(workflow_params.map{ [it, null]})

       
    emit:
        results
}

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

    error = null

    if (params.containsKey("minimap_index_opts")) {
        error = "`--minimap_index_opts` parameter is deprecated. Use parameter `--minimap2_index_opts` instead."
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
        ref_genome = OPTIONAL_FILE
    }
    if (params.containsValue("denovo")) {
        error = "Denovo transcriptome source is no longer supported. Please use the reference-guided or precomputed options."
    }

    if (params.ref_annotation){
        ref_annotation = file(params.ref_annotation, type: "file")

        if (!ref_annotation.exists()) {
            error = "--ref_annotation: File doesn't exist, check path."
        }
        use_ref_ann = true
    }else{
        ref_annotation= OPTIONAL_FILE
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
    ref_transcriptome = OPTIONAL_FILE
    if (params.ref_transcriptome){
        log.info("Reference Transcriptome provided will be used for differential expression.")
        ref_transcriptome = file(params.ref_transcriptome, type:"file")
    }
    if (params.de_analysis){
        if (!params.ref_annotation){
            error = "When running in --de_analysis mode you must provide a reference annotation."
        }
        if (!params.sample_sheet){
            error = "You must provide a sample_sheet with at least alias and condition columns."
        }
        if (params.containsKey("condition_sheet")) {
        error = "Condition sheets have been deprecated. Please add a 'condition' column to your sample sheet instead. Check the quickstart for more information."
        }
    } else{
        if (!params.ref_annotation){
            log.info("Warning: As no --ref_annotation was provided, the output transcripts will not be annotated.")
        }
    }
    if (error){
        throw new Exception(error)
    }

    if (params.fastq) {
        samples = fastq_ingress([
            "input":params.fastq,
            "sample":params.sample,
            "sample_sheet":params.sample_sheet,
            "analyse_unclassified":params.analyse_unclassified,
            "stats": true,
            "fastcat_extra_args": "",
            "per_read_stats": true])
    } else {
        samples = xam_ingress([
            "input":params.bam,
            "sample":params.sample,
            "sample_sheet":params.sample_sheet,
            "analyse_unclassified":params.analyse_unclassified,
            "keep_unaligned": true,
            "return_fastq": true,
            "stats": true,
            "per_read_stats": true])
    }

    pipeline(samples, ref_genome, ref_annotation,
        jaffal_refBase, params.jaffal_genome, params.jaffal_annotation,
        ref_transcriptome, use_ref_ann)

    output(pipeline.out.results)

}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
