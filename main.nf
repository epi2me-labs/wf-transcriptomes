#!/usr/bin/env nextflow

/* This workflow is a adapted from two previous pipeline written in Snakemake:
- https://github.com/nanoporetech/pipeline-nanopore-ref-isoforms
*/

import groovy.json.JsonBuilder;
import nextflow.util.BlankSeparatedList;
import java.util.ArrayList;
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { configure_igv } from './lib/common'
include { reference_assembly } from './subworkflows/reference_assembly'
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
    cpus params.threads
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
    memory "4 GB"

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
    -o  ${prefix}.gff -l ${prefix} ${bam}
     """
}


process merge_gff_bundles{
    /*
    Merge gff bundles into a single gff file per sample, and get summary statistics
    */
    label 'isoforms'
    cpus params.threads
    memory "2 GB"

    input:
        tuple val(sample_id), path (gff_bundle)
    output:
        tuple val(sample_id), path("${sample_id}.gff"), emit: gff
        tuple val(sample_id), path("transcriptome_summary.pickle"), emit: summary
    script:
    def merged_gff = "${sample_id}.gff"
    """
    echo '##gff-version 2' >> $merged_gff;
    echo '#pipeline-nanopore-isoforms: stringtie' >> $merged_gff;

    for fn in ${gff_bundle};
    do
        grep -v '#' \$fn >> $merged_gff

    done

    workflow-glue summarise_gff \
        $merged_gff \
        $sample_id \
        transcriptome_summary.pickle
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

        mv *.tmap "${out_dir}"
        mv *.refmap "${out_dir}"
        cp "${out_dir}/str_merged.annotated.gtf" "${sample_id}_annotated.gtf"

        workflow-glue parse_gffcompare \
            --sample_id "${sample_id}" \
            --gffcompare_dir "${out_dir}" \
            --isoform_table_out "${sample_id}_transcripts_table.tsv" \
            --tracking $out_dir/str_merged.tracking \
            --annotation ${ref_annotation}
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
            tuple val(sample_id), path("transcripts.gff"), path(gffcompare_dir), path("reference.fa")
        output:
            tuple val(sample_id), path("*transcriptome.fas"), emit: transcriptome
        script:
        def transcriptome = "${sample_id}_transcriptome.fas"
        def merged_transcriptome = "${sample_id}_merged_transcriptome.fas"
        // if no ref_annotation gffcmp_dir will be optional file
        // so skip getting transcriptome FASTA from the annotated files.
        if (params.ref_annotation){
        """
        gffread -F -g reference.fa -w ${merged_transcriptome} $gffcompare_dir/str_merged.annotated.gtf
        """
        } else {
        """
        gffread -g reference.fa -w ${transcriptome} "transcripts.gff"
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
    stringtie --rf --merge -G $ref_annotation -p ${task.cpus} -o stringtie.gtf query_annotations/*
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

    label "wf_common"
    cpus 2
    memory "4 GB"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "wf-transcriptomes-report.html"
    input:
        val metadata
        path stats, stageAs: "stats_*"
        path versions
        val wf_version
        path "params.json"
        path "transcriptome_aln_stats/*"
        path pychopper, stageAs: "pychopper_report/*"
        path aln_stats, stageAs: "aln_stats/*"
        path gffcmp_dir, stageAs: "gffcmp_dir/*"
        path gff_annotation, stageAs: "gff_annotation/*"
        path de_report, stageAs: "de_report/*"
        path isoforms_table, stageAs: "isoforms_table/*"
        path transcriptome_summary, stageAs: "transcriptome_summary/summary_*.pkl"

    output:
        path ("wf-transcriptomes-*.html"), emit: report
        path ("results_dge.tsv"), emit: results_dge, optional: true
        path ("unfiltered_tpm_transcript_counts.tsv"), emit: tpm, optional: true
        path ("unfiltered_transcript_counts_with_genes.tsv"), emit: unfiltered, optional: true
        path ("filtered_transcript_counts_with_genes.tsv"), emit: filtered, optional: true
        path ("all_gene_counts.tsv"), emit: gene_counts, optional: true
    script:
        String report_name = "wf-transcriptomes-report.html"
        String metadata = new JsonBuilder(metadata).toPrettyString()
        String gff_opts = gff_annotation.fileName.name == OPTIONAL_FILE.name ? "" :  "--gff_annotation gff_annotation/"
        String de_report_opts = de_report.fileName.name == OPTIONAL_FILE.name ? "" : "--de_report de_report/ --de_stats transcriptome_aln_stats/"
        String gffcmp_opts = gffcmp_dir.fileName.name == OPTIONAL_FILE.name  ? "" :  "--gffcompare_dir gffcmp_dir/"
        String aln_stats_opts = aln_stats.fileName.name == OPTIONAL_FILE.name ? "" :  "--alignment_stats aln_stats/"
        String pychop_opts = pychopper.fileName.name == OPTIONAL_FILE.name ? "" :  "--pychop_report pychopper_report/"
        String iso_table_opts = isoforms_table.fileName.name == OPTIONAL_FILE.name ? "" :  "--isoform_table isoforms_table/"
        String tr_summary_opts = transcriptome_summary.fileName.name == OPTIONAL_FILE.name ? "" :  "--transcriptome_summary transcriptome_summary/"
    """
    echo '${metadata}' > metadata.json
    workflow-glue report \
        --report $report_name \
        --versions $versions \
        --wf_version $wf_version \
        --params params.json \
        $aln_stats_opts \
        $pychop_opts \
        --stats $stats \
        --metadata metadata.json \
        $gff_opts \
        $iso_table_opts \
        $gffcmp_opts \
        --isoform_table_nrows ${params.isoform_table_nrows} \
        $de_report_opts \
        $tr_summary_opts
    """
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
process publish_results {
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


// Check ref_annotation transcript strand column for "." if in de_analysis mode
process check_annotation_strand {
    label "isoforms"
    cpus 1
    memory "2 GB"
    input:
        path "ref_annotation.gtf"
    output:
        tuple stdout, path("ref_annotation.gtf")   
    """
    awk '{if (\$3=="transcript" && \$7 != "+" && \$7 != "-") print \$3, \$7}' "ref_annotation.gtf"
    """
}


// Process to create the faidx index
process faidx {
    // If the input file is gzipped, we need to emit the indexes for the input gzip file
    // only. Therefore, this become redundant to be emitted as it won't be used by the
    // IGV configuration, but only by internal processes together with the decompressed
    // FASTA file. To avoid unnecessary emissions, we enable only if the input file is
    // decompressed.
    publishDir "${params.out_dir}/igv_reference", mode: 'copy', pattern: "*", enabled: !params.ref_genome.toLowerCase().endsWith("gz")
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        path(ref)
    output:
        path("${ref}.fai")
    script:
    """
    samtools faidx ${ref}
    """
}


// Process to create the faidx indexes for a gzipped reference
process gz_faidx {
    publishDir "${params.out_dir}/igv_reference", mode: 'copy', pattern: "*"
    label "wf_common"
    cpus 1
    memory 4.GB
    // If a user provides a non-bgzipped file, the process won't
    // generate the indexes. We should tolerate that, still avoid emitting
    // the reference and simply have a broken IGV file.
    // The gzi is not required to operate the workflow, so we actually tolerate any failure.
    errorStrategy 'ignore'
    input:
        path(ref)
    output:
        tuple path("${ref}.fai"), path("${ref}.gzi")
    script:
    """
    samtools faidx ${ref}
    """
}



// workflow module
workflow pipeline {
    take:
        reads
        ref_genome
        ref_annotation
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
        | collectFastqIngressResultsInDir
        // fastq_ingress doesn't have the index; add one extra null for compatibility.
        // We do not use variable name as assigning variable name with a tuple
        // not matching (e.g. meta, bam, bai, stats <- [meta, bam, stats]) causes
        // the workflow to crash.
        reads = reads
        .map{
            it.size() == 4 ? it : [it[0], it[1], null, it[2]]
        }
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
        // Define BAM output Directory
        String publish_bams = "BAMS"
        software_versions = getVersions()
        workflow_params = getParams()
        input_reads = reads.map{ meta, samples, index, stats -> [meta, samples]}
        sample_ids = input_reads.flatMap({meta,samples -> meta.alias})

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
            transcriptome_summary = merge_gff_bundles.out.summary.map {it[1]}.collect()

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
            // Output BAMS in a dedicated directory
            bam_results = assembly.bam.map{
                sample_id, bam, bai -> [bam, bai]}.flatten().map{ [it, publish_bams] }
        }
        else{
            gff_compare = OPTIONAL_FILE
            isoforms_table = OPTIONAL_FILE
            merge_gff = OPTIONAL_FILE
            assembly_stats = OPTIONAL_FILE
            transcriptome_summary = OPTIONAL_FILE
            use_ref_ann = false
        }
        if (params.de_analysis){
            sample_sheet = file(params.sample_sheet, type:"file")
            // check ref annotation contains only + or - strand as DE analysis will error on .
            check_annotation_strand(ref_annotation).map { stdoutput, annotation ->
            // check if there was an error message
            if (stdoutput) error "In ref_annotation, transcript features must have a strand of either '+' or '-'."
                    stdoutput
                }
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
            de_outputs = de.de_outputs
            count_transcripts_file = de.count_transcripts
        } else{
            de_report = OPTIONAL_FILE
            count_transcripts_file = OPTIONAL_FILE
        }

        // get metadata and stats files, keeping them ordered (could do with transpose I suppose)
        reads.multiMap{ meta, path, index, stats ->
            meta: meta
            stats: stats
        }.set { for_report }
        metadata = for_report.meta.collect()
        stats = for_report.stats.collect()

        makeReport(
            metadata,
            stats,
            software_versions,
            workflow.manifest.version,
            workflow_params,
            count_transcripts_file,
            pychopper_report,
            assembly_stats,
            gff_compare,
            merge_gff,
            de_report,
            isoforms_table,
            transcriptome_summary)

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

       results = results.map{ [it, null] }.concat(fastq_ingress_results.map { [it, "fastq_ingress_results"] })

        if (params.de_analysis){
           de_results = report.concat(
            transcriptome, de_outputs.flatten(),
            makeReport.out.results_dge,  makeReport.out.tpm,
            makeReport.out.filtered,  makeReport.out.unfiltered,
            makeReport.out.gene_counts)
            // Output de_analysis results in the dedicated directory.
            results = results.concat(de_results.map{ [it, "de_analysis"] })
        }


        results.concat(workflow_params.map{ [it, null]})
        // IGV config
        if (params.transcriptome_source == "precomputed" && params.igv){
            log.warn("IGV configuration does not work if transcriptome sources is set to `precomputed`.")
        }
        if (params.transcriptome_source != "precomputed" && params.igv){
            is_compressed = file("${params.ref_genome}").extension == "gz"
            String publish_ref = "igv_reference"
            reference_genome = Channel.fromPath("${params.ref_genome}")
            igv_ref = reference_genome | flatten | map { it -> "${it.toUriString()}" }
            if (is_compressed){
                // Define indexes names.
                String input_fai_index = "${params.ref_genome}.fai"
                String input_gzi_index = "${params.ref_genome}.gzi"
                // Check whether the input gzref is indexed. If so, pass these as indexes.
                // Otherwise, generate the gzip + fai indexes for the compressed reference.
                if (file(input_fai_index).exists() && file(input_gzi_index).exists()){
                    gzindexes = Channel.fromPath(input_fai_index)
                    | mix(
                        Channel.fromPath(input_gzi_index)
                    )
                    gz_igv = gzindexes | flatten | map { it -> "${it.toUriString()}" }
                } else {
                    gz_igv = gz_faidx(Channel.fromPath("${params.ref_genome}"))
                    | flatten
                    | map  { it -> "$publish_ref/${it.Name}" }
                    gz_igv | ifEmpty{
                        if (params.containsKey("igv") && params.igv){
                            log.warn """\
                                The input reference is compressed but not with bgzip, which is required to create an index.
                                The workflow will proceed but it will not be possible to load the reference in the IGV Viewer.
                                To use the IGV Viewer, provide an uncompressed, or bgzip compressed version of the input reference next time you run the workflow.
                                """.stripIndent()
                        }
                    }
                }
            } else {
                gzindexes = Channel.empty()
                gz_igv = Channel.empty()
            }

            // Generate fai index if the file is either compressed, or if fai doesn't exists
            if (!is_compressed && file("${params.ref_genome}.fai").exists()){
                ref_idx = Channel.fromPath("${params.ref_genome}.fai")
                igv_index = ref_idx | flatten | map { it -> "${it.toUriString()}" }
            } else {
                ref_idx = faidx(reference_genome)
                igv_index = ref_idx | map { it -> "$publish_ref/${it.Name}" }
            }

            // get list of file names
            // Absolute paths required for directories
            igv_files = reads
            | map { meta, sample, index, stats -> meta.alias }
            | toSortedList
            | map { list -> list.collect{
                [
                 "$publish_bams/${it}_reads_aln_sorted.bam",
                 "$publish_bams/${it}_reads_aln_sorted.bam.bai"
                ]
            } }
            | concat (igv_ref)
            | flatten
            | concat ( igv_index)
            | concat (gz_igv)
            | flatten
            | collectFile(name: "file-names.txt", newLine: true, sort: false)

            // configure IGV
            igv_conf = configure_igv(
                igv_files,
                Channel.of(null), // igv locus
                [displayMode: "SQUISHED", colorBy: "strand"], // bam extra opts
                Channel.of(null), // vcf extra opts
                )

            results = results.concat(igv_conf.map{ [it, null]})
            results = results.concat(bam_results)
      }

       
    emit:
        results
}

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

    error = null

    if (params.containsValue("jaffal_refBase")) {
        error = "JAFFAL fusion detection has been removed from this workflow."
    }
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

    pipeline(samples, ref_genome, ref_annotation, ref_transcriptome, use_ref_ann)
    publish_results(pipeline.out.results)
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
