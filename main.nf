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

process validate_ref_annotation {
    label "isoforms"
    cpus 1
    memory "2 GB"
    input:
        path "annotation.gtf"
        path "reference.fasta"
    output: 
        stdout
    // Checks for overlap between seq_id column in annotation gtf and fasta reference ID's
    // If no overlap is found exit
    // Partial overlap (eg. user supplies genes/contigs of interest in annotation but the genome sequence) - warning
    script:
    """
    grep -v '^#' annotation.gtf | cut -f1  | sort -u > seq_ids.txt
    awk '/^>/ {print substr(\$1,2)}' reference.fasta | sort -u > ref_ids.txt
    matches=\$(comm -12 seq_ids.txt ref_ids.txt)
    only_in_annotation=\$(comm -23 seq_ids.txt ref_ids.txt)
    only_in_reference=\$(comm -13 seq_ids.txt ref_ids.txt)
    if [[ -z "\$matches" ]]; then
        echo "
        ERROR: Seqid mismatch found between the provided ref_annotation (GTF/GFF)
        file and ref_genome (FASTA).
        For the reference guided differential expression subworkflow they must overlap.
        " >&2
        echo "Annotation ID examples:"
        head -n 5 seq_ids.txt
        echo "Reference ID examples:"
        head -n 5 ref_ids.txt
        echo "We recommend getting both files from the same source.
        eg. both from Ensembl or both from NCBI.
        Alternatively provide a pre-computed transcriptome using the ref_transcriptome parameter
        See the README for more details on which inputs are supported."
        exit 78 
    fi
    if [[ -n "\$only_in_annotation" ]]; then
        echo "Warning: Some sequence IDs are only present in the reference annotation and not the
             reference genome so will not be used in downstream analysis eg."
        echo "\$only_in_annotation" | head -n 5
    fi
    if [[ -n "\$only_in_reference" ]]; then
        echo "Warning: Some FASTA reference IDs are only present in the reference genome
            and not the reference annotation so will not be used in downstream analysis eg."
        echo "\$only_in_reference" | head -n 5
    fi
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

    */

    label 'isoforms'
    cpus params.threads
    memory "15 GB"

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
    cpus Math.min(params.threads, 6)
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
        tuple val(sample_id), path ('gff_bundles/annotation*.gff')
    output:
        tuple val(sample_id), path("${sample_id}.gff"), emit: gff
        tuple val(sample_id), path("transcriptome_summary.pickle"), emit: summary
    script:
    def merged_gff = "${sample_id}.gff"
    """
    echo '##gff-version 2' >> $merged_gff;
    echo '#pipeline-nanopore-isoforms: stringtie' >> $merged_gff;

    find -L gff_bundles -type f -name "*.gff" \
        -exec awk '!/^#/ {print}' {} \\; >> "${sample_id}.gff"
    if ! [ -s "${sample_id}.gff" ]; then
        echo "No transcripts found for ${sample_id}"
        # This is unlikely to ever happen, but if it does, we should fail the workflow.
        exit 70
    fi

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
    cpus Math.min(params.threads, 6)
    memory "2 GB"
    input:
        path "query_annotations/*"
        path ref_annotation
        path ref_genome
    output:
        path "final_non_redundant_transcriptome.fasta", emit: fasta
        path "stringtie.gtf", emit: gtf
    """
    stringtie --merge -G "${ref_annotation}" -p ${task.cpus} -o stringtie.gtf query_annotations/*
    gffread -g "${ref_genome}" -w "final_non_redundant_transcriptome.fasta" "stringtie.gtf"
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
        raw_reads_ch // from initial fastq_ingress or xam_ingress in the entrypoint workflow
        ref_genome_file
        ref_annotation_file
        ref_transcriptome_file
        use_ref_ann_bool
    main:
        // Define input mode
        def is_bam_input = params.bam_input_dir && params.sample_sheet
        def is_fastq_input = params.fastq && params.sample_sheet

        if (is_bam_input && is_fastq_input) {
            error "Both --bam_input_dir and --fastq are specified with --sample_sheet. Please provide only one input type."
        } else if (is_bam_input && params.bam) {
            error "--bam_input_dir and --bam cannot be used together. Use --sample_sheet with --bam_input_dir for multiple BAMs."
        } else if (is_fastq_input && params.bam) {
             error "--fastq and --bam cannot be used together. For single BAM, just use --bam. For multiple BAMs via sheet, use --bam_input_dir."
        } else if (!is_bam_input && !is_fastq_input && !params.bam) { //
            error "No valid input provided. Please specify --fastq and --sample_sheet, or --bam_input_dir and --sample_sheet, or a single --bam."
        }


        // Initialize conditional output channels
        Channel PYCHOPPER_REPORT_CH = Channel.value(OPTIONAL_FILE)
        Channel ASSEMBLY_STATS_CH = Channel.value(OPTIONAL_FILE)
        Channel BAMS_FOR_SPLITTING = Channel.empty()
        Channel FULL_LEN_READS_FOR_DE = Channel.empty() // For DE analysis input, might be FASTQ or BAM

        // Process ref_genome and ref_annotation first
        Channel ref_genome_ch
        if (ref_genome_file && ref_genome_file.extension == "gz") {
            ref_genome_ch = decompress_ref(ref_genome_file)
        } else {
            ref_genome_ch = Channel.fromPath(ref_genome_file)
        }

        Channel ref_annotation_ch
        if (ref_annotation_file && ref_annotation_file.extension == "gz") {
            decompress_annot_ch = decompress_annotation(ref_annotation_file)
            ref_annotation_ch = preprocess_ref_annotation(decompress_annot_ch)
        } else {
            ref_annotation_ch = preprocess_ref_annotation(ref_annotation_file)
        }

        // Build minimap index if ref_genome is provided
        MINIMAP_INDEX_CH = Channel.empty()
        if (params.ref_genome) {
            MINIMAP_INDEX_CH = build_minimap_index(ref_genome_ch)
        } else if (params.transcriptome_source != "precomputed") {
            // This case should have been caught by entrypoint validation, but double check
            error "Reference genome (--ref_genome) is required when transcriptome_source is not 'precomputed'."
        }


        // This channel will hold [meta, file1, file2_or_null, stats_file_or_null] from ingress
        Channel INGRESS_RESULTS_CH = Channel.empty()

        if (is_bam_input) {
            log.info "Running in BAM input mode."
            // xam_ingress for BAMs via sample sheet
            // Assuming xam_ingress is updated or suitable for bam_input_dir
            // For now, we'll use the existing `reads` channel structure from the entrypoint for `xam_ingress`
            // if it's already configured for bam_input_dir there.
            // The key is that `raw_reads_ch` from the entry point needs to be the output of xam_ingress
            // when params.bam_input_dir is set.

            INGRESS_RESULTS_CH = raw_reads_ch // This comes from xam_ingress in WorkflowMain

            // Map BAMs for splitting and DE
            // Input: [meta, bam, bai, stats]
            BAMS_FOR_SPLITTING = INGRESS_RESULTS_CH
                .map { meta, bam, bai, stats ->
                    if (bam == null) return null // Should not happen if ingress is correct
                    [meta.alias, bam]
                }
                .filter { it != null }

            // Attempt to pass BAMs to DE. Subworkflow might need changes.
            FULL_LEN_READS_FOR_DE = BAMS_FOR_SPLITTING.map { id, bam -> [[alias:id], bam] }

            // Pychopper and assembly stats are not generated directly in this path
            PYCHOPPER_REPORT_CH = Channel.value(OPTIONAL_FILE)
            ASSEMBLY_STATS_CH = Channel.value(OPTIONAL_FILE)

        } else if (is_fastq_input) {
            log.info "Running in FASTQ input mode."
            // This is the original logic for FASTQ input
            INGRESS_RESULTS_CH = raw_reads_ch // This comes from fastq_ingress in WorkflowMain

            // Prepare input for preprocess_reads or direct use
            // input_reads_for_processing: [meta, fastq_path]
            input_reads_for_processing = INGRESS_RESULTS_CH.map { meta, fastq, index, stats -> [meta, fastq] }


            if (!params.direct_rna){
                preprocess_reads(input_reads_for_processing)
                FULL_LEN_READS_FOR_DE = preprocess_reads.out.full_len_reads // [sample_id, fastq_path]
                PYCHOPPER_REPORT_CH = preprocess_reads.out.report.collectFile(keepHeader: true)
                // Collect pychopper output directories for results
                pychopper_results_dir = preprocess_reads.out.pychopper_output.map{ it -> it[1]}
                results = results.concat(pychopper_results_dir)
            } else {
                // If direct RNA, use raw FASTQs for DE (after mapping to [sample_id, fastq_path])
                FULL_LEN_READS_FOR_DE = input_reads_for_processing.map{ meta, fastq -> [meta.alias, fastq]}
                PYCHOPPER_REPORT_CH = Channel.value(OPTIONAL_FILE)
            }

            if (params.transcriptome_source != "precomputed"){
                if (!MINIMAP_INDEX_CH.isEmpty()) { // Ensure index is available
                    assembly = reference_assembly(MINIMAP_INDEX_CH, ref_genome_ch, FULL_LEN_READS_FOR_DE, publish_prefix_bams)
                    ASSEMBLY_STATS_CH = assembly.stats.map{ it -> it[1]}.collect()
                    BAMS_FOR_SPLITTING = assembly.bam.map {sample_id, bam, bai -> [sample_id, bam]}
                } else {
                     // This should not happen if params.ref_genome was required
                    log.warn "Minimap index not available, skipping reference assembly for FASTQ input."
                    ASSEMBLY_STATS_CH = Channel.value(OPTIONAL_FILE) // ensure it's OPTIONAL_FILE
                    BAMS_FOR_SPLITTING = Channel.empty()
                }
            } else {
                // FASTQ input with precomputed transcriptome
                log.info "FASTQ input with precomputed transcriptome. BAMs for splitting will not be generated from reads."
                ASSEMBLY_STATS_CH = Channel.value(OPTIONAL_FILE) // ensure it's OPTIONAL_FILE
                BAMS_FOR_SPLITTING = Channel.empty()
            }

        } else if (params.bam) {
            // Single BAM file input (original simple BAM mode, not via bam_input_dir)
            // This mode might need further integration with the new channel structure if it's to be fully supported
            // alongside the sample sheet methods. For now, this subtask focuses on bam_input_dir.
            // We can assume `raw_reads_ch` from `WorkflowMain` handles this and provides a compatible structure.
            log.info "Running in single BAM input mode (params.bam)."
            INGRESS_RESULTS_CH = raw_reads_ch // This comes from xam_ingress (single BAM) in WorkflowMain

            BAMS_FOR_SPLITTING = INGRESS_RESULTS_CH
                .map { meta, bam, bai, stats -> // Assuming xam_ingress for single BAM also gives this structure
                    if (bam == null) return null
                    [meta.alias, bam]
                }
                .filter { it != null }
            FULL_LEN_READS_FOR_DE = BAMS_FOR_SPLITTING.map { id, bam -> [[alias:id], bam] }
            PYCHOPPER_REPORT_CH = Channel.value(OPTIONAL_FILE)
            ASSEMBLY_STATS_CH = Channel.value(OPTIONAL_FILE)

        } else {
            // Should have been caught by initial checks
            error "Pipeline started with an undetermined input configuration."
        }


        // Harmonize sample_ids and collect input files for report
        // INGRESS_RESULTS_CH is [meta, file1, file2_or_null, stats_file_or_null]
        // file1 is fastq for FASTQ mode, bam for BAM mode.
        // file2_or_null is bam_index for BAM mode.
        input_files_for_report = INGRESS_RESULTS_CH
            | map { meta, f1, f2, stats -> [meta, f1, stats] } // take meta, primary file, and stats
            | collectFastqIngressResultsInDir // Renamed for clarity, but re-used logic

        // Extract sample_ids for downstream processes
        // Ensure sample_ids are derived correctly regardless of input
        sample_ids = INGRESS_RESULTS_CH.flatMap({meta, f1, f2, stats -> meta.alias})


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
        String publish_prefix_bams = "BAMS" // Used by reference_assembly and IGV
        software_versions = getVersions()
        workflow_params = getParams()

        // `input_reads` and `full_len_reads` are effectively replaced by `FULL_LEN_READS_FOR_DE`
        // and `BAMS_FOR_SPLITTING` or direct BAM inputs.
        // `pychopper_report` is replaced by `PYCHOPPER_REPORT_CH`
        // `assembly_stats` is replaced by `ASSEMBLY_STATS_CH`

        if (params.transcriptome_source != "precomputed"){
            // Reference-guided assembly path
            // This block will run if BAMS_FOR_SPLITTING has items, regardless of original input type (BAM or FASTQ-derived BAMs)
            // If BAMS_FOR_SPLITTING is empty (e.g. FASTQ input with precomputed transcriptome), these steps will be skipped.

            log.info("Transcriptome source is not precomputed. Attempting transcript assembly.")

            // split_bam now takes BAMS_FOR_SPLITTING
            // Ensure BAMS_FOR_SPLITTING is not empty before calling split_bam if it's a possibility
            SPLIT_BAM_BUNDLES_CH = Channel.empty()
            if (!BAMS_FOR_SPLITTING.isEmpty()) {
                 SPLIT_BAM_BUNDLES_CH = split_bam(BAMS_FOR_SPLITTING).bundles
            } else if (params.transcriptome_source == "reference-guided") {
                // This implies an issue if we expected BAMs (e.g. from FASTQ alignment) but didn't get them
                log.warn "BAMS_FOR_SPLITTING is empty but transcriptome_source is reference-guided. Assembly will be skipped."
            }


            // assemble_transcripts takes output of split_bam
            // It also needs ref_annotation_ch and use_ref_ann_bool
            ASSEMBLED_GFF_BUNDLES_CH = Channel.empty()
            if (!SPLIT_BAM_BUNDLES_CH.isEmpty()) {
                 ASSEMBLED_GFF_BUNDLES_CH = assemble_transcripts(
                    SPLIT_BAM_BUNDLES_CH.flatMap(map_sample_ids_cls).combine(ref_annotation_ch),
                    use_ref_ann_bool
                ).gff_bundles
            } else if (params.transcriptome_source == "reference-guided") {
                 log.warn "Split BAM bundles are empty. Transcript assembly will be skipped."
            }


            MERGED_GFF_CH = Channel.empty()
            TRANSCRIPTOME_SUMMARY_CH = Channel.value(OPTIONAL_FILE) // Initialize as optional

            if (!ASSEMBLED_GFF_BUNDLES_CH.isEmpty()) {
                merge_gff_out = merge_gff_bundles(ASSEMBLED_GFF_BUNDLES_CH.groupTuple())
                MERGED_GFF_CH = merge_gff_out.gff
                TRANSCRIPTOME_SUMMARY_CH = merge_gff_out.summary.map {it[1]}.collect()
            } else if (params.transcriptome_source == "reference-guided") {
                log.warn "Assembled GFF bundles are empty. GFF merging and downstream steps will be skipped."
                // Ensure MERGED_GFF_CH remains empty and TRANSCRIPTOME_SUMMARY_CH is OPTIONAL_FILE
            }


            GFFCOMPARE_DIR_CH = Channel.empty()
            GFFCOMPARE_GTF_CH = Channel.empty() // for merge_transcriptomes
            ISOFORMS_TABLE_CH = Channel.value(OPTIONAL_FILE)
            GFF_TUPLE_FOR_TRANSCRIPTOME_CH = Channel.empty()


            // only run gffcompare if ref annotation provided.
            if (params.ref_annotation && !MERGED_GFF_CH.isEmpty()){
                gffcompare_results = run_gffcompare(MERGED_GFF_CH, ref_annotation_ch)
                GFFCOMPARE_DIR_CH = gffcompare_results.gffcmp_dir
                GFFCOMPARE_GTF_CH = gffcompare_results.gtf // Used by merge_transcriptomes
                ISOFORMS_TABLE_CH = gffcompare_results.isoforms_table.map{ it -> it[1]}.collect()
                // create per sample gff tuples with gff compare directories
                GFF_TUPLE_FOR_TRANSCRIPTOME_CH = MERGED_GFF_CH.join(GFFCOMPARE_DIR_CH)
            } else if (params.ref_annotation && MERGED_GFF_CH.isEmpty() && params.transcriptome_source == "reference-guided") {
                log.warn "Merged GFF is empty, GFFCompare will be skipped."
                ISOFORMS_TABLE_CH = Channel.value(OPTIONAL_FILE) // Ensure it's optional
                GFFCOMPARE_DIR_CH = Channel.value(OPTIONAL_FILE).collect() // makeReport expects a list-like structure for gff_compare
            }
            else if (!params.ref_annotation && !MERGED_GFF_CH.isEmpty()) {
                 // No ref annotation, but have assembled GFFs
                optional_file_ch = Channel.fromPath(OPTIONAL_FILE)
                GFF_TUPLE_FOR_TRANSCRIPTOME_CH = MERGED_GFF_CH.combine(optional_file_ch) // Use OPTIONAL_FILE for gffcompare_dir input to get_transcriptome
                ISOFORMS_TABLE_CH = Channel.value(OPTIONAL_FILE)
                GFFCOMPARE_DIR_CH = Channel.value(OPTIONAL_FILE).collect() // makeReport expects a list-like structure for gff_compare
            } else {
                 // No ref annotation and no merged GFF (e.g. precomputed or failed assembly)
                ISOFORMS_TABLE_CH = Channel.value(OPTIONAL_FILE)
                GFFCOMPARE_DIR_CH = Channel.value(OPTIONAL_FILE).collect()
            }


            // For reference based assembly, there is only one reference
            // So map this reference to all sample_ids
            // Ensure sample_ids is not empty
            TRANSCRIPTOME_FROM_ASSEMBLY_CH = Channel.empty()
            if (!GFF_TUPLE_FOR_TRANSCRIPTOME_CH.isEmpty() && !sample_ids.isEmpty()){
                seq_for_transcriptome_build = sample_ids.flatten().combine(ref_genome_ch)
                TRANSCRIPTOME_FROM_ASSEMBLY_CH = get_transcriptome(
                    GFF_TUPLE_FOR_TRANSCRIPTOME_CH.join(seq_for_transcriptome_build)
                ).transcriptome
            } else if (params.transcriptome_source == "reference-guided") {
                 log.warn "Inputs for get_transcriptome are not available. Transcriptome generation from assembly will be skipped."
            }


            MERGED_GFF_COLLECTED_CH = MERGED_GFF_CH.map{ it -> it[1]}.collect{ it ?: OPTIONAL_FILE }
            // Ensure it emits OPTIONAL_FILE if empty to satisfy makeReport
            if (MERGED_GFF_COLLECTED_CH.isEmpty()){
                MERGED_GFF_COLLECTED_CH = Channel.value(OPTIONAL_FILE)
            }


        } else { // params.transcriptome_source == "precomputed"
            log.info("Transcriptome source is precomputed. Skipping assembly steps.")
            // Set all assembly-derived channels to optional or empty as appropriate
            ISOFORMS_TABLE_CH = Channel.value(OPTIONAL_FILE)
            MERGED_GFF_COLLECTED_CH = Channel.value(OPTIONAL_FILE)
            ASSEMBLY_STATS_CH = Channel.value(OPTIONAL_FILE) // Already default, but reaffirm
            TRANSCRIPTOME_SUMMARY_CH = Channel.value(OPTIONAL_FILE)
            GFFCOMPARE_DIR_CH = Channel.value(OPTIONAL_FILE).collect() // makeReport expects a list-like structure
            // FULL_LEN_READS_FOR_DE should be populated based on input type (FASTQ or BAM)
            // BAMS_FOR_SPLITTING will be empty if precomputed and FASTQ input.
            // If BAM input and precomputed, BAMS_FOR_SPLITTING is populated but assembly is skipped.
        }


        // Differential Expression
        Channel DE_REPORT_CH = Channel.value(OPTIONAL_FILE)
        Channel DE_ALIGNMENT_STATS_CH = Channel.value(OPTIONAL_FILE)
        Channel FINAL_TRANSCRIPTOME_FOR_DE_CH = Channel.empty()
        Channel GTF_FOR_DE_CH = Channel.empty()
        Channel DE_OUTPUTS_CH = Channel.empty()


        if (params.de_analysis){
            sample_sheet_file = file(params.sample_sheet, type:"file")
            // check ref annotation contains only + or - strand as DE analysis will error on .
            check_annotation_strand(ref_annotation_ch).map { stdoutput, annotation ->
            // check if there was an error message
            if (stdoutput) error "In ref_annotation, transcript features must have a strand of either '+' or '-'."
                    stdoutput
                }
            if (!params.ref_transcriptome){ // Transcriptome needs to be generated
                if (params.transcriptome_source == "precomputed") {
                    error "Differential expression (--de_analysis) requires a reference transcriptome (--ref_transcriptome) when --transcriptome_source is 'precomputed'."
                }
                // Use assembly GTF if available, otherwise error (or handle if GTF can be optional for merge_transcriptomes)
                // GFFCOMPARE_GTF_CH should be from run_gffcompare if ref_annotation was used, or could be from other source if assembly done without ref_ann
                // For now, assume run_gffcompare's GTF is the one needed if ref_annotation was provided.
                // If no ref_annotation, merge_transcriptomes might need a different GTF or handle its absence.
                // This part might need refinement based on exact logic for GTF source.
                // For now, we rely on GFFCOMPARE_GTF_CH which is populated if params.ref_annotation and MERGED_GFF_CH is not empty.

                if (GFFCOMPARE_GTF_CH.isEmpty() && params.ref_annotation) {
                     log.warn "GFFCompare GTF is not available, which might be needed for merge_transcriptomes if no ref_transcriptome is provided. DE might fail."
                     // Potentially set FINAL_TRANSCRIPTOME_FOR_DE_CH and GTF_FOR_DE_CH to empty to prevent DE from running
                }

                // Validate ref_annotation against ref_genome
                validate_ref_annotation(ref_annotation_ch, ref_genome_ch).map { stdoutput ->
                    if (stdoutput) {
                        log.warn(stdoutput)
                    }
                }
                // Collect GFFCOMPARE_GTF_CH, ensuring it's a single file or OPTIONAL_FILE
                collected_gtfs_for_merge = GFFCOMPARE_GTF_CH.collect()
                // merge_transcriptomes expects a list of GTFs from query_annotations/*, here we use the collected GTF from gffcompare
                // This assumes GFFCOMPARE_GTF_CH contains the GTFs that would have gone into query_annotations

                merged_trans_results = merge_transcriptomes(collected_gtfs_for_merge, ref_annotation_ch, ref_genome_ch)
                FINAL_TRANSCRIPTOME_FOR_DE_CH = merged_trans_results.fasta
                GTF_FOR_DE_CH = merged_trans_results.gtf
            }
            else { // User-provided reference transcriptome
                FINAL_TRANSCRIPTOME_FOR_DE_CH =  Channel.fromPath(ref_transcriptome_file)
                if (ref_transcriptome_file.extension == "gz") {
                    FINAL_TRANSCRIPTOME_FOR_DE_CH = decompress_transcriptome(FINAL_TRANSCRIPTOME_FOR_DE_CH)
                }
                FINAL_TRANSCRIPTOME_FOR_DE_CH = preprocess_ref_transcriptome(FINAL_TRANSCRIPTOME_FOR_DE_CH)
                GTF_FOR_DE_CH = ref_annotation_ch // Use the main ref_annotation as GTF
            }

            // FULL_LEN_READS_FOR_DE is already prepared:
            // - FASTQ mode: output of preprocess_reads or raw FASTQs if direct_rna
            // - BAM mode: input BAMs mapped to [[alias:id], bam_path]
            // The differential_expression subworkflow needs to handle either FASTQ or BAM.
            if (!FULL_LEN_READS_FOR_DE.isEmpty() && !FINAL_TRANSCRIPTOME_FOR_DE_CH.isEmpty() && !GTF_FOR_DE_CH.isEmpty()) {
                de_workflow_inputs = FULL_LEN_READS_FOR_DE.map{ sample_id_or_meta, reads_path ->
                    // Ensure it's always [[alias:id], path]
                    def meta_map = (sample_id_or_meta instanceof Map) ? sample_id_or_meta : [alias:sample_id_or_meta]
                    [meta_map, reads_path]
                }

                de_results = differential_expression(FINAL_TRANSCRIPTOME_FOR_DE_CH, de_workflow_inputs, sample_sheet_file, GTF_FOR_DE_CH)
                DE_REPORT_CH = de_results.all_de
                DE_OUTPUTS_CH = de_results.de_outputs // Collect this for results
                DE_ALIGNMENT_STATS_CH = de_results.de_alignment_stats
            } else {
                 log.warn "Differential expression skipped due to missing inputs (reads, transcriptome, or GTF)."
                 DE_REPORT_CH = Channel.value(OPTIONAL_FILE)
                 DE_ALIGNMENT_STATS_CH = Channel.value(OPTIONAL_FILE)
                 DE_OUTPUTS_CH = Channel.empty()
            }

        } else{ // Not params.de_analysis
            DE_REPORT_CH = Channel.value(OPTIONAL_FILE)
            DE_ALIGNMENT_STATS_CH = Channel.value(OPTIONAL_FILE)
            DE_OUTPUTS_CH = Channel.empty()
        }


        // Prepare inputs for makeReport
        // `stats` for makeReport comes from INGRESS_RESULTS_CH
        INGRESS_RESULTS_CH.multiMap{ meta, f1, f2, stats ->
            meta: meta
            stats: stats
        }.set { for_report_map }

        metadata_for_report = for_report_map.meta.collect()
        stats_for_report = for_report_map.stats.collect()
        // Ensure stats_for_report is OPTIONAL_FILE if empty, as makeReport expects a path
        if (stats_for_report.isEmpty()) {
            stats_for_report = Channel.value(OPTIONAL_FILE)
        }


        makeReport(
            metadata_for_report,
            stats_for_report,
            software_versions,
            workflow.manifest.version,
            workflow_params,
            DE_ALIGNMENT_STATS_CH, // from de_analysis or OPTIONAL_FILE
            PYCHOPPER_REPORT_CH,   // from FASTQ path or OPTIONAL_FILE
            ASSEMBLY_STATS_CH,     // from FASTQ assembly path or OPTIONAL_FILE
            GFFCOMPARE_DIR_CH.collect{ it ?: OPTIONAL_FILE }, // from assembly path or OPTIONAL_FILE, ensure collected
            MERGED_GFF_COLLECTED_CH, // from assembly path or OPTIONAL_FILE
            DE_REPORT_CH,          // from de_analysis or OPTIONAL_FILE
            ISOFORMS_TABLE_CH,     // from assembly path or OPTIONAL_FILE
            TRANSCRIPTOME_SUMMARY_CH) // from assembly path or OPTIONAL_FILE

       report = makeReport.out.report
       results = results.concat(report) // Initial result is the report itself

       // Collect other results based on what ran
       if (use_ref_ann_bool && params.transcriptome_source != "precomputed"){
            // These results are only expected if reference annotation was used and assembly happened
            if(!GFFCOMPARE_DIR_CH.isEmpty()){
                 results = results.concat(GFFCOMPARE_DIR_CH.map { it[1] }) // Assuming it's [[id, path], ...]
            }
            if(!ASSEMBLY_STATS_CH.isEmpty() && ASSEMBLY_STATS_CH.getVal() != OPTIONAL_FILE){ // Check if it's not the placeholder
                 results = results.concat(ASSEMBLY_STATS_CH)
            }
            if(!ISOFORMS_TABLE_CH.isEmpty() && ISOFORMS_TABLE_CH.getVal() != OPTIONAL_FILE){
                 results = results.concat(ISOFORMS_TABLE_CH)
            }
            if(!TRANSCRIPTOME_FROM_ASSEMBLY_CH.isEmpty()){
                 results = results.concat(TRANSCRIPTOME_FROM_ASSEMBLY_CH.flatMap(map_sample_ids_cls).map {it -> it[1]})
            }
       } else if (!use_ref_ann_bool && params.transcriptome_source == "reference-guided"){
            // Results if no ref annotation but reference-guided assembly
            if(!ASSEMBLY_STATS_CH.isEmpty() && ASSEMBLY_STATS_CH.getVal() != OPTIONAL_FILE){
                 results = results.concat(ASSEMBLY_STATS_CH)
            }
            if(!TRANSCRIPTOME_FROM_ASSEMBLY_CH.isEmpty()){
                results = results.concat(TRANSCRIPTOME_FROM_ASSEMBLY_CH.flatMap(map_sample_ids_cls).map {it -> it[1]})
            }
       }

       // Add input files collected earlier (harmonized across BAM/FASTQ)
       results = results.map{ [it, null] }.concat(input_files_for_report.map { [it, "input_files"] })


        if (params.de_analysis){
           de_results_to_publish = report.concat(
            FINAL_TRANSCRIPTOME_FOR_DE_CH.ifEmpty(Channel.empty()), // ensure it doesn't fail if empty
            DE_OUTPUTS_CH.flatten().ifEmpty(Channel.empty()),
            makeReport.out.results_dge.ifEmpty(Channel.empty()),
            makeReport.out.tpm.ifEmpty(Channel.empty()),
            makeReport.out.filtered.ifEmpty(Channel.empty()),
            makeReport.out.unfiltered.ifEmpty(Channel.empty()),
            makeReport.out.gene_counts.ifEmpty(Channel.empty()))
            // Output de_analysis results in the dedicated directory.
            results = results.concat(de_results_to_publish.map{ [it, "de_analysis"] })
        }

        results = results.concat(workflow_params.map{ [it, null]})

        // IGV config
        // This part needs careful review if BAMS_FOR_SPLITTING is the source of BAMs for IGV
        // or if it should use assembly.bam from the reference_assembly process when run.
        // For now, assume BAMS_FOR_SPLITTING can be used if it contains the relevant BAMs.
        // The original logic used `reads` channel which was the input fastq/bam.
        // The `publish_prefix_bams` is "BAMS". IGV needs `sample_id_reads_aln_sorted.bam`
        // If BAMS_FOR_SPLITTING is from xam_ingress, they are likely already named correctly by sample_id.
        // If BAMS_FOR_SPLITTING is from reference_assembly, they are also named by sample_id.

        if (params.transcriptome_source == "precomputed" && params.igv){
            log.warn("IGV configuration might not be fully functional if transcriptome source is 'precomputed' and BAMs are not directly provided or generated in a standard way for IGV.")
        }

        if (params.igv && params.ref_genome && (is_fastq_input || is_bam_input || params.bam) ){ // Ensure there's a ref and some form of input
            // Determine source of BAMs for IGV:
            // If assembly ran (FASTQ -> BAM), use assembly.bam.
            // If BAM input, use those BAMs.
            // This logic needs to be robust. For now, let's assume BAMS_FOR_SPLITTING holds the necessary BAMs.
            // However, original IGV logic used `reads | map { meta, sample, index, stats -> meta.alias }`
            // which implies it needs the initial ingress channel structure.
            // We'll use INGRESS_RESULTS_CH for meta and assume BAMs are in publish_prefix_bams.

            is_compressed = ref_genome_file.extension == "gz"
            String publish_ref_dir = "igv_reference" // Directory where reference and indexes will be published by faidx/gz_faidx

            // Reference genome path for IGV (local path within container/execution)
            igv_ref_uri_ch = ref_genome_ch | flatten | map { file(it).toUriString() }


            Channel igv_fai_idx_ch = Channel.empty()
            Channel igv_gzi_idx_ch = Channel.empty() // For .gz.gzi

            if (is_compressed){
                String input_fai_index = "${ref_genome_file}.fai"
                String input_gzi_index = "${ref_genome_file}.gzi"
                if (file(input_fai_index).exists() && file(input_gzi_index).exists()){
                    igv_fai_idx_ch = Channel.fromPath(input_fai_index) | map { file(it).toUriString() }
                    igv_gzi_idx_ch = Channel.fromPath(input_gzi_index) | map { file(it).toUriString() }
                } else {
                    gz_faidx_out = gz_faidx(ref_genome_ch) // Use the channel for ref_genome
                    igv_fai_idx_ch = gz_faidx_out.map { fai, gzi -> "$publish_ref_dir/${file(fai).getName()}" }
                    igv_gzi_idx_ch = gz_faidx_out.map { fai, gzi -> "$publish_ref_dir/${file(gzi).getName()}" }
                    igv_fai_idx_ch | ifEmpty{ log.warn "Failed to generate .fai for compressed ref; IGV might be affected." }
                    igv_gzi_idx_ch | ifEmpty{ log.warn "Failed to generate .gzi for compressed ref; IGV might be affected." }
                }
            } else { // Not compressed
                if (file("${ref_genome_file}.fai").exists()){
                    igv_fai_idx_ch = Channel.fromPath("${ref_genome_file}.fai") | map { file(it).toUriString() }
                } else {
                    fai_out = faidx(ref_genome_ch) // Use the channel for ref_genome
                    igv_fai_idx_ch = fai_out.map { fai -> "$publish_ref_dir/${file(fai).getName()}" }
                    igv_fai_idx_ch | ifEmpty{ log.warn "Failed to generate .fai for uncompressed ref; IGV might be affected." }
                }
            }


            // Get list of sample aliases for BAM paths
            // Assuming BAMs are published to `params.out_dir + / + publish_prefix_bams / + {sample_id}_reads_aln_sorted.bam`
            // This was the pattern in reference_assembly. If input BAMs are used, their names need to match this
            // or the paths need to be constructed from BAMS_FOR_SPLITTING.
            // For simplicity, we use sample_ids and assume the published path structure.
            igv_bam_paths_ch = sample_ids
            | toSortedList
            | map { list -> list.collect{ id ->
                [
                 // This path needs to be the final path where sorted indexed BAMs are expected for IGV
                 // If direct BAM input, these might not be `_reads_aln_sorted.bam` unless explicitly named so.
                 // This part is tricky if input BAMs are not processed by `reference_assembly`.
                 // For now, assume this naming convention is achieved by upstream processes or direct input naming.
                 "$publish_prefix_bams/${id}_reads_aln_sorted.bam",
                 "$publish_prefix_bams/${id}_reads_aln_sorted.bam.bai"
                ]
            } }
            | flatten


            igv_all_files_ch = igv_bam_paths_ch
            | concat (igv_ref_uri_ch.ifEmpty(Channel.empty()))
            | concat (igv_fai_idx_ch.ifEmpty(Channel.empty()))
            | concat (igv_gzi_idx_ch.ifEmpty(Channel.empty()))
            | flatten
            | collectFile(name: "igv-file-names.txt", newLine: true, sort: false)


            igv_conf_ch = configure_igv(
                igv_all_files_ch,
                Channel.of(null), // igv locus
                [displayMode: "SQUISHED", colorBy: "strand"], // bam extra opts
                Channel.of(null), // vcf extra opts
                )
            results = results.concat(igv_conf_ch.map{ [it, null]})
      } else if (params.igv) {
            log.warn "IGV session generation skipped due to missing reference genome or input files."
      }

    emit:
        results
}

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

    error = null

    // --- Parameter Validation ---
    if (params.bam && params.fastq) {
        error = "Parameters --bam and --fastq are mutually exclusive. Please specify only one."
    }
    if (params.bam_input_dir && params.fastq) {
        error = "Parameters --bam_input_dir and --fastq are mutually exclusive. Please specify only one."
    }
    if (params.bam_input_dir && params.bam) {
        error = "Parameters --bam_input_dir and --bam are mutually exclusive. Use --bam_input_dir with --sample_sheet."
    }
    if ((params.bam_input_dir || params.fastq) && !params.sample_sheet) {
        error = "A sample sheet (--sample_sheet) is required when using --bam_input_dir or --fastq."
    }
    if (!params.bam && !params.fastq && !params.bam_input_dir) {
        error = "No input specified. Please provide --bam, or --fastq with --sample_sheet, or --bam_input_dir with --sample_sheet."
    }


    if (params.containsValue("jaffal_refBase")) {
        error = "JAFFAL fusion detection has been removed from this workflow."
    }
    if (params.containsKey("minimap_index_opts")) {
        error = "`--minimap_index_opts` parameter is deprecated. Use parameter `--minimap2_index_opts` instead."
    }

    if (params.transcriptome_source == "precomputed" && !params.ref_transcriptome){
        error = "As transcriptome source parameter is precomputed you must include a ref_transcriptome parameter"
    }
    // ref_genome is now optional if transcriptome_source is precomputed.
    // However, if de_analysis is true and no ref_transcriptome is given, then ref_genome is needed for merge_transcriptomes.
    // And if igv is true, ref_genome is needed.
    if (params.transcriptome_source == "reference-guided" && !params.ref_genome){
        error = "As transcriptome source is reference guided you must include a ref_genome parameter"
    }
    if (params.de_analysis && !params.ref_transcriptome && !params.ref_genome){
        error = "For de_analysis without a precomputed ref_transcriptome, a ref_genome is required to build one."
    }
    if (params.igv && !params.ref_genome){
        error = "IGV session generation requires a reference genome (--ref_genome)."
    }


    Channel ref_genome_file_ch
    if (params.ref_genome){
        ref_genome_file_ch = file(params.ref_genome, type: "file", checkIfExists: true)
    }else {
        ref_genome_file_ch = OPTIONAL_FILE // Pass as optional file if not provided
    }

    if (params.containsValue("denovo")) {
        error = "Denovo transcriptome source is no longer supported. Please use the reference-guided or precomputed options."
    }

    boolean use_ref_ann_flag = false
    Channel ref_annotation_file_ch
    if (params.ref_annotation){
        ref_annotation_file_ch = file(params.ref_annotation, type: "file", checkIfExists: true)
        use_ref_ann_flag = true
    }else{
        ref_annotation_file_ch= OPTIONAL_FILE
        use_ref_ann_flag = false
    }

    Channel ref_transcriptome_file_ch = OPTIONAL_FILE
    if (params.ref_transcriptome){
        log.info("Reference Transcriptome provided will be used for differential expression.")
        ref_transcriptome_file_ch = file(params.ref_transcriptome, type:"file", checkIfExists: true)
    }

    if (params.de_analysis){
        if (!params.ref_annotation){
            error = "When running in --de_analysis mode you must provide a reference annotation (--ref_annotation)."
        }
        if (!params.sample_sheet){
            error = "You must provide a sample_sheet (--sample_sheet) with at least alias and condition columns for DE analysis."
        }
        if (params.containsKey("condition_sheet")) {
        error = "Condition sheets have been deprecated. Please add a 'condition' column to your sample sheet instead. Check the quickstart for more information."
        }
    } else{
        if (!params.ref_annotation && params.transcriptome_source != "precomputed"){
            log.info("Warning: As no --ref_annotation was provided and transcriptome_source is not 'precomputed', the output transcripts will not be annotated against a reference annotation.")
        }
    }
    if (error){
        throw new Exception(error)
    }

    // Input data ingress
    // This `samples_ch` will be passed to the main pipeline workflow
    // It should consistently provide [meta, file1, file2_or_null, stats_file_or_null]
    // file1: fastq or bam
    // file2_or_null: bam_index (bai) if file1 is bam
    // stats_file_or_null: stats from ingress
    Channel samples_ch

    if (params.fastq) { // FASTQ input with sample sheet
        samples_ch = fastq_ingress([
            "input":params.fastq,
            "sample":params.sample, // sample param might be redundant if sheet is comprehensive
            "sample_sheet":params.sample_sheet,
            "analyse_unclassified":params.analyse_unclassified,
            "stats": true,
            "fastcat_extra_args": "",
            "per_read_stats": true])
        // fastq_ingress emits: [meta, concat_seqs, fastcat_stats]
        // We need to adapt it to the [meta, file1, null, stats] structure for consistency
        samples_ch = samples_ch.map { meta, concat_seqs, fastcat_stats ->
            [meta, concat_seqs, null, fastcat_stats]
        }
    } else if (params.bam_input_dir) { // BAM input with sample sheet
        samples_ch = xam_ingress([
            "input":params.bam_input_dir, // xam_ingress needs to handle directory input based on sample sheet
            "sample":params.sample,
            "sample_sheet":params.sample_sheet,
            "analyse_unclassified":params.analyse_unclassified,
            "keep_unaligned": true, // Or make this configurable
            "return_fastq": false,  // Ensure BAMs are output
            "stats": true,
            "per_read_stats": true])
        // xam_ingress emits: [meta, bam, bai, stats] - this matches the target structure.
    } else if (params.bam) { // Single BAM file input
         samples_ch = xam_ingress([
            "input":params.bam, // Single BAM path
            "sample":params.sample, // Optional sample name for this single BAM
            "sample_sheet":null, // No sample sheet for single BAM
            "analyse_unclassified":params.analyse_unclassified,
            "keep_unaligned": true,
            "return_fastq": false,
            "stats": true,
            "per_read_stats": true])
        // xam_ingress emits: [meta, bam, bai, stats]
    }


    pipeline(samples_ch, ref_genome_file_ch, ref_annotation_file_ch, ref_transcriptome_file_ch, use_ref_ann_flag)
    publish_results(pipeline.out.results)
}

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
                validate_ref_annotation(ref_annotation, ref_genome).map { stdoutput ->
                    if (stdoutput) {
                        log.warn(stdoutput)
                    }
                }
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
            de = differential_expression(transcriptome, full_len_reads.map{ sample_id, reads -> [[alias:sample_id], reads]}, sample_sheet, gtf)
            de_report = de.all_de
            de_outputs = de.de_outputs
            de_alignment_stats = de.de_alignment_stats
        } else{
            de_report = OPTIONAL_FILE
            de_alignment_stats = OPTIONAL_FILE
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
            de_alignment_stats,
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
                 "$publish_prefix_bams/${it}_reads_aln_sorted.bam",
                 "$publish_prefix_bams/${it}_reads_aln_sorted.bam.bai"
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
