// This subworkflow pre-processes an input reference
// genome by preparing the appropriate indexes and cache.
// By default, the workflow takes an input reference,
// decompresses it if it is gzipped, and either imports
// an existing index (if available) or generates one
// using `samtools faidx`.
// If requested by the user, the workflow can also:
// 1. Generate a CRAM cache, using the `REF_PATH`
//    environment variable. IMPORTANT: If the ingress workflow
//    outputs CRAM files (when output_xam_fmt='cram'), a CRAM cache
//    MUST be generated here and passed to downstream processes.
//    Even though users cannot currently input CRAMs, the ingress
//    workflow can generate them, and downstream processes require
//    the cache to decompress CRAM files.
// 2. Generate a minimap2 `.mmi` index for faster alignment.
Map colors = NfcoreTemplate.logColours(params.monochrome_logs)

// Argument parser
Map parse_reference(Map arguments) {
    ArgumentParser parser = new ArgumentParser(
        args:[
        ],
        kwargs:[
            "output_cache": false,
            "output_mmi": false,
            "mmi_opts": "-x lr:hq",
            "mmi_memory": ["8GB", "15GB", "31GB"],
        ],
        name: "reference_ingress")
    return parser.parse_args(arguments)
}

// Process to generate the CRAM cache and
// create the REF_PATH variable.
// NOTE: The CRAM cache is essential for any downstream processes that
// need to decompress or work with CRAM files. Even though users cannot
// currently input CRAMs, if output_xam_fmt='cram' is used in ingress,
// downstream processes will require this cache to properly handle the
// CRAM output from ingress.
process cram_cache {
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        path reference
    output:
        tuple path("ref_cache/"), env(REF_PATH), emit: ref_cache
    script:
    """
    # Invoke from binary installed to container PATH
    seq_cache_populate.pl -root ref_cache/ "${reference}"
    REF_PATH="ref_cache/%2s/%2s/%s"
    """
}

// Process to create the faidx index
// main.nf currently handles publishing the faidx if required
// maintain naming so it matches input ref when output
process faidx {
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        path(ref)
    output:
        path("${ref}.fai")
    script:
    """
    samtools faidx "${ref}"
    """
}

// Process to create the faidx indexes for a gzipped reference
process gz_faidx {
    label "wf_common"
    cpus 1
    memory 4.GB
    // If a user provides a non-bgzipped file indexing will fail, which should be tolerated.
    // With error_strategy set to 'ignore', this will result in an empty channel 
    // that can be checked for downstream using `ifEmpty()` 
    errorStrategy 'ignore'
    input:
        path(ref)
    output:
        tuple path(ref), path("${ref}.fai"), path("${ref}.gzi")
    script:
    """
    samtools faidx "${ref}"
    """
}

// Decompress the reference genome
// NOTE -f required to decompress symlink
process decompress_ref {
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        path "ref.fa.gz"
    output:
        path "ref.fa", emit: decompressed_ref
    """
    gzip -df ref.fa.gz
    """
}

// Prepare minimap2 .mmi index
process make_mmi {
    memory {  margs["mmi_memory"][task.attempt - 1] }
    maxRetries { margs["mmi_memory"].size() - 1 }
    errorStrategy 'retry'
    label "wf_common"
    cpus 3
    input:
        path("ref.fa")
        val margs
    output:
        path("ref.mmi")
    script:
    """
    minimap2 -t ${task.cpus} ${margs["mmi_opts"]} -d ref.mmi ref.fa
    """
}


// Workflow to prepare the reference genome and its indexes.
/**
 * This workflow accepts either a string file path or a Nextflow channel as the reference input,
 * making it flexible for both direct file inputs and outputs from upstream processes.
 * 
 * IMPORTANT: If your workflow uses CRAM output from ingress (output_xam_fmt='cram'),
 * you must set "output_cache": true to generate the CRAM cache. This cache is required
 * by downstream processes to properly handle CRAM files, even though users cannot currently
 * input CRAMs directly.
 * 
 * inputs:
 *  - reference_path_or_channel: (String or Channel) Path to reference FASTA file or a Nextflow channel emitting a reference file
 *  - arguments: (Map) Optional configuration parameters:
 *    - "output_cache": (boolean) If true, generate a CRAM cache for the reference (default: false)
 *                      REQUIRED if using CRAM output from ingress for downstream processing.
 *    - "output_mmi": (boolean) If true, generate a minimap2 .mmi index for the reference (default: false)
 *    - "mmi_opts": (string) Options string passed to minimap2 for index creation (default: "-x lr:hq")
 *    - "mmi_memory": (list) Memory options to pass to minimap2 index; creation will be retried with next in list if it fails due to memory issues (default: ["8GB", "15GB", "31GB"])
 *
 * emits:
 *  - "ref_tuple": Tuple of (reference, fai_index)
 *  - "ref_cache": CRAM cache for the reference genome if "output_cache" is true.
 *                 Required for downstream processes when CRAM files are used.
 *  - "ref_mmi": Minimap2 index for the reference genome if "output_mmi" is true
 *  - "ref_gzidx": Gzip index for the reference genome (useful for IGV)
 */
workflow prepare_reference {
    take:
        reference_path_or_channel
        arguments
    main:
        Map margs = parse_reference(arguments)
        
        // Determine if input is a string path or channel
        Boolean is_string = reference_path_or_channel in String
        
        // Base ref channel
        if (!is_string) {
            // It's a channel - use it directly
            ref = reference_path_or_channel
            is_compressed = false  // Assume uncompressed when coming from channel
            gzindexes = Channel.empty()
        } else {
            // It's a string path
            if (file(reference_path_or_channel).exists()){
                ref = Channel.fromPath(reference_path_or_channel)
            } else {
                throw new Exception(colors.red + "File ${reference_path_or_channel} not found." + colors.reset)
            }
            
            // Check if compressed
            is_compressed = reference_path_or_channel.toLowerCase().endsWith("gz")
            
            if (is_compressed) {
                // Define indexes names.
                String input_fai_index = "${reference_path_or_channel}.fai"
                String input_gzi_index = "${reference_path_or_channel}.gzi"

                // Decompress reference genome.
                ref = decompress_ref(ref)
                // Check whether the input gzref is indexed. If so, pass these as indexes.
                // Otherwise, generate the gzip + fai indexes for the compressed reference.
                if (file(input_fai_index).exists() && file(input_gzi_index).exists()){
                    gzindexes = Channel.fromPath(reference_path_or_channel)
                    | mix(
                        Channel.fromPath(input_fai_index),
                        Channel.fromPath(input_gzi_index)
                    )
                } else {
                    gzindexes = gz_faidx(Channel.fromPath(reference_path_or_channel))
                    gzindexes.ifEmpty{
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
            }
        }

        // Generate fai index if the file is either compressed, or if fai doesn't exists
        if (!is_compressed) {
            if (is_string && file("${reference_path_or_channel}.fai").exists()){
                ref_idx = Channel.fromPath("${reference_path_or_channel}.fai")
            } else {
                ref_idx = faidx(ref)
            }
        } else {
            ref_idx = faidx(ref)
        }

        // Combine ref and its fai index
        ref_tuple = ref.combine(ref_idx)

        // Generate CRAM cache
        if (margs.output_cache){
            cram_cache(ref)
            ref_cache = cram_cache.out.ref_cache
        } else {
            ref_cache = null
        }

        // Generate mmi index
        if (margs.output_mmi){
            ref_mmi = make_mmi(ref, margs)
        } else {
            ref_mmi = null
        }

    emit:
        ref_tuple = ref_tuple
        ref_cache = ref_cache
        ref_mmi = ref_mmi
        ref_gzidx = gzindexes
}
