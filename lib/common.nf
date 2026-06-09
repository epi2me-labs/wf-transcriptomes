import groovy.json.JsonBuilder

process getParams {
    label "wf_common"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "params.json"
    cache false
    cpus 1
    memory "2 GB"
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString().replaceAll("'", "'\\\\''")
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process configure_igv {
    publishDir "${params.out_dir}/", mode: 'copy', pattern: 'igv.json', enabled: params.containsKey("igv") && params.igv
    label "wf_common"
    cpus 1
    memory "2 GB"
    input:
        // the python script will work out what to do with all the files based on their
        // extensions
        path "file-names.txt"
        val locus_str
        val aln_extra_opts
        val var_extra_opts
        val keep_track_order
    output: path "igv.json"
    script:
    // the locus argument just makes sure that the initial view in IGV shows something
    // interesting
    String locus_arg = locus_str ? "--locus '$locus_str'" : ""
    // extra options for alignment tracks
    def aln_opts_json_str = \
        aln_extra_opts ? new JsonBuilder(aln_extra_opts).toPrettyString() : ""
    String aln_extra_opts_arg = \
        aln_extra_opts ? "--extra-alignment-opts extra-aln-opts.json" : ""
    // extra options for variant tracks
    def var_opts_json_str = \
        var_extra_opts ? new JsonBuilder(var_extra_opts).toPrettyString() : ""
    String var_extra_opts_arg = \
        var_extra_opts ? "--extra-vcf-opts extra-var-opts.json" : ""
    String keep_track_order_arg = \
        keep_track_order ? "--keep-track-order" : ""
    """
    # write out JSON files with extra options for the alignment and variant tracks
    echo '$aln_opts_json_str' > extra-aln-opts.json
    echo '$var_opts_json_str' > extra-var-opts.json

    workflow-glue configure_igv \
        --fofn file-names.txt \
        $locus_arg \
        $aln_extra_opts_arg \
        $var_extra_opts_arg \
        $keep_track_order_arg \
    > igv.json
    """
}


process minimap2_alignment {
    label "wf_common"
    cpus { bsargs["alignment_threads"] }
    memory { bsargs["minimap2_memory"][task.attempt - 1] }
    maxRetries { bsargs["minimap2_memory"].size() - 1 }
    errorStrategy "retry"  // retry on any error as we may not exit with appropriate codes when exceeding memory
    input:
        tuple val(meta), path(reads, stageAs: 'reads/*'), path(reference), path(ref_index)
        tuple val(align_ext), val(index_ext) // either [bam, bai] or [cram, crai]
        val bsargs
    
    output:
        tuple val(meta),
              path("reads.${align_ext}"), 
              path("reads.${align_ext}.${index_ext}"),
              path("bamstats_results"),
              emit: alignment
    script:
        String reset_cmd_body = "samtools reset -x tp,cm,s1,s2,NM,MD,AS,SA,ms,nn,ts,cg,cs,dv,de,rl,zd"
        String fastq_cmd_body = "samtools fastq -T '*'"
        // Default required threads is 6
        // Samtools x3 and bamstats will all be single threaded
        def minimap2_threads = Math.max(task.cpus - 4 , 1)
        def per_read_stats_arg = bsargs["per_read_stats"] ? "| bgzip > bamstats_results/bamstats.readstats.tsv.gz" : " > /dev/null"
        def bam_input_cmd = """samtools view -H --no-PG reads/\"\$(ls reads | head -n 1)\" \
            > reads.header && samtools cat reads/* \
            | ${reset_cmd_body} --no-PG - -o -  \
            | ${fastq_cmd_body} - """
        def fastq_input_cmd = "touch reads.header && fastcat --reheader reads/*"
        def bam_or_fastq_input = params.bam ? bam_input_cmd : fastq_input_cmd
        def minimap2_opts = bsargs["minimap2_opts"]
        def sort_threads = bsargs["sort_threads"]
    """
    rm -rf bamstats_results
    mkdir bamstats_results
    ${bam_or_fastq_input} \
        | minimap2 -y -t ${minimap2_threads} -a ${minimap2_opts} --cap-kalloc 100m --cap-sw-mem 50m \
            ${reference} - \
        | workflow-glue reheader_samstream reads.header \
            --insert \$'@PG\\tID:reset\\tPN:samtools\\tCL:${reset_cmd_body}' \
            --insert \$'@PG\\tID:fastq\\tPN:samtools\\tCL:${fastq_cmd_body}' \
        | samtools sort -@ ${sort_threads} -u -O BAM - \
        | tee >(samtools view --reference ${reference} \
            --write-index -o reads.${align_ext}##idx##reads.${align_ext}.${index_ext} -) \
        | bamstats -s ${meta.alias} -u \
            -f bamstats_results/bamstats.flagstat.tsv  \
            -i bamstats_results/bamstats.runids.tsv \
            -l bamstats_results/bamstats.basecallers.tsv \
            --histograms histograms - \
            ${per_read_stats_arg}

        # post-pipe bamstats tidying
        mv histograms/* bamstats_results/

        # get n_seqs from flagstats - need to sum them up
        awk 'NR==1{for (i=1; i<=NF; i++) {ix[\$i] = i}} NR>1 {c+=\$ix["total"]} END{print c}' \
            bamstats_results/bamstats.flagstat.tsv > bamstats_results/n_seqs
       
        # get unique run IDs (we add `-F '\\t'` as `awk` uses any stretch of whitespace
        # as field delimiter otherwise and thus ignore empty columns)
        awk -F '\\t' '
            NR==1 {for (i=1; i<=NF; i++) {ix[\$i] = i}}
            # only print run_id if present
            NR>1 && \$ix["run_id"] != "" {print \$ix["run_id"]}
        ' bamstats_results/bamstats.runids.tsv | sort | uniq > bamstats_results/run_ids
        
        # get unique basecall models
        awk -F '\\t' '
           NR==1 {for (i=1; i<=NF; i++) {ix[\$i] = i}}
           # only print model if present
           NR>1 && \$ix["basecaller"] != "" {print \$ix["basecaller"]}
        ' bamstats_results/bamstats.basecallers.tsv | sort | uniq > bamstats_results/basecallers ;

    """
}
