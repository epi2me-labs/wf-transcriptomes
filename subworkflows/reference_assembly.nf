process map_reads{
    /*
    Map reads to reference using minimap2.
    Filter reads by mapping quality.
    Filter internally-primed reads.
    */
    label "isoforms"
    cpus params.threads

    //added (AS 29v2023)
    publishDir params.mappedOut, mode:'copy'



    input:
       path index
       path reference
       tuple val(sample_id), path (fastq_reads)

    output:
       tuple val(sample_id), path("${sample_id}_reads_aln_sorted.bam"), emit: bam
       tuple val(sample_id), path("${sample_id}_read_aln_stats.tsv"), emit: stats
    script:
        def ContextFilter = """AlnContext: { Ref: "${reference}", LeftShift: -${params.poly_context},
        RightShift: ${params.poly_context}, RegexEnd: "[Aa]{${params.max_poly_run},}",
        Stranded: True,Invert: True, Tsv: "internal_priming_fail.tsv"} """
    """
    minimap2 -t ${params.threads} -ax splice ${params.minimap2_opts} ${index} ${fastq_reads}\
        | samtools view -q ${params.minimum_mapping_quality} -F 2304 -Sb -\
        | seqkit bam -j ${params.threads} -x -T '${ContextFilter}' -\
        | samtools sort -@ ${params.threads} -o "${sample_id}_reads_aln_sorted.bam" - ;
    ((cat "${sample_id}_reads_aln_sorted.bam" | seqkit bam -s -j ${params.threads} - 2>&1)  | tee ${sample_id}_read_aln_stats.tsv ) || true

    if [[ -s "internal_priming_fail.tsv" ]];
        then
            tail -n +2 "internal_priming_fail.tsv" | awk '{print ">" \$1 "\\n" \$4 }' - > "context_internal_priming_fail_start.fasta"
            tail -n +2 "internal_priming_fail.tsv" | awk '{print ">" \$1 "\\n" \$6 }' - > "context_internal_priming_fail_end.fasta"
    fi
    """
}

workflow reference_assembly {
    take:
       index
       reference
       fastq_reads
    main:
        map_reads(index, reference, fastq_reads)
    emit:
       bam = map_reads.out.bam
       stats = map_reads.out.stats
}
