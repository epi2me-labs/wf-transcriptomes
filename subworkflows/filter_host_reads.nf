/*
Subworkflow for filtering host reads from xenograft samples
Map to both genomes separateply (use previously mapped reads for graft - product of subworkflow map_reads_all_genome)
Filter using XenofileR
bam to fastq using samtools
*/


params.filteredFastqOut="${params.out_dir}/fastq_filtered_graft"
params.scriptDir="${projectDir}/bin"


process build_minimap_index_host{
    /*
    Build minimap index from host reference genome
    */
    label "isoforms"
    cpus params.threads

    input:
        path reference_host
    output:
        path "genome_host_index.mmi", emit: index_host
    script:
    """
    minimap2 -t ${params.threads} ${params.minimap_index_opts} -I 1000G -d "genome_host_index.mmi" ${reference_host}
    """
}



process map_reads_unfilt_host{
    /*
    Map reads to host reference using minimap2.
    Filter reads by mapping quality.
    Filter internally-primed reads.
    */
    label "isoforms"
    cpus params.threads

    input:
       path index_host
       tuple val(sample_id), path (full_len_reads_host_graft)

    output:
       tuple val(sample_id), path("${sample_id}_all_alns.host.minimap2.bam"), emit: bam_host

    script:
    """
    minimap2 -t ${params.threads} -ax splice -uf ${index_host} ${full_len_reads_host_graft}\
        | samtools view -hbo -\
        | samtools sort -@ ${params.threads} -o "${sample_id}_all_alns.host.minimap2.bam" - 

    """
}

process map_reads_unfilt_graft{
    /*
    Map reads to graft reference using minimap2.
    Filter reads by mapping quality.
    Filter internally-primed reads.
    */
    label "isoforms"
    cpus params.threads

    input:
       path index
       tuple val(sample_id), path(full_len_reads_host_graft)

    output:
       tuple val(sample_id), path("${sample_id}_all_alns.graft.minimap2.bam"), emit: bam_graft

    script:
    """
    minimap2 -t ${params.threads} -ax splice -uf ${index} ${full_len_reads_host_graft}\
        | samtools view -hbo -\
        | samtools sort -@ ${params.threads} -o "${sample_id}_all_alns.graft.minimap2.bam" - 

    """
}



process filter_host_reads_bam{
    /*
    Filter host reads using XenofilteR
    */

	label "xenofilter"
    cpus 8

    input:
    tuple val(sample_id), path(mapped_graft_bam)
    tuple val(sample_id), path(mapped_host_bam)

    output:
    tuple val(sample_id), path("Filtered_bams/${sample_id}_Filtered.bam"), path("Filtered_bams/${sample_id}_Filtered.bam.bai"), emit: bam_filtered_graft

    script:
    """
    Rscript $params.scriptDir/run_xenofilter.R ${sample_id} ${mapped_graft_bam} ${mapped_host_bam}
    """

}


process convert_graft_reads{
    /*
    bam to fastq conversion
    */

    label "isoforms"
    cpus params.threads

    publishDir params.filteredFastqOut, mode:'copy'

    input:
    tuple val(sample_id), path(bam_filtered_graft), path(bam_filtered_graft_bai)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.graft.fastq"), emit: fastq_graft
    tuple val(sample_id), path("${sample_id}.host_filtering_stats.txt"), emit: stats_fastq_filt

    script:
    """
    samtools fastq -F 0x4 ${bam_filtered_graft} > ${sample_id}.filtered.graft.fastq
    wc -l ${sample_id}.filtered.graft.fastq >${sample_id}.host_filtering_stats.txt

    """

}



workflow filter_host_reads {
    take:
       index
       fastq_reads
       reference_host
    main:
    	build_minimap_index_host(reference_host)
    	index_host = build_minimap_index_host.out.index_host
        map_reads_unfilt_host(index_host, fastq_reads)
        map_reads_unfilt_graft(index, fastq_reads)

        filter_host_reads_bam(map_reads_unfilt_graft.out.bam_graft, map_reads_unfilt_host.out.bam_host)
        convert_graft_reads(filter_host_reads_bam.out.bam_filtered_graft)

    emit:
       fastq_graft = convert_graft_reads.out.fastq_graft
       stats_filt = convert_graft_reads.out.stats_fastq_filt

}




