params.mappedAllOut="${params.out_dir}/bam_minimap_genome_all"


process map_reads_all{
    /*
    Map reads to reference using minimap2.
    Filter reads by mapping quality.
    Filter internally-primed reads.
    */
    label "isoforms"
    cpus params.threads


    //added (AS 29v2023)
    publishDir params.mappedAllOut, mode:'copy'



    input:
       path index
       path reference
       tuple val(sample_id), path(fastq_reads)

    output:
       tuple val(sample_id), path("${sample_id}_all_alns.minimap2.bam"), emit: bam
       tuple val(sample_id), path("${sample_id}_all_alns.minimap2.mapped.txt"), emit: stats_map
       tuple val(sample_id), path("${sample_id}_all_alns.minimap2.unmapped.txt"), emit: stats_unmap
       tuple val(sample_id), path("${sample_id}_all_alns.minimap2.mapstats.txt"), emit: stats_mapstats

    script:
    """
    minimap2 -t ${params.threads} -ax splice -uf ${index} ${fastq_reads}\
        | samtools view -hbo -\
        | samtools sort -@ ${params.threads} -o "${sample_id}_all_alns.minimap2.bam" - 

    samtools stats -f 4 ${sample_id}_all_alns.minimap2.bam > ${sample_id}_all_alns.minimap2.unmapped.txt

    samtools stats -F 4 ${sample_id}_all_alns.minimap2.bam > ${sample_id}_all_alns.minimap2.mapped.txt

    echo "alignment statistics for file ${sample_id}_all_alns.minimap2.bam" >>${sample_id}_all_alns.minimap2.mapstats.txt
    echo "all primary alignments -F 0x904" >>${sample_id}_all_alns.minimap2.mapstats.txt
    samtools view -F 0x904 ${sample_id}_all_alns.minimap2.bam | wc -l >>${sample_id}_all_alns.minimap2.mapstats.txt

    echo "unmapped -f 4" >>${sample_id}_all_alns.minimap2.mapstats.txt
    samtools view -f 4 ${sample_id}_all_alns.minimap2.bam | wc -l >>${sample_id}_all_alns.minimap2.mapstats.txt

    echo "all primary alignments at MAPQ 60 -F 0x904 -q 60" >>${sample_id}_all_alns.minimap2.mapstats.txt
    samtools view -F 0x904 -q 60 ${sample_id}_all_alns.minimap2.bam | wc -l >>${sample_id}_all_alns.minimap2.mapstats.txt

    """
}





workflow map_reads_all_genome {
    take:
       index
       reference
       fastq_reads
    main:
        map_reads_all(index, reference, fastq_reads)
    emit:
       bam = map_reads_all.out.bam
       stats_map = map_reads_all.out.stats_map
       stats_unmap = map_reads_all.out.stats_unmap
       stats_mapstats = map_reads_all.out.stats_mapstats
}


