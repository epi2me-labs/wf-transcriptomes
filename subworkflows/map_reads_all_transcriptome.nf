params.mappedAllTrxOut="${params.out_dir}/bam_minimap_transcriptome_filt"
params.salmonOut="${params.out_dir}/salmon"


process map_reads_trx_all{
    /*
    Map reads to reference using minimap2.
    Filter reads by mapping quality.
    */
    label "isoforms"
    cpus params.threads


    //added (AS 29v2023)
    publishDir params.mappedAllTrxOut, mode:'copy'



    input:
       path index_trx
       tuple val(sample_id), path(fastq_reads)

    output:
       tuple val(sample_id), path("${sample_id}_filt_alns.transcriptome.minimap2.bam"), emit: bam_trx
       tuple val(sample_id), path("${sample_id}_alns.minimap2.mapped.transcriptome.txt"), emit: stats_map

    script:
    """
    minimap2 -t ${params.threads} -ax splice -uf ${index_trx} ${fastq_reads}\
        | samtools view -q ${params.minimum_mapping_quality} -F 2304 -hbo ${sample_id}_filt_alns.transcriptome.minimap2.bam

    samtools stats -F 4 ${sample_id}_filt_alns.transcriptome.minimap2.bam > ${sample_id}_alns.minimap2.mapped.transcriptome.txt

    """
}

process salmon{
    /*
    Alignment based quantification using salmon
    */
    label "isoforms"
    cpus params.threads


    //added (AS 29v2023)
    publishDir params.salmonOut, mode:'copy'

    input:
       tuple val(sample_id), path(bam)
       path(ref_transcriptome_fa)

    output:
       path("${sample_id}_quant_salmon/*"), emit: quantout

    script:
    """
    salmon quant --noErrorModel -p ${params.threads} -t "${ref_transcriptome_fa}" -l SF -a "${bam}" -o "${sample_id}_quant_salmon"
    """
}




workflow map_reads_all_transcriptome {
    take:
       index_trx
       ref_transcriptome_fa
       fastq_reads
    main:
        map_reads_trx_all(index_trx, fastq_reads)
        salmon(map_reads_trx_all.out.bam_trx, ref_transcriptome_fa)
    emit:
       bam = map_reads_trx_all.out.bam_trx
       stats_map = map_reads_trx_all.out.stats_map
       salmon_quant=salmon.out.quantout
}


