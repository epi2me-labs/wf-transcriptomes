
process jaffal{
    label "isoforms"
    input:
        tuple val(sample_id), path(fastq)
        path refBase
        val genome
        val annotation
    output:
        tuple val(sample_id), path("jaffal_output_$sample_id"), emit: results
        tuple val(sample_id), path("jaffal_output_$sample_id/*jaffa_results.csv"), emit: results_csv
    script:
    """
    JAFFAOUT=jaffal_output_$sample_id
    $params.jaffal_dir/tools/bin/bpipe run \
        -n $params.threads \
        -p jaffa_output="\$JAFFAOUT/" \
        -p refBase=$refBase \
        -p genome=$genome \
        -p annotation=$annotation \
        -p fastqInputFormat="*.fastq" \
        $params.jaffal_dir/JAFFAL.groovy \
        $fastq
    mv "\$JAFFAOUT/jaffa_results.csv" "\$JAFFAOUT/${sample_id}_jaffa_results.csv"

    # Add sample id column
    sed "s/\$/,${sample_id}/" \$JAFFAOUT/${sample_id}_jaffa_results.csv > tmp1
    # Add header
    sed "1 s/${sample_id}/sample_id/" tmp1 > tmp2
    mv tmp2 \$JAFFAOUT/${sample_id}_jaffa_results.csv
    """
}

// workflow module
workflow gene_fusions {
    take:
        fastq
        refBase
        genome
        annotation
    main:
        jaffal(fastq, refBase, genome, annotation)
    emit:
        results_csv = jaffal.out.results_csv
        results = jaffal.out.results
}
