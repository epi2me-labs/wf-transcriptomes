
process jaffal{
    label "isoforms"
    cpus params.threads
    memory "31 GB"
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
    # Jaffa requires a virtual environment
    # so override JAVA_TOOL_OPTIONS env variable.
    JAVA_TOOL_OPTIONS=""
    JAFFAOUT=jaffal_output_$sample_id

    # JAFFAL exists with status code 1 when there's 0 fusion hits. Prevent this with '||:'
    $params.jaffal_dir/tools/bin/bpipe run \
        -n "${task.cpus}" \
        -p jaffa_output="\$JAFFAOUT/" \
        -p refBase=$refBase \
        -p genome=$genome \
        -p annotation=$annotation \
        -p fastqInputFormat="*.fastq" \
        $params.jaffal_dir/JAFFAL.groovy \
        $fastq || :

    summary="\$JAFFAOUT/all/all.summary"

    if [ -f \$summary ]; then
        # The summary is writtten so assume JAFFAL completed.
        if [ ! -s \$summary ]; then
            echo "JAFFAL failed to find any fusion transcripts for ${sample_id}"
            touch "\$JAFFAOUT/${sample_id}_jaffa_results.csv"
        else
            echo JAFFAL found fusion transcripts for ${sample_id}
            mv "\$JAFFAOUT/jaffa_results.csv" "\$JAFFAOUT/${sample_id}_jaffa_results.csv"
            # Add sample id column and header
            sed "s/\$/,${sample_id}/" \$JAFFAOUT/${sample_id}_jaffa_results.csv \
                | sed "1 s/${sample_id}/sample_id/" > tmp
            mv tmp \$JAFFAOUT/${sample_id}_jaffa_results.csv
        fi
    else
        echo JAFFAL encountered an error while prosessing ${sample_id}
    fi
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
