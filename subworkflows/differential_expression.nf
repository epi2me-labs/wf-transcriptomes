nextflow.enable.dsl = 2


process runDifferentialAnalysis {
    label "wf_transcriptomes"
    cpus 4
    memory "32 GB"
    input:
        path transcript_rds
        path gene_rds
        path sample_sheet
    output:
        path "de_analysis", emit: dir
    script:
        String covariates_arg = params.covariates ? "--covariates '${params.covariates}'" : ""
        String reference_arg = params.reference_level ? "--reference_level '${params.reference_level}'" : ""
    """
    supeRglue de_analysis \
        --transcript_rds "${transcript_rds}" \
        --gene_rds "${gene_rds}" \
        --sample_sheet "${sample_sheet}" \
        --condition_column "${params.condition_column}" \
        ${covariates_arg} \
        ${reference_arg} \
        --out_dir de_analysis
    """
}


workflow differential_expression {
    take:
        transcript_rds
        gene_rds
        sample_sheet
    main:
        results = runDifferentialAnalysis(transcript_rds, gene_rds, sample_sheet)
    emit:
        dir = results.dir
}
