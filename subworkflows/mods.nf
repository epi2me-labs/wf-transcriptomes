// Check that the bam has modifications
process validate_modbam {
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        tuple val(meta),
            path(alignment),
            path(alignment_index),
            val(alignment_stats)
    output:
        tuple val(meta),
            path(alignment),
            path(alignment_index),
            val(alignment_stats),
            env(valid)

    script:
    """
    valid=0
    workflow-glue check_valid_modbam ${alignment} || valid=\$?

    # Allow EX_OK and EX_DATAERR, otherwise explode
    if [ \$valid -ne 0 ] && [ \$valid -ne 65 ]; then
        exit 1
    fi
    """
}

process runModkitPileup {
    label "modkit"
    cpus { params.threads ?: 4 }
    memory "16 GB"
    input:
        tuple val(alias),
            val(meta),
            path(alignment),
            path(alignment_index),
            val(alignment_stats),
            val(mod_codes)
        tuple path(reference),
            path(reference_index)
    output:
        tuple val(alias),
            path("${alias}.mods.bedmethyl.gz"),
            emit: bedmethyl
    publishDir "${params.out_dir}/${output_key}/mods", mode: 'copy'

    script:
    output_key = alias == "cohort" ? "cohort" : "samples/${alias}"  // nodef
    String modified_bases_arg = "--modified-bases " + mod_codes
        .split(',')
        .join(' ')
    """
    modkit pileup \
        "${alignment}" \
        "${alias}.mods.bedmethyl.gz" \
        ${modified_bases_arg} \
        --reference "${reference}" \
        --threads ${task.cpus} \
        --bgzf
    """
}

process modkit_tobigwig {
    label "modkit"
    cpus 4
    memory "2 GB"
    input:
        tuple path(reference),
            path(reference_index)
        tuple val(alias),
            path(bedmethyl),
            val(mod_codes)
        path mod_code_labels
    output:
        tuple val(alias),
            path("${alias}.mods.*.bw"),
            emit: bigwig
    publishDir "${params.out_dir}/${output_key}/mods", mode: 'copy'

    script:
    output_key = alias == "cohort" ? "cohort" : "samples/${alias}"  // nodef
    String mod_code_args = mod_codes
        .split(',')
        .join(' ')
    """
    for mod_code in ${mod_code_args}; do
        mod_code_value="\${mod_code#*:}"
        mod_label=\$(mod_code_label "\${mod_code}" "${mod_code_labels}")
        zcat "${bedmethyl}" | \
            modkit bedmethyl tobigwig \
                --sizes "${reference_index}" \
                --nthreads ${task.cpus} \
                --mod-codes "\${mod_code_value}" \
                - \
                "${alias}.mods.\${mod_label}.bw"
    done
    """
}

process summariseModkitBedmethyl {
    label "wf_common"
    cpus 1
    memory "2 GB"
    input:
        tuple val(alias),
            path(bedmethyl),
            val(mod_codes)
        path mod_code_labels
    output:
        tuple val(alias),
            path("${alias}.mods.summary.tsv"),
            emit: summary
    publishDir "${params.out_dir}/${output_key}/mods", mode: 'copy'

    script:
    output_key = alias == "cohort" ? "cohort" : "samples/${alias}"  // nodef
    """
    workflow-glue summarise_modkit_bedmethyl \
        "${bedmethyl}" \
        "${alias}" \
        "${mod_codes}" \
        "${mod_code_labels}" \
        "${alias}.mods.summary.tsv"
    """
}

process inferModkitBases {
    label "modkit"
    cpus 1
    memory "4 GB"
    input:
        tuple val(alias),
            val(meta),
            path(alignment),
            path(alignment_index),
            val(alignment_stats)
    output:
        tuple val(alias),
            env(mod_codes)

    script:
    """
    modkit modbam check-tags "${alignment}" --num-reads 10000 --mapped-only --out-dir check_tags
    infer_modkit_codes check_tags/modified_bases.tsv > mod_codes.txt
    if [ ! -s mod_codes.txt ]; then
        echo "Failed to infer modified base codes from ${alignment}" >&2
        exit 1
    fi

    read -r mod_codes < mod_codes.txt
    """
}

workflow mods {
    take:
        xams
        ref_genome
    main:
        // Check inputs have modtags, we'll early abort mod analysis for any that don't
        validate_modbam(xams)
            .branch {
                nomods: it[-1] == '65'
                    return it[0].alias
                mods: it[-1] == '0'
                    return [it[0].alias] + it[0..-2]  // prepend alias for joining and drop exit_code marker
            }
            .set{xams_with}

        // warn for samples without mods
        xams_with.nomods.subscribe {
            log.warn "Input ${it} does not contain modified base tags. Was a modified basecalling model selected when basecalling this data?"
        }

        // determine what mods to ask modkit pileup for
        sample_modcodes = params.mod_codes
            ? xams_with.mods.map { [it[0], params.mod_codes.trim()] }  // cross all aliases with user mod_codes
            : inferModkitBases(xams_with.mods)                         // otherwise infer per-sample from modbam
        mod_samples = xams_with.mods.join(sample_modcodes)

        pileup = runModkitPileup(mod_samples, ref_genome)
        sample_summaries = summariseModkitBedmethyl(
            pileup.bedmethyl.join(sample_modcodes),
            file("$projectDir/data/mod_code_labels.tsv")
        )

        sample_bigwigs = modkit_tobigwig(
            ref_genome,
            pileup.bedmethyl.join(sample_modcodes),
            file("$projectDir/data/mod_code_labels.tsv")
        )

    emit:
        bedmethyl = pileup.bedmethyl
        summary = sample_summaries.summary
        bigwig = sample_bigwigs.bigwig
}
