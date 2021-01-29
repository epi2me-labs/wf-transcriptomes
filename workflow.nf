#!/usr/bin/env extflow

// Developer notes
// 
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion.
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

nextflow.enable.dsl = 2

params.help = ""

if(params.help) {
    log.info ''
    log.info 'Workflow template'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run workflow.nf [options]'
    log.info ''
    log.info 'Script Options: '
    log.info '    --fastq        FILE    Path to FASTQ file'
    log.info '    --out_dir      DIR     Path for output'
    log.info ''

    return
}


process readSeqs {
    // Just write a file with sequence lengths
    label "containerCPU"
    input:
        file reads
    output:
        file "seqs.txt"

    """
    #!/usr/bin/env python
    import pysam
    with open("seqs.txt", 'w') as fh:
        for rec in pysam.FastxFile("$reads"):
            fh.write("{}\\t{}\\n".format(rec.name, len(rec.sequence)))
    """
}



// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}


// workflow module
workflow pipeline {
    take:
        reads
    main:
        seqs = readSeqs(reads)
    emit:
        seqs
}

// entrypoint workflow
workflow {
    reads = channel.fromPath(params.reads, checkIfExists:true)
    results = pipeline(reads)
    output(results)
}
