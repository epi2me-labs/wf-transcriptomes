#!/usr/bin/env nextflow

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


def helpMessage(){
    log.info """
Workflow template'

Usage:
    nextflow run epi2melabs/wf-template [options]

Script Options:
    --fastq        DIR     Path to directory containing FASTQ files (required)
    --out_dir      DIR     Path for output (default: $params.out_dir)
"""
}


process concatFastq {
    // concatenate fastq and fastq.gz in a dir

    label "pysam"
    cpus 1
    input:
        file "input"
    output:
        file "reads.fastq.gz"

    shell:
    '''
#!/usr/bin/env python
from glob import glob
import gzip
import itertools
import os
import pysam

# we use pysam just because it will read both fastq and fastq.gz
# and we don't have to worry about having a combination or not
with gzip.open("reads.fastq.gz", "wt") as fh:
    files = itertools.chain(
        glob("input/*.fastq"), glob("input/*.fastq.gz"))
    records = itertools.chain.from_iterable(
        pysam.FastxFile(fn) for fn in files) 
    for rec in records:
        annot = " {}".format(rec.comment) if rec.comment else ""
        qual = rec.quality if rec.quality else "+"*len(rec.sequence)
        fh.write("@{}{}\\n{}\\n+\\n{}\\n".format(rec.name, annot, rec.sequence, qual))
    '''
}


process readSeqs {
    // Just write a file with sequence lengths
    label "pysam"
    input:
        file reads
    output:
        file "seqs.txt"

    """
    read_lengths.py $reads seqs.txt
    """
}


process makeReport {
    label "pysam"
    input:
        file "seqs.txt"
    output:
        file "report.html"
    """
    report.py report.html seqs.txt
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
        reads = concatFastq(reads)
        summary = readSeqs(reads)
        report = makeReport(summary)
    emit:
        summary.concat(report)
}

// entrypoint workflow
workflow {

    if (params.help) {
        helpMessage()
        exit 1
    }

    if (!params.fastq) {
        helpMessage()
        println("")
        println("`--fastq` is required")
        exit 1
    }


    reads = file("$params.fastq/*.fastq*", type: 'file', maxdepth: 1)
    if (reads) {
        reads = Channel.fromPath(params.fastq, type: 'dir', checkIfExists: true)
        results = pipeline(reads)
        output(results)
    } else {
        println("No .fastq(.gz) files found under `${params.fastq}`.")
    }
}
