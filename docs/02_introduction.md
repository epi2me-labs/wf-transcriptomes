This workflow analyses Oxford Nanopore long-read RNA sequencing data. It uses
[`bambu`](https://bioconductor.org/packages/bambu/) to build and quantify
transcript models, can optionally run
[`SQANTI3`](https://github.com/ConesaLab/SQANTI3) for transcript classification
and QC, and can optionally run
[`DESeq2`](https://bioconductor.org/packages/DESeq2/) and
[`DEXSeq`](https://bioconductor.org/packages/DEXSeq/) for differential
analysis.

The workflow supports:

+ transcript identification from either cDNA or direct RNA reads
+ transcript discovery guided by a supplied genome and annotation
+ quantification against a supplied reference annotation
+ optional transcript classification and QC with `SQANTI3`
+ differential gene expression with `DESeq2`
+ differential transcript usage with `DEXSeq`

The main transcriptome result is a shared `bambu` model built from all samples
together. The workflow also produces separate per-sample transcriptomes, so
each sample has its own GTF, FASTA, count tables, and optional `SQANTI3`
summary alongside the shared results.

For users familiar with earlier transcriptome workflows, the main change is
that transcript discovery, quantification, and optional differential analysis
now use the shared `bambu` outputs rather than the older
StringTie/GffCompare/Salmon-based approach. The rest of this README explains
the current workflow in plain terms, while the `FAQ` and `Troubleshooting`
sections call out the main differences from the previous workflow version.
