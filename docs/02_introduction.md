This workflow analyses Oxford Nanopore long-read transcript sequencing data. It uses
[`bambu`](https://bioconductor.org/packages/bambu/) to build and quantify
transcript models,
[`SQANTI3`](https://github.com/ConesaLab/SQANTI3) for transcript classification
and QC, can optionally run
[`DESeq2`](https://bioconductor.org/packages/DESeq2/) and
[`DEXSeq`](https://bioconductor.org/packages/DEXSeq/) for differential
analysis, and will run
[`modkit`](https://github.com/nanoporetech/modkit) for base modification
pileups on aligned reads if relevant tags are present.

The workflow supports:

+ Transcript identification from either cDNA or direct RNA reads
+ Transcript discovery guided by a supplied genome and annotation
+ Quantification against a supplied reference annotation
+ Transcript classification and QC with `SQANTI3`
+ Differential gene expression with `DESeq2`
+ Differential transcript usage with `DEXSeq`
+ Base modification pileups with `modkit`

The main transcriptome result is a shared `bambu` model built from all samples
together. The workflow also produces separate per-sample transcriptomes, so
each sample has its own GTF, FASTA, count tables, and `SQANTI3`
summary alongside the shared results.

<figure>
<img src="docs/images/wf-transcriptomes.drawio.svg" alt="wf-transcriptomes overview schematic."/>
<figcaption>Schematic depicting wf-transcriptomes workflow.</figcaption>
</figure>

For users familiar with earlier transcriptome workflows, the main change is
that transcript discovery, quantification, and optional differential analysis
now use the shared `bambu` outputs rather than the older
StringTie/GffCompare/Salmon-based approach. The rest of this README explains
the current workflow in plain terms, while the `FAQ` and `Troubleshooting`
sections call out the main differences from the previous workflow version.
For `bambu` implementation details, see the
[`bambu` GitHub repository](https://github.com/GoekeLab/bambu). For
`SQANTI3` classification categories, see the
[`SQANTI3` isoform classification documentation](https://github.com/ConesaLab/SQANTI3/wiki/).
