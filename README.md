# Workflow Transcriptomes

Transcriptome analysis including assembly and annotation of cDNA and direct RNA sequencing data and differential expression analysis.



## Introduction

This workflow can be used for the following:

+ Identify RNA transcripts using either cDNA or direct RNA reads.
+ Reference aided transcriptome assembly.
+ Annotation of assembled transcripts.
+ Differential gene expression analysis using a pre-computed or assembled reference transcriptome.
+ Differential transcript usage analysis using a precomputed or assembled reference transcriptome.



## Compute requirements

Recommended requirements:

+ CPUs = 16
+ Memory = 32GB

Minimum requirements:

+ CPUs = 8
+ Memory = 32GB

Approximate run time: 15 minutes per sample, with 1 million reads and recommended resources.

ARM processor support: False




## Install and run


These are instructions to install and run the workflow on command line.
You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop)
or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
to provide isolation of the required software.
Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.
This is controlled by the
[`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
parameter as exemplified below.

It is not required to clone or download the git repository
in order to run the workflow.
More information on running EPI2ME workflows can
be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow.
This will pull the repository in to the assets folder of
Nextflow and provide a list of all parameters
available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-transcriptomes --help
```
To update a workflow to the latest version on the command line use
the following command:
```
nextflow pull epi2me-labs/wf-transcriptomes
```

A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-transcriptomes/wf-transcriptomes-demo.tar.gz
tar -xzvf wf-transcriptomes-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using:
```
nextflow run epi2me-labs/wf-transcriptomes \
	--de_analysis \
	--direct_rna \
	--fastq 'wf-transcriptomes-demo/differential_expression_fastq' \
	--minimap2_index_opts '-k15' \
	--ref_annotation 'wf-transcriptomes-demo/gencode.v22.annotation.chr20.gtf' \
	--ref_genome 'wf-transcriptomes-demo/hg38_chr20.fa' \
	--sample_sheet 'wf-transcriptomes-demo/sample_sheet.csv' \
	-profile standard
```

For further information about running a workflow on
the command line see https://labs.epi2me.io/wfquickstart/




## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts either FASTQ or BAM files as input.

The FASTQ or BAM input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ or BAM file; (ii) the path to a top-level directory containing FASTQ or BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ or BAM files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

```
(i)                     (ii)                 (iii)    
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```



## Input parameters

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| bam | string | BAM or unaligned BAM (uBAM) files to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| transcriptome_source | string | Select how the transcriptome used for analysis should be prepared. | For differential expression analysis, use of an existing transcriptome may be preferred and so 'precomputed' should be selected. In this case the 'ref_transcriptome' parameter should be specified. To create a reference transcriptome using an existing reference genome, select 'reference guided' and specify the 'ref_genome' parameter. | reference-guided |
| ref_genome | string | Path to reference genome sequence [.fa/.fq/.fa.gz/fq.gz]. Required for reference-based workflow. | A reference genome is required for reference-based assembly of a transcriptome. |  |
| ref_transcriptome | string | Transcriptome reference file. Required for precomputed transcriptome calculation and for differential expression analysis. | A reference transcriptome related to the sample under study. Must be supplied when the 'Transcriptome source' parameter has been set to 'precomputed' or to perform differential expression. |  |
| ref_annotation | string | A reference annotation in GFF2 or GFF3 format (extensions .gtf(.gz), .gff(.gz), .gff3(.gz)). Only annotation files from [Encode](https://www.encodeproject.org), [Ensembl](https://www.ensembl.org/index.html) and [NCBI](https://www.ncbi.nlm.nih.gov/) are supported. | This will be used for guiding the transcriptome assembly and to label transcripts with their corresponding gene identifiers. Note: If in de_analysis mode transcript strands must be only + or -. |  |
| direct_rna | boolean | Set to true for direct RNA sequencing. |  Omits the pychopper step. | False |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all user-facing files. |  | output |
| igv | boolean | Visualize outputs in the EPI2ME IGV visualizer. | Enabling this option will visualize the output alignment files in the EPI2ME Desktop App IGV visualizer. | False |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. If you are running the differential expression workflow, there must be an additional column `condition` with two labels, one of which must be `control` (e.g. `control` and `treated`). Control will indicate which samples will be used as the reference. There should be at least 3 repeats for each condition. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Options for reference-based workflow

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| plot_gffcmp_stats | boolean | Create a PDF of plots from showing gffcompare results | If set to true, a PDF file containing detailed gffcompare reults will be output | True |
| gffcompare_opts | string | Extra command-line options to give to gffcompare -r | For a list of possible options see [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml). | -R |
| minimap2_index_opts | string | Extra command-line options for minimap2 indexing. | See [minimap2 index options](https://lh3.github.io/minimap2/minimap2.html#4) for more information. These will only be relevant in the reference based transcriptome assembly. | -k14 |
| minimap2_opts | string | Additional command-line options for minimap2 alignment. | See [minimap2 options](https://lh3.github.io/minimap2/minimap2.html#5) for further information. These will only be relevant in the reference based transcriptome assembly. | -uf |
| minimum_mapping_quality | integer | filter aligned reads by MAPQ quality. | Reads that do not meet this mapping quality after minimap2 alignment, will be filtered out. | 40 |
| stringtie_opts | string | Extra command-line options for stringtie transcript assembly. | For additional String tie options see [here](https://github.com/gpertea/stringtie#stringtie-options). | --conservative |


### Differential Expression Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| de_analysis | boolean | Run DE anaylsis | Running this requires you to provide at least two replicates for a control and treated sample as well as a sample sheet param. | False |
| min_gene_expr | integer | The minimum number of total mapped sequence reads required for a gene to be considered in differential transcript usage analysis. | Filtering at the gene level ensures that the observed transcript ratios are calculated with a minimum number of counts per gene. | 10 |
| min_feature_expr | integer | The minimum number of reads assigned to a transcript for it to be considered in differential transcript usage analysis. | Filter out transcripts that do not have this minimum number of transcript expression, reducing noise. | 3 |
| min_samps_gene_expr | integer | Set the minimum number of samples in which a gene is expressed to be included in the differential transcript usage analysis. | A gene must be expressed in at least this number of samples for the gene be included in the differential transcript usage analysis. Filtering at the gene level improves the reliability of the observed transcript ratios. | 3 |
| min_samps_feature_expr | integer | Set the minimum number of samples in which a transcript is expressed to be included in the differential transcript usage analysis. | A transcript must expressed in at least this minimum number of samples to be included in the analysis. Should be equal to the number of replicates per sample you have. | 1 |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Number of CPU threads. | Only provided to processes including alignment and and assembly that benefit from multiple threads. | 4 |
| cdna_kit | string | If cDNA reads are used, select the kit used. | This will be used by pychopper to preprocess the reads for downstream analysis. | SQK-PCS109 |
| pychopper_backend | string | Pychopper can use one of two available backends for identifying primers in the raw reads | 'edlib' is set by default due to its high performance. However, it may be less sensitive than 'phmm'. | edlib |
| pychopper_opts | string | Extra pychopper opts | See available options (here)[https://github.com/epi2me-labs/pychopper#usage] |  |
| bundle_min_reads | integer | Minimum size of bam bundle for parallel processing. |  | 50000 |
| isoform_table_nrows | integer | Maximum rows to dispay in the isoform report table |  | 5000 |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | wf-transcriptomes-report.html | a HTML report document detailing the primary findings of the workflow | aggregated |
| Per file read stats | fastq_ingress_results/{{ alias }}//reads/fastcat_stats/per-file-stats.tsv | A TSV with per file read stats, including all samples. | aggregated |
| Read stats | fastq_ingress_results/{{ alias }}//reads/fastcat_stats/per-read-stats.tsv | A TSV with per read stats, including all samples. | aggregated |
| Run ID's | fastq_ingress_results/{{ alias }}//reads/fastcat_stats/run_ids | List of run IDs present in reads. | aggregated |
| Meta map json | fastq_ingress_results/{{ alias }}//reads/metamap.json | Metadata used in workflow presented in a JSON. | aggregated |
| Concatenated sequence data | fastq_ingress_results/{{ alias }}//reads/{{ alias }}.fastq.gz | Per sample reads concatenated in to one FASTQ file. | per-sample |
| Assembled transcriptome | {{ alias }}_transcriptome.fas | Per sample assembled transcriptome.  Not output if a reference annotation was supplied | per-sample |
| Annotated assembled transcriptome | {{ alias }}_merged_transcriptome.fas | Per sample annotated assembled transcriptome. Only output if a reference annotation was supplied | per-sample |
| Alignment summary statistics | {{ alias }}_read_aln_stats.tsv | Per sample alignment summary statistics. | per-sample |
| GFF compare results. | {{ alias }}_gffcompare | All GFF compare output files. | per-sample |
| Differential gene expression results | de_analysis/results_dge.tsv | This is a gene-level result file that describes genes and their probability of showing differential expression between experimental conditions. | aggregated |
| Differential gene expression report | de_analysis/results_dge.pdf | Summary report of differential gene expression analysis as a PDF. | aggregated |
| Differential transcript usage gene TSV | de_analysis/results_dtu_gene.tsv | This is a gene-level result file from DEXSeq that lists annotated genes and their probabilities of differential expression. | aggregated |
| Differential transcript usage report | de_analysis/results_dtu.pdf | Summary report of differential transcript usage results as a PDF. | aggregated |
| Differential transcript usage TSV | de_analysis/results_dtu_transcript.tsv | This is a transcript-level result file from DEXSeq that lists annotated genes and their probabilities of differential expression. | aggregated |
| Differential transcript usage stageR TSV | de_analysis/results_dtu_stageR.tsv  | This is the output from StageR and it shows both gene and transcript probabilities of differential expression | aggregated |
| Differential transcript usage DEXSeq TSV | de_analysis/results_dexseq.tsv | The complete output from the DEXSeq-analysis, shows both gene and transcript probabilities of differential expression. | aggregated |
| Gene counts | de_analysis/all_gene_counts.tsv | Raw gene counts created by the Salmon tool, before filtering. | aggregated |
| Gene counts per million | de_analysis/cpm_gene_counts.tsv | This file shows counts per million (CPM) of the raw gene counts to facilitate comparisons across samples. | aggregated |
| Transcript counts | de_analysis/unfiltered_transcript_counts_with_genes.tsv | Raw transcript counts created by the Salmon tool, before filtering. Includes reference to the associated gene ID. | aggregated |
| Transcript per million counts | de_analysis/unfiltered_tpm_transcript_counts.tsv | This file shows transcripts per million (TPM) of the raw counts to facilitate comparisons across samples. | aggregated |
| Transcript counts filtered | de_analysis/filtered_transcript_counts_with_genes.tsv | Filtered transcript counts, used for differential transcript usage analysis. Includes a reference to the associated gene ID. | aggregated |
| Transcript info table | {{ alias }}_transcripts_table.tsv | This file details each isoform that was reconstructed from the input reads. It contains a subset of columns from the .tmap output from [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) | per-sample |
| Final non redundant transcriptome | de_analysis/final_non_redundant_transcriptome.fasta | Transcripts that were used for differential expression analysis including novel transcripts with the identifiers used for DE analysis. Only applicable when the ref_transcriptome parameter is not provided. | aggregated |
| Index of reference FASTA file | igv_reference/{{ ref_genome file }}.fai | Reference genome index of the FASTA file required for IGV config. | aggregated |
| GZI index of the reference FASTA file | igv_reference/{{ ref_genome file }}.gzi | GZI Index of the reference FASTA file. | aggregated |
| JSON configuration file for IGV browser | igv.json | JSON configuration file to be loaded in IGV for visualising alignments against the reference. | aggregated |
| BAM file (minimap2) | BAMS/{{ alias }}.reads_aln_sorted.bam | BAM file generated from mapping input reads to the reference. | per-sample |
| BAM index file (minimap2) | BAMS/{{ alias }}.reads_aln_sort.bam.bai | Index file generated from mapping input reads to the reference. | per-sample |




## Pipeline overview

### 1. Concatenate input files and generate per read stats.
The [fastcat](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow. It will also output per read stats including average read lengths and qualities.

### 2. Preprocess cDNA.
If input sequences are cDNA [Pychopper](https://github.com/epi2me-labs/pychopper) is used to orient, trim and rescue full length cDNA reads and associated statistics. If the `direct_rna` parameter is selected this step will be skipped.

### 3. Build transcriptome.
If the `transcriptome_source` parameter is "reference-guided" a transcriptome will be built for each sample as outlined below. If the `transcriptome_source` is "precomputed" and the `reference_transcriptome` parameter is provided the workflow will skip step 3.

#### 3.1 Align reads with reference genome.
The reference genome will be indexed and aligned using [Minimap2](https://github.com/lh3/minimap2). The output is sorted and converted to a BAM file using [Samtools](https://www.htslib.org/). Alignment stats are created from these using [Seqkit BAM](https://bioinf.shenwei.me/seqkit/usage/#bam).

Additionally, the workflow will generate an IGV configuration file if `--igv` is selected. This file allows the user to view the aligned BAM in the EPI2ME Desktop Application in the Viewer tab.

#### 3.2 Chunk BAM
The aligned BAMs are split into chunks using the bundle_min_reads parameter (default: 50000).

#### 3.3 Assemble transcripts
[StringTie](https://ccb.jhu.edu/software/stringtie/) is then used to assemble the transcripts using the aligned segments in the chunked BAM files. The assembled transcript will be output as a [GFF file](https://www.ensembl.org/info/website/upload/gff3.html). If a `ref_annotation` file is provided this will also be included in the GFF.

#### 3.4 Merge Chunks
Transcript GFF files from the chunks with the same sample aliases will then be merged.

#### 3.5 Annotate
[GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.html) is then used to compare query and reference annotations, merging records where appropriate and then annotating them. This also creates estimates of accuracy of the GFF files output in a stats file per sample.

#### 3.6 Create transcriptomes
[Gffread](https://github.com/gpertea/gffread) is used to create a transcriptome FASTA file from the final GFF as well as a merged transcriptome that includes annotations in the FASTA headers where available.

### 4. Differential expression analysis

Differential gene expression (DGE) and differential transcript usage (DTU) analyses aim to identify genes and transcripts that show statistically altered expression patterns.

Differential Expression requires at least 2 replicates of each sample to compare (but we recommend three). You can see an example sample_sheet.csv below.

#### Sample sheet condition column
The sample sheet should be a comma separated values file (.csv) and include at least three columns named `barcode`, `alias` and `condition`.
- Each `barcode` should refer to a directory of the same name in the input FASTQ directory (in the example below `barcode01` to `barcode06` reflect the `test_data` directory).
- The `alias` column allows you to rename each barcode to an alias that will be used in the report and other output files.
- The condition column will need to contain one of two keys to indicate the two samples being compared. Control must be one of the keys, used to indicate which samples will be used as the reference in the differential expression analysis.

eg. sample_sheet.csv
```
barcode,alias,condition
barcode01,sample01,control
barcode02,sample02,control
barcode03,sample03,control
barcode04,sample04,treated
barcode05,sample05,treated
barcode06,sample06,treated
```

#### 4.1 Merge cross sample transcriptomes
If a `ref_transcriptome` is not provided, the transcriptomes created by the workflow will be used for DE analysis. To do this, the GFF outputs of GffCompare are merged using StringTie. A final non redundant FASTA file of the transcripts is created using the merged GFF file and the reference genome using seqkit.

#### 4.2 Create a final non redundant transcriptome
The reads from all the samples will be aligned with the final non redundant transcriptome using Minimap2 in a splice aware manner.

#### 4.3 Count genes and transcripts
[Salmon](https://github.com/COMBINE-lab/salmon) is used for transcript quantification, giving gene and transcript counts.

#### 4.4 edgeR based differential expression analysis
A statistical analysis is first performed using [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) to identify the subset of differentially expressed genes using the gene counts as input. A normalisation factor is calculated for each sequence library using the default TMM method (see [McCarthy et al. (2012)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3378882/) for further details). The defined experimental design is used to calculate estimates of dispersion for each of the gene features. Statistical tests are calculated using the contrasts defined in the experimental design. The differentially expressed genes are corrected for false discovery (FDR) using the method of Benjamini & Hochberg ([Benjamini and Hochberg (1995)](https://www.jstor.org/stable/2346101))

#### 4.5 Pre-filtering of quantitative data using DRIMSeq
[DRIMSeq](https://bioconductor.org/packages/release/bioc/html/DRIMSeq.html) is used to filter the transcript count data from the Salmon analysis for differential transcript usage (DTU) analysis. The filter step will be used to select for genes and transcripts that satisfy rules for the number of samples in which a gene or transcript must be observed, and minimum threshold levels for the number of observed reads. The parameters used for filtering are `min_samps_gene_expr`, `min_samps_feature_expr`, `min_gene_expr`, and `min_feature_expr`. By default, any transcripts with zero expression or one transcript in all samples are filtered out at this stage.

#### 4.6 Differential transcript usage using DEXSeq
Differential transcript usage analysis is performed using the R [DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html) package ([Anders et al. (2012)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3460195/)). Similar to the edgeR package, DEXSeq estimates the variance between the biological replicates and applies generalised linear models for the statistical testing. The key difference is that the DEXSeq method looks for differences at the exon count level. DEXSeq uses the filtered transcript count data prepared earlier in this analysis. 

#### 4.7 StageR stage-wise analysis of DGE and DTU
The final component of this isoform analysis is a stage-wise statistical test using the R software package [stageR](https://bioconductor.org/packages/release/bioc/html/stageR.html)([Van den Berge and Clement (2018)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1277-0)). stageR uses (1) the raw p-values for DTU from the DEXSeq analysis in the previous section and (2) a false-discovery corrected set of p-values from testing whether individual genes contain at least one exon showing DTU. A hierarchical two-stage statistical testing evaluates the set of genes for DTU.




## Troubleshooting

+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).
+ Renaming, moving or deleting the input BAM, reference genome or the output directory from the location provided at runtime will stop IGV in the EPI2ME Desktop app from loading.



## FAQ's

*Does the workflow support de novo assembly?* - Currently the workflow does not have a *de novo* mode.

*Why is the IGV panel not showing?* - The workflow expects either an uncompressed or [`bgzip`](https://www.htslib.org/doc/bgzip.html)-compressed reference. If the user provides a reference compressed not with `bgzip`, the workflow will run to completion, but won't be able to generate the necessary indexes to visualize the outputs in IGV.

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-transcriptomes/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).



## Related blog posts

+ [How to align your data](https://labs.epi2me.io/how-to-align/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



