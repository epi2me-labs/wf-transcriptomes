### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| bam | string | BAM or unaligned BAM (uBAM) files to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| transcriptome_source | string | Select how the transcriptome used for analysis should be prepared. | To analyse only gene fusions and differential expression use of an existing transcriptome may be preferred and so 'precomputed' should be selected. In this case the 'ref_transcriptome' parameter should be specified. To create a reference transcriptome using an existing reference genome, select 'reference guided' and specify the 'ref_genome' parameter. | reference-guided |
| ref_genome | string | Path to reference genome sequence [.fa/.fq/.fa.gz/fq.gz]. Required for reference-based workflow. | A reference genome is required for reference-based assembly of a transcriptome. |  |
| ref_transcriptome | string | Transcriptome reference file. Required for precomputed transcriptome calculation and for differential expression analysis. | A reference transcriptome related to the sample under study. Must be supplied when the 'Transcriptome source' parameter has been set to 'precomputed' or to perform differential expression. |  |
| ref_annotation | string | A reference annotation in GFF2 or GFF3 format (extensions .gtf(.gz), .gff(.gz), .gff3(.gz)). Only annotation files from [Encode](https://www.encodeproject.org), [Ensembl](https://www.ensembl.org/index.html) and [NCBI](https://www.ncbi.nlm.nih.gov/) are supported. | This will be used for guiding the transcriptome assembly and to label transcripts with their corresponding gene identifiers. |  |
| direct_rna | boolean | Set to true for direct RNA sequencing. |  Omits the pychopper step. | False |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all user-facing files. |  | output |


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


### Gene Fusion Detection Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| jaffal_refBase | string | JAFFAl reference genome directory. | JAFFAL human hg38 reference data directory can be downloaded from here: https://figshare.com/ndownloader/files/25410494 or see the README for alternative instructions. If custom gemome files are required, see the instructions here: https://github.com/Oshlack/JAFFA/wiki/FAQandTroubleshooting#how-can-i-generate-the-reference-files-for-a-non-supported-genome. |  |
| jaffal_genome | string | Genome reference prefix. e.g. hg38. | JAFFAL reference files are prefixed with the genome reference file name and need to be supplied . If using the human reference data provided by JAFFAL, this can be left at `hg38`. | hg38 |
| jaffal_annotation | string | Annotation suffix. | JAFFAL reference files are suffixed with the annotation filename and this needs to be supplied. For the human hg38 reference data supplied by JAFFAL, this is `genCode22`. | genCode22 |


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


