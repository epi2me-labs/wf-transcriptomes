## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop),
[Singularity](https://sylabs.io/singularity/) to provide isolation of
the required software. Each method is automated out-of-the-box provided
either docker or singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).


### Workflow options

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-transcriptomes --help
```

to see the options for the workflow.

**Download demonstration data**

A small test dataset is provided for the purposes of testing the workflow software. It consists of reads, reference,
and annotations from human chromosome 20 only.
It can be downloaded using:
```shell
wget -O test_data.tar.gz https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-isoforms/wf-isoforms_test_data.tar.gz 
tar -xzvf  test_data.tar.gz
```

**Example execution of a workflow for reference-based transcript assembly and fusion detection**
```
OUTPUT=~/output;
nexflow run epi2me-labs/wf-transcriptomes \
  --fastq ERR6053095_chr20.fastq \
  --ref_genome chr20/hg38_chr20.fa \
  --ref_annotation chr20/gencode.v22.annotation.chr20.gtf \
  --jaffal_refBase chr20/ \
  --jaffal_genome hg38_chr20 \
  --jaffal_annotation "genCode22" \
  --out_dir outdir -w workspace_dir
```

**Example workflow for denovo transcript assembly**
```
OUTPUT=~/output
nextflow run . --fastq test_data/fastq \
  --denovo \
  --ref_genome test_data/SIRV_150601a.fasta \
  --out_dir ${OUTPUT} \
  -w ${OUTPUT}/workspace \
  --sample sample_id
```
A full list of options can be seen in nextflow_schema.json. 
Parameters can be specified either in a config like `parameter = value` or on the command line like `--parameter value`.
Below are some commonly used parameters in the format used in config files.

Select how the transcriptome used for analysis should be prepared:

- To create a reference transcriptome using an existing reference genome `--transcriptome_source reference-guided` (default)
- Use a a supplied transcriptome `--transcriptome_source precomputed"`
- Gnerate transcriptome via the denovo pipeline `--transcriptome_source denovo"` 


To run the workflow with direct RNA reads `--direct_rna false` (this just skips the pychopper step).

Pychopper and minimap2 can take options via `--minimap2_opts` and `--pychopper_opts`, for example:

- When using the SIRV synthetic test data  
  - `--minimap2_opts '-uf --splice-flank=no'`
- pychopper needs to know which cDNA synthesis kit used, which can be specified with
  - SQK-PCS109: `--pychopper_opts '-k PCS109'` (default)
  - SQK-PCS110: `--pychopper_opts '-k PCS110'`
  - SQK-PCS111: `--pychopper_opts '-k PCS111'`
- pychopper can use one of two available backends for identifying primers in the raw reads
  - nhmmscan `--pychopper opts '-m phmm'` 
  - edlib `--pychopper opts '-m edlib'`

__Note__: edlib is set by default in the config as it's quite a lot faster. However, it may be less sensitive than nhmmscan. 

### Fusion detection

JAFFAL from the [JAFFA](https://github.com/Oshlack/JAFFA)
package is used to identify potential fusion transcripts.  

In order to use JAFFAL, reference files must first be downloaded.
To use pre-processed hg38 genome and GENCODE v22 annotation files (as used in the JAFFAL paper)
do:
```shell
mkdir jaffal_data_dir
cd jaffal_data_dir/
sh path/to/wf-transcriptomes/subworkflows/JAFFAL/download_jaffal_references.sh
````
Then the path to the directory containing the downloaded reference data must be specified with 
`--jaffal_refBase`.


**Using alternative genome and annotation files**

These should be prepared as described
[here](https://github.com/Oshlack/JAFFA/wiki/FAQandTroubleshooting#how-can-i-generate-the-reference-files-for-a-non-supported-genome).

The resulting JAFFAL reference files will look something like `hg38_genCode22.fa`. The following options enable JAFFAL to find these
files:

`--jaffal_genome reference_genome_name` optional (default: `hg38`)  
`--jaffal_annotation jaffal_annotation_prefix` optional (default: `genCode22`) 


__Note__: JAFFAL is not currently working on Mac M1 (osx-arm64 architecture).

### Differential Expression

Differential Expression requires at least 2 replicates of each sample to compare (but we recommend three). You can see an example sample_sheet.csv below.

**Example workflow for differential expression transcript assembly**

#### Sample sheet condition column
The sample sheet should be a comma separated values file (.csv) and include at least three columns named `barcode`, `alias` and `condition`.
- Each `barcode` should refer to a directory of the same name in the input FASTQ directory (in the example below `barcode01` to `barcode06` reflect the `test_data` directory).
- The `alias` column allows you to rename each barcode to an alias that will be used in the report and other output files.
- The condition column will need to contain one of two keys to indicate the two samples being compared. for example: treated/untreated, sample/control etc.

In the default `sample_sheet.csv` available in the test_data directory we have used the following.

eg. sample_sheet.csv
```
barcode,alias,condition
barcode01,sample01,untreated
barcode02,sample02,untreated
barcode03,sample03,untreated
barcode04,sample04,treated
barcode05,sample05,treated
barcode06,sample06,treated
```

You will also need to provide a reference genome and a reference annotation file.
Here is an example cmd to run the workflow. First you will need to download the data with wget. 
eg.
```
wget -O differential_expression.tar.gz https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-isoforms/differential_expression.tar.gz && tar -xzvf differential_expression.tar.gz
OUTPUT=~/output;
nextflow run epi2me-labs/wf-transcriptomes \
  --fastq  differential_expression/differential_expression_fastq \
  --de_analysis \
  --ref_genome differential_expression/hg38_chr20.fa \
  --ref_annotation differential_expression/gencode.v22.annotation.chr20.gtf \
  --direct_rna --minimap_index_opts \
  -k15
```
You can also run the differential expression section of the workflow on its own by providing a reference transcriptome and setting the transcriptome assembly parameter to false.
eg.
```
nextflow run epi2me-labs/wf-transcriptomes \
  --fastq  differential_expression/differential_expression_fastq \
  --de_analysis \
  --ref_genome differential_expression/hg38_chr20.fa \
  --ref_annotation differential_expression/gencode.v22.annotation.chr20.gtf \
  --direct_rna --minimap_index_opts \
  -k15 \
  --ref_transcriptome differential_expression/ref_transcriptome.fasta
```

## Workflow outputs
* an HTML report document detailing the primary findings of the workflow.
* for each sample:
  * [gffcomapre](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) output directories
  * read_aln_stats.tsv - alignment summary statistics
  * transcriptome.fas - the assembled transcriptome
  * merged_transcritptome.fas - annotated, assembled transcriptome
  * [jaffal](https://github.com/Oshlack/JAFFA) ooutput directories
  
### Fusion detection outputs
in `${out_dir}/jaffal_output_${sample_id}` you will find:
* jaffa_results.csv - the csv results summary file 
* jaffa_results.fasta - fusion transcritpt sequences

### Differential Expression outputs
* `de_analysis/results_dge.tsv` and `de_analysis/results_dge.pdf`- results of `edgeR` differential gene expression analysis.
* `de_analysis/results_dtu_gene.tsv`, `de_analysis/results_dtu_transcript.tsv` and `de_analysis/results_dtu.pdf` - results of differential transcript usage by `DEXSeq`.
* `de_analysis/results_dtu_stageR.tsv` - results of the `stageR` analysis of the `DEXSeq` output.
* `de_analysis/dtu_plots.pdf` - DTU results plot based on the `stageR` results and filtered counts.

### References

* Benjamini, Yoav, and Yosef Hochberg. 1995. “Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.” Journal of the Royal Statistical Society. Series B (Methodological) 57 (1): 289–300. http://www.jstor.org/stable/2346101.
* McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. “Differential Expression Analysis of Multifactor Rna-Seq Experiments with Respect to Biological Variation.” Nucleic Acids Research 40 (10): 4288–97.
* Nowicka, Malgorzata, and Mark D. Robinson. 2016. “DRIMSeq: A Dirichlet-Multinomial Framework for Multivariate Count Outcomes in Genomics [Version 2; Referees: 2 Approved].” F1000Research 5 (1356). https://doi.org/10.12688/f1000research.8900.2.
* Patro, Robert, Geet Duggal, Michael I Love, Rafael A Irizarry, and Carl Kingsford. 2017. “Salmon Provides Fast and Bias-Aware Quantification of Transcript Expression.” Nature Methods 14 (March). https://doi.org/10.1038/nmeth.4197.
* Robinson, Mark D, Davis J McCarthy, and Gordon K Smyth. 2010. “EdgeR: A Bioconductor Package for Differential Expression Analysis of Digital Gene Expression Data.” Bioinformatics 26 (1): 139–40.
* Love, Michael I., et al. Swimming Downstream: Statistical Analysis of Differential Transcript Usage Following Salmon Quantification. 7:952, F1000Research, 14 Sept. 2018. f1000research.com, https://f1000research.com/articles/7-952