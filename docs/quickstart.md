## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop),
[Singularity](https://sylabs.io/singularity/) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Each method is automated out-of-the-box provided
either docker, singularity or conda is installed.

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
nexflow run epi2me-labs/wf-transcriptomes --fastq ERR6053095_chr20.fastq --ref_genome chr20/hg38_chr20.fa --ref_annotation chr20/gencode.v22.annotation.chr20.gtf \
      --jaffal_refBase chr20/ --jaffal_genome hg38_chr20 --jaffal_annotation genCode22" --out_dir outdir -w workspace_dir -profile conda -resume
```

**Example workflow for denovo transcript assembly**
```
OUTPUT=~/output
nextflow run . --fastq test_data/fastq --denovo --ref_genome test_data/SIRV_150601a.fasta  -profile local --out_dir ${OUTPUT} -w ${OUTPUT}/workspace \
--sample sample_id -resume
```
A full list of options can be seen in nextflow_schema.json. Below are some commonly used ones.

- Threshold for including isoforms into interactive table `transcript_table_cov_thresh = 50`
- Run the denovo pipeline `denovo = true` (default false)
- To run the workflow with direct RNA reads `--direct_rna` (this just skips the pychopper step).


Pychopper and minimap2 can take options via `minimap2_opts` and `pychopper_opts`, for example:


- When using the SIRV synthetic test data  
  - `minimap2_opts = '-uf --splice-flank=no'`
- pychopper needs to know which cDNA synthesis kit used
  - SQK-PCS109: use `pychopper_opts = '-k PCS109'` (default)
  - SQK-PCS110: use `pychopper_opts = '-k PCS110'`
  - SQK-PCS11:  use `pychopper_opts = '-k PCS111'`
- pychopper can use one of two available backends for identifying primers in the raw reads
  - nhmmscan `pychopper opts = '-m phmm'` 
  - edlib `pychopper opts = '-m edlib'`

__Note__: edlib is set by default in the config as it's quite a lot faster. However, it may be less sensitive than nhmmscan. 

### Fusion detection

JAFFAL from the [JAFFA](https://github.com/Oshlack/JAFFA)
package is used to identify potential fusion transcripts. To get this this working, there are a couple of things that need doing first.

**Install JAFFA**

to install JAFFA and it's dependencies run the folllowing:
```shell
cd wf-transcriptomes/
./subworkflows/JAFFAL/install_jaffa.sh
```

**Prepare JAFFAL reference data**

To use pre-processed reference files for the hg38 genome and GENCODE v22 annotation (as used in the JFFAAL paper),
do:
```shell
mkdir jaffal_data_dir
cd jaffal_data_dir/
wf-transcriptomes/subworkflows/JAFFAL/load_jaffal_references.sh
````

To use alternative genome and annotation files, they should be prepared as described
[here](https://github.com/Oshlack/JAFFA/wiki/FAQandTroubleshooting#how-can-i-generate-the-reference-files-for-a-non-supported-genome)

**Specifying the location of the JAFFA code and reference directories**

`--jaffal_dir`
Full path to the directory made by running install_jaffa.sh as shown above. eg: /home/wf-trnascriptomes/JAFFA

`--jaffal_refBase`
The directory containing the reference data prepared for use with JAFFAL


**JAFFAL annotation and genome files**

The prepared JAFFAL reference files will look something like `hg38_chr20_genCode22.fa`. To enable JAFFAL to find these
files `--jaffal_genome` should be set to `hg38_chr20` and `--jaffal_annotation` to `genCode22`


__JAFFAL Notes__: 
g++ must be installed. JAFFAL is not currently working on Mac M1 (osx-arm64 architecture). If there are no fusion transcripts
detected, the workflow will terminate with an error at the JAFFAL stage. If this happens, 
skip the JAFFAL stage by omitting ` --jaffal_refBase`

### Differential Expression

Differential Expression requires at least 2 replicates of each sample to compare. You can see an example condition_sheet.tsv in test_data.  

**Example workflow for differential expression transcript assembly**

#### Condition sheet
The condition sheet should be a .tsv with two columns.
- The sample column will need to match the 6 directories in the input fastq directory.
- The condition column will need to contain one of two keys to indicate the two samples being compared.

In the default `condition_sheet.tsv` available in the test_data directory we have used the following.

eg. condition_sheet.tsv
```
sample,condition
barcode01,untreated
barcode02,untreated
barcode03,untreated
barcode04,treated
barcode05,treated
barcode06,treated
```

You will also need to provide a reference genome and a reference annotation file.
Here is an example cmd to run the workflow using the test_data provided. 

```
OUTPUT=~/output;
nexflow run epi2me-labs/wf-transcriptomes --fastq test_data/differential_expression_fastq \
  --de_analysis \
  --ref_genome test_data/hg38_chr20.fa \
  --ref_annotation test_data/gencode.v22.annotation.chr20.gtf \
  --direct_rna
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