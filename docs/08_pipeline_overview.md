### 1. Concatenate input files and generate per read stats.
The [fastcat](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow. It will also output per read stats including average read lengths and qualities.

### 2. Preprocess cDNA.
If input sequences are cDNA [Pychopper](https://github.com/epi2me-labs/pychopper) is used to orient, trim and rescue full length cDNA reads and associated statistics. If the `direct_rna` parameter is selected this step will be skipped.

### 3. Build transcriptome.
If the `transcriptome_source` parameter is "reference-guided" a transcriptome will be built for each sample as outlined below. If the `transcriptome_source` is "precomputed" and the `reference_transcriptome` parameter is provided the workflow will skip step 3.

#### 3.1 Align reads with reference genome.
The reference genome will be indexed and aligned using [Minimap2](https://github.com/lh3/minimap2). The output is sorted and converted to a BAM file using [Samtools](https://www.htslib.org/). Alignment stats are created from these using [Seqkit BAM](https://bioinf.shenwei.me/seqkit/usage/#bam).

#### 3.2 Chunk BAM
The aligned BAMs are split into chunks using the bundle_min_reads parameter (default: 50000).

#### 3.3 Assemble transcripts
[StringTie](https://ccb.jhu.edu/software/stringtie/) is then used to assemble the transcripts using the aligned segments in the chunked BAM files. The assembled transcript will be output as a [GFF file](https://www.ensembl.org/info/website/upload/gff3.html). If a `ref_annotation` file is provided this will also be included in the GFF.

#### 3.4 Merge Chunks
Transcript GFF files from the chunks with the same sample aliases will then be merged.

#### 3.5 Annnotate
[GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.html) is then used to compare query and reference annotations, merging records where appropriate and then annotating them. This also creates estimates of accuracy of the GFF files output in a stats file per sample.

#### 3.6 Create transcriptomes
[Gffread](https://github.com/gpertea/gffread) is used to create a transcriptome FASTA file from the final GFF as well as a merged transcriptome that includes annotations in the FASTA headers where available.

### 4. Find gene fusions
If gene fusion options are provided, fusion gene detection is performed using [JAFFA](https://github.com/Oshlack/JAFFA), with the JAFFAL extension. To enable this provide the gene fusion detection options: `jaffal_refBase`, `jaffal_genome` and `jaffal_annotation`.

### 5. Differential expression analysis

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

#### 5.1 Merge cross sample transcriptomes
If a `ref_transcriptome` is not provided, the transcriptomes created by the workflow will be used for DE analysis. To do this, the GFF outputs of GffCompare are merged using StringTie. A final non redundant FASTA file of the transcripts is created using the merged GFF file and the reference genome using seqkit.

#### 5.2 Create a final non redundant transcriptome
The reads from all the samples will be aligned with the final non redundant transcriptome using Minimap2 in a splice aware manner.

#### 5.3 Count genes and transcripts
[Salmon](https://github.com/COMBINE-lab/salmon) is used for transcript quantification, giving gene and transcript counts.

#### 5.4 Pre-filtering of quantitative data using DRIMSeq
[DRIMSeq](https://bioconductor.org/packages/release/bioc/html/DRIMSeq.html) is used to filter the transcript count data from the Salmon analysis. The filter step will be used to select for genes and transcripts that satisfy rules for the number of samples in which a gene or transcript must be observed, and minimum threshold levels for the number of observed reads. The parameters used for filtering are `min_samps_gene_expr`, `min_samps_feature_expr`, `min_gene_expr`, and `min_feature_expr`. 

#### 5.5 edgeR based differential expression analysis
A statistical analysis is first performed using [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) to identify the subset of differentially expressed genes. The filtered list of gene counts is used as input. A normalisation factor is calculated for each sequence library (using the default TMM method described by [McCarthy et al. (2012)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3378882/) for further details). The defined experimental design is used to calculate estimates of dispersion for each of the gene features. Statistical tests are calculated using the contrasts defined in the experimental design. The differentially expressed genes are corrected for false discovery (fdr) using the method of Benjamini & Hochberg ([Benjamini and Hochberg (1995)](https://www.jstor.org/stable/2346101))

#### 5.6 Differential transcript usage using DEXSeq
Differential transcript usage analysis is performed using the R [DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html) package ([Anders et al. (2012)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3460195/)). Similar to the edgeR package, DEXSeq estimates the variance between the biological replicates and applies generalised linear models for the statistical testing. The key difference is that the DEXSeq method looks for differences at the exon count level. DEXSeq uses the filtered transcript count data prepared earlier in this analysis. 

#### 5.7 StageR stage-wise analysis of DGE and DTU
The final component of this isoform analysis is a stage-wise statistical test using the R software package [stageR](https://bioconductor.org/packages/release/bioc/html/stageR.html)([Van den Berge and Clement (2018)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1277-0)). stageR uses (1) the raw p-values for DTU from the DEXSeq analysis in the previous section and (2) a false-discovery corrected set of p-values from testing whether individual genes contain at least one exon showing DTU. A hierarchical two-stage statistical testing evaluates the set of genes for DTU.



