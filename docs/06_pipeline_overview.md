### 0. Background.

The methodology implemented within the wf-transcriptomics workflow follows from the largest independent long-read RNA benchmark to date.
The [Systematic assessment of long-read RNA-seq methods for transcript identification and quantification](https://www.nature.com/articles/s41592-024-02298-3) concluded that, in well-annotated genomes, reference-based methods perform best.
Our own previous research, benchmarking, and support of community members has shown that an automated, hands-off de-novo discovery pipeline to be bothersome for many use cases.
The wf-transcriptomes workflow therefore focuses on a reference-guided approach rather than a novelty-first one.

The benchmark paper above explicitly recommends `bambu` for identifying sample-specific transcriptomes in well-annotated organisms when only limited novelty is expected.
The paper also names `bambu` as one of the best options when quantification is important, which supports using it as the core engine for downstream DGE and DTU analyses.

In spike-in evaluations, `bambu` generally showed high precision and was among the better F1 performers.
This is an acceptable tradeoff for a production workflow where false transcript calls might be confound downstream analysis.
Users interested more in novel discovery may wish to amend the parameters of the workflow away from their defaults.
`bambu` also performed especially well on long non-spliced SIRVs, which supports its use on long-read datasets where transcript-end definition matters.

The workflow's choice of SQANTI3 as a companion QC and annotation layer matches the benchmark paper, which used SQANTI3 categories and metrics as its transcript assessment framework; so our outputs align with the field’s standard reporting language.



### 1. Getting your files into the workflow

The shared EPI2ME input handling collects FASTQ or BAM inputs, works out
whether you have a single sample or a multiplexed run, and produces per-sample
FASTQ files plus read statistics. These files are used in the downstream report.

### 2. Sample sheet formulation

The sample sheet is optional for simple single-sample runs, but it becomes the
main source of sample names for multiplexed runs and is required for
`--de_analysis`.

+ Every row must contain `barcode` and `alias`.
+ `barcode` must use the usual ONT-style naming such as `barcode01`,
  `barcode02`, and the values must be unique.
+ `alias` is the user-facing sample name, must be unique, and must not begin
  with the word `barcode`.
+ If a `type` column is present, it must use one of:
  `test_sample`, `positive_control`, `negative_control`, or
  `no_template_control`.
+ If an `analysis_group` column is present, every row must have a value.
+ For `--de_analysis`, the sheet must also contain the primary condition
  column, `condition` by default, plus any columns named in `--covariates`.

When multiplexed input folders are named by barcode, the workflow matches those
folder names against the `barcode` column. If the folders are named by alias,
the workflow can match them against `alias`, but the sample sheet still needs a
`barcode` column because the shared validator expects it.

Example sample sheet:

```csv
barcode,alias,type,condition,batch
barcode01,control_rep1,test_sample,control,b1
barcode02,control_rep2,test_sample,control,b2
barcode03,treated_rep1,test_sample,treated,b1
barcode04,treated_rep2,test_sample,treated,b2
```

This example is suitable for a multiplexed run and also satisfies the minimum
requirements for a two-group DE/DTU comparison.

### 3. Genome alignment

Each sample is aligned to the supplied reference genome with
[`minimap2`](https://github.com/lh3/minimap2), then sorted and indexed with
[`samtools`](https://www.htslib.org/). The aligned BAMs under
`samples/<alias>/alignment/` are the main alignment files used for transcriptome
analysis, optional [`SQANTI3`](https://github.com/conesalab/SQANTI3) QC, and optional IGV viewing.

### 4. Cohort transcriptome construction

All aligned samples are analysed together with `bambu` to produce the primary
cohort transcriptome, transcript counts, gene counts, and the `RDS` objects used
for downstream differential analysis. This shared model is the main cohort-level
result and is published under `cohort/`.
Before writing outputs, transcript filtering removes only transcripts with zero
total transcript counts across samples; it does not use `fullLengthCounts` for
this quantification filter.

### 5. Independent per-sample transcriptomes

Each sample is also processed separately with `bambu` so the workflow produces
sample-specific GTF, FASTA, count tables, and metadata under
`samples/<alias>/`. These per-sample outputs are useful for inspecting sample
specific transcript models without changing the shared cohort transcriptome used
for DE/DTU.

### 6. Transcript sequence generation and QC

Transcript FASTA files are derived from GTF plus genome using `gffread`.
When `--skip_sqanti` is not set, `SQANTI3` classifies the cohort and per-sample
transcriptomes and produces structural QC summaries. The cohort `SQANTI3`
results live under `cohort/sqanti/`, while per-sample `SQANTI3`
directories are published under `samples/<alias>/sqanti/`.

### 7. Optional DE and DTU analysis

When `--de_analysis` is enabled, the workflow checks the experimental design,
runs `DESeq2` for differential gene expression, and runs `DEXSeq` for
differential transcript usage. These analyses use the shared `bambu` outputs
and the design columns in the sample sheet, and each comparison is written to
its own subdirectory under `de_analysis/<contrast>/`.

### 8. What you need to provide

The workflow's analysis is controlled by a user provided genome, annotation, and
`bambu` mode.

* use `--transcriptome_mode` to choose between `discover` and
  `fixed_annotation`
* `--transcriptome_source` has been removed; use `--transcriptome_mode` instead
* both `--ref_genome` and `--ref_annotation` are required in both modes
* `--ref_transcriptome` has been removed; if you want annotation-based
  quantification, use `--transcriptome_mode fixed_annotation` together with
  `--ref_genome` and `--ref_annotation`
* when `--de_analysis` is enabled, the sample sheet must contain `alias`, the
  primary condition column, and any requested columns named in `--covariates`

### 9. How to read the output folder

The published outputs are organised around a small number of top-level
directories:

+ `cohort/` contains the primary joint `bambu` transcriptome, count tables, and optional cohort `SQANTI3` outputs
+ `samples/<alias>/` contains alignments, independent per-sample `bambu` outputs and
  optional per-sample `SQANTI3` outputs
+ `de_analysis/<contrast>/` contains DE and DTU results for each contrast when
  differential analysis is enabled
+ `igv_reference/` contains the published reference indexes used for IGV
