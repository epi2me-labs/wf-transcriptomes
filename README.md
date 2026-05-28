# Transcriptomes

Long-read transcript discovery, quantification, differential expression, QC and mod counting.



## Introduction

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




## Compute requirements

Recommended requirements:

+ CPUs = 32
+ Memory = 96GB

Minimum requirements:

+ CPUs = 12
+ Memory = 64GB

Approximate run time: Varies with read depth and sample count; a small (1-5 million reads) single-sample run may finish in under 60 minutes with the recommended resources. Compute requirements are influenced most by total read count, reference complexity, number of samples, and whether optional steps such as DE/DTU, and modified-base summarisation are enabled.

ARM processor support: False




## Install and run


These are instructions to install and run the workflow on command line.
You can also access the workflow via the
[EPI2ME Desktop application](https://epi2me.nanoporetech.com/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using either
[Docker](https://docs.docker.com/get-started/)
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
be found in the
[documentation](https://epi2me.nanoporetech.com/epi2me-docs/wfquickstart/).

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




## Pipeline overview

### Background

The methodology implemented within the wf-transcriptomes workflow follows from the largest independent long-read RNA benchmark to date.
The [Systematic assessment of long-read RNA-seq methods for transcript identification and quantification](https://www.nature.com/articles/s41592-024-02298-3) concluded that, in well-annotated genomes, reference-based methods perform best.
Our own previous research, benchmarking, and support of community members has shown that an automated, hands-off de-novo discovery pipeline to be bothersome for many use cases.
The wf-transcriptomes workflow therefore focuses on a reference-guided approach rather than a novelty-first one.

The benchmark paper above explicitly recommends `bambu` for identifying sample-specific transcriptomes in well-annotated organisms when only limited novelty is expected.
The paper also names `bambu` as one of the best options when quantification is important, which supports using it as the core engine for downstream DGE and DTU analyses.
For method details and current implementation notes, see the [`bambu` GitHub repository](https://github.com/GoekeLab/bambu).

In spike-in evaluations, `bambu` generally showed high precision and was among the better F1 performers.
This is an acceptable tradeoff for a production workflow where false transcript calls might be confound downstream analysis.
Users interested more in novel discovery may wish to amend the parameters of the workflow away from their defaults.
`bambu` also performed especially well on long non-spliced SIRVs, which supports its use on long-read datasets where transcript-end definition matters.

The workflow's choice of SQANTI3 as a companion QC and annotation layer matches the benchmark paper, which used SQANTI3 categories and metrics as its transcript assessment framework; so our outputs align with the field’s standard reporting language.
This step is useful for structural isoform classification and provides standard QC summaries for reporting, method comparison, or deeper transcript model review.
See the [`SQANTI3` repository](https://github.com/ConesaLab/SQANTI3) and the
[`SQANTI3` isoform classification documentation](https://github.com/ConesaLab/SQANTI3/wiki/)
for category definitions and usage details.



### 1. Getting your files into the workflow

The shared EPI2ME input handling collects FASTQ or BAM inputs, works out
whether you have a single sample or a multiplexed run, and produces per-sample
FASTQ files plus read statistics. These files are used in the downstream report.

### 2. Sample sheet formulation

The sample sheet is optional for simple single-sample runs, but it allows
sample aliases to be mapped to barcodes in multiplexed runs and is required for
`--de_analysis`.

+ Every row must contain `barcode` and `alias`.
+ `barcode` must use the usual ONT-style naming such as `barcode01`,
  `barcode02`, and the values must be unique.
+ `alias` is the user-facing sample name, must be unique, must not begin
  with the word `barcode` and may contain only letters, numbers, `.`, `_` or `-`.
+ If a `type` column is present, it must use one of:
  `test_sample`, `positive_control`, `negative_control`, or
  `no_template_control`.
+ If an `analysis_group` column is present, every row must have a value.
+ For `--de_analysis`, the sheet must also contain the primary condition
  column (`condition` by default, overridable with `--condition_column`),
  plus any columns named in `--covariates`.

When multiplexed input folders are named by barcode, the workflow matches those
folder names against the `barcode` column. If the folders are named by alias,
the workflow can match them against `alias`, but the sample sheet still needs a
`barcode` column because the shared validator expects it.

Example sample sheets:

#### Example sample sheet for a simple multiplexed analysis
This example contains the minimal `barcode` and `alias` columns:

```csv
barcode,alias
barcode01,rep1
barcode02,rep2
barcode03,rep3
barcode04,rep4
```


#### Example sample sheet for a differential expression analysis
This example is suitable for a multiplexed run and also satisfies the minimum
requirements for a two-group DE/DTU comparison; containing the `barcode`, `alias` 
and `condition` columns.

```csv
barcode,alias,condition
barcode01,control_rep1,control
barcode02,control_rep2,control
barcode03,control_rep3,control
barcode04,treated_rep1,treated
barcode05,treated_rep2,treated
barcode06,treated_rep3,treated
```
Additional columns to use as contrast facets may be named with the `--covariates` parameter.

### 3. Genome alignment

Each sample is aligned to the supplied reference genome with
[`minimap2`](https://github.com/lh3/minimap2) in a splice aware mode, then sorted and indexed with
[`samtools`](https://www.htslib.org/). The aligned BAMs under
`samples/<alias>/alignment/` are the main alignment files used for transcriptome
analysis, [`SQANTI3`](https://github.com/conesalab/SQANTI3) QC, and optional IGV viewing.

### 4. Optional modified base summarisation

When aligned BAMs contain modified base tags (`MM` and `ML`), the workflow also
runs `modkit` on each sample alignment. It first checks which modified base
codes are present in the BAM, then runs `modkit pileup` to produce a per-sample
bedMethyl file, a simple per-sample summary table, and one bigWig track per
requested or inferred modification under `samples/<alias>/mods/`.

If `--mod_codes` is set, those codes are passed directly to `modkit pileup`.
If it is omitted, the workflow infers the available `primary_base:mod_code`
pairs from the aligned BAM with `modkit modbam check-tags`. These outputs are
also included in the optional IGV configuration when `--igv` is enabled.

### 5. Cohort transcriptome construction

All aligned samples are analysed together with `bambu` to produce the primary
cohort transcriptome, transcript counts, gene counts, and the `RDS` objects used
for downstream differential analysis. This shared model is the main cohort-level
result and is published under `cohort/`.
Before writing outputs, transcript filtering removes only transcripts with zero
total transcript counts across samples; it does not use `fullLengthCounts` for
this quantification filter.

### 6. Independent per-sample transcriptomes

Each sample is also processed separately with `bambu` so the workflow produces
sample-specific GTF, FASTA, count tables, and metadata under
`samples/<alias>/`. These per-sample outputs are useful for inspecting sample
specific transcript models without changing the shared cohort transcriptome used
for DE/DTU.

### 7. Transcript sequence generation and QC

Transcript FASTA files are derived from GTF plus genome using `gffread`.
`SQANTI3` classifies the cohort and per-sample
transcriptomes and produces structural QC summaries. The cohort `SQANTI3`
results live under `cohort/sqanti/`, while per-sample `SQANTI3`
directories are published under `samples/<alias>/sqanti/`.

### 8. Optional DE and DTU analysis

When `--de_analysis` is enabled, the workflow checks the experimental design,
runs `DESeq2` for differential gene expression, and runs `DEXSeq` for
differential transcript usage. These analyses use the shared `bambu` outputs
and the design columns in the sample sheet, and each comparison is written to
its own subdirectory under `de_analysis/<contrast>/`.

### 9. What you need to provide

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





## Input parameters

### Main Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ reads to analyse. | You can provide a single FASTQ, a folder of FASTQs, or a multiplexed folder containing one sub-folder per sample or barcode. |  |
| bam | string | BAM or uBAM reads to analyse. | You can provide a single BAM or uBAM, a folder of BAMs, or a multiplexed folder containing one sub-folder per sample or barcode. |  |
| ref_genome | string | Reference genome FASTA. | Required in both discover and fixed_annotation modes. |  |
| ref_annotation | string | Reference transcript annotation in GTF or GFF format. | Required in both discover and fixed_annotation modes. |  |
| transcriptome_mode | string | How bambu should prepare the transcriptome model. | Use discover for reference-guided transcript discovery and quantification, or fixed_annotation for quantification only against the supplied annotation. | discover |
| direct_rna | boolean | Set this for direct RNA sequencing libraries. |  | False |


### Read Filtering Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| analyse_unclassified | boolean | Include unclassified reads from multiplexed input directories. |  | False |
| analyse_fail | boolean | Include fail reads from bam_fail and fastq_fail folders found in sample folders on the input path. |  | False |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | CSV file describing barcodes, aliases, and optional experimental design columns. | For multiplexed runs, the sample sheet should contain both barcode and alias. Additionally, for differential analysis, it must also contain the condition column, and any extra columns named in `--covariates`. |  |
| sample | string | Single sample name for singleplexed input or to restrict multiplexed analysis to one sample. |  |  |


### Differential Expression Analysis Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| de_analysis | boolean | Run differential gene expression and differential transcript usage analyses. |  | False |
| condition_column | string | Main comparison column in the sample sheet. |  | condition |
| covariates | string | Comma-separated extra sample-sheet columns to adjust for, for example batch. | Each listed name must exist as a column in the sample sheet. |  |
| reference_level | string | Baseline group for the main comparison column. | If omitted, the workflow will use control when that level exists. |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for user-facing workflow outputs. |  | output |
| igv | boolean | Generate an IGV configuration file for the aligned BAM outputs. |  | False |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Thread count to use for the core workflow processes. |  | 4 |
| mod_codes | string | Comma-separated modified base codes to pass to modkit pileup. | Provide values accepted by `modkit pileup --modified-bases`, for example `A:a,C:m`. If omitted, the workflow infers `primary_base:mod_code` pairs from the BAM with `modkit modbam check-tags`. |  |
| force_alignment | boolean | Force re-alignment of input BAM files. | Read alignment is skipped if the existing sequence names in the aligned BAM match the provided reference. Enable this if the existing alignments used incorrect minimap2 presets (e.g. missing --splice or direct RNA settings). | False |
| ndr | number | Optional bambu novel discovery rate override. | Lower values are more conservative (higher precision), while higher values are more permissive (higher novel-discovery sensitivity). See the [`bambu` repository](https://github.com/GoekeLab/bambu) for method details. |  |
| sqanti_skip_orf | boolean | Skip ORF prediction during SQANTI3 QC. |  | True |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | wf-transcriptomes-report.html | HTML report summarising transcript discovery, quantification, SQANTI3 classification, and optional differential analysis results. | aggregated |
| Aligned BAM | samples/{{ alias }}/alignment/reads.bam | Genome-aligned BAM used for bambu, SQANTI3 QC, and IGV. | per-sample |
| Aligned BAM index | samples/{{ alias }}/alignment/reads.bam.bai | Index for the aligned BAM. | per-sample |
| Alignment summary | samples/{{ alias }}/alignment/bamstats.flagstat.tsv | bamstats flagstat summary for the aligned BAM. | per-sample |
| Modified base pileup | samples/{{ alias }}/mods/{{ alias }}.mods.bedmethyl.gz | Per-sample modkit bedMethyl pileup generated from the aligned BAM when MM and ML tags are present. | per-sample |
| Modified base summary | samples/{{ alias }}/mods/{{ alias }}.mods.summary.tsv | Per-sample global modification-percent summary aggregated from the modkit bedMethyl pileup, with one row per modification code. | per-sample |
| Modified base bigWig | samples/{{ alias }}/mods/{{ alias }}.mods.*.bw | Per-sample modkit bigWig tracks generated from the aligned BAM, with one file per requested or inferred modification code. | per-sample |
| Reference and annotation preparation summary | cohort/reference/annotation_reference_summary.json | Summary of reference and annotation preparation, including seqname overlap, build/provider hints, and excluded unstranded annotation counts. | aggregated |
| Excluded unstranded annotation records | cohort/reference/unstranded_annotation.gtf | Full set of annotation records excluded because their strand was not '+' or '-'. Present only when unstranded records are found. | aggregated |
| Cohort transcriptome GTF | cohort/transcripts.gtf | Joint bambu transcript model used as the primary cohort transcriptome. | aggregated |
| Cohort transcriptome FASTA | cohort/cohort.transcriptome.fa | Transcript sequences derived from the joint cohort GTF. | aggregated |
| Cohort transcript counts | cohort/transcript_counts.tsv | Transcript-level count matrix produced by bambu. | aggregated |
| Cohort gene counts | cohort/gene_counts.tsv | Gene-level count matrix derived from bambu output. | aggregated |
| Cohort transcript metadata | cohort/transcript_metadata.tsv | Transcript annotations and bambu transcript classes for the cohort model. | aggregated |
| Cohort SQANTI3 summary | cohort/sqanti/classification_summary.tsv | SQANTI3 classification summary for the cohort transcriptome. | aggregated |
| Per-sample transcriptome GTF | samples/{{ alias }}/transcripts.gtf | Independent bambu transcript model for an individual sample. | per-sample |
| Per-sample transcriptome FASTA | samples/{{ alias }}/{{ alias }}.transcriptome.fa | Transcript sequences derived from the per-sample GTF. | per-sample |
| Per-sample transcript counts | samples/{{ alias }}/transcript_counts.tsv | Transcript-level abundance estimates for the per-sample bambu model. | per-sample |
| Per-sample gene counts | samples/{{ alias }}/gene_counts.tsv | Gene-level abundance estimates for the per-sample bambu model. | per-sample |
| Per-sample transcript metadata | samples/{{ alias }}/transcript_metadata.tsv | Transcript annotations and bambu transcript classes for the per-sample model. | per-sample |
| Per-sample SQANTI3 summary | samples/{{ alias }}/sqanti/classification_summary.tsv | SQANTI3 classification summary for the per-sample transcriptome. | per-sample |
| Differential gene expression results | de_analysis/{{ contrast }}/results_dge.tsv | DESeq2 gene-level differential expression results for one contrast. | aggregated |
| Differential gene expression plots | de_analysis/{{ contrast }}/results_dge.pdf | PDF plots generated during DESeq2 analysis for one contrast. | aggregated |
| Differential transcript usage results | de_analysis/{{ contrast }}/results_dtu_transcript.tsv | Transcript-level DTU results for one contrast. | aggregated |
| Differential transcript usage gene summary | de_analysis/{{ contrast }}/results_dtu_gene.tsv | Gene-level DTU summary for one contrast. | aggregated |
| DEXSeq results | de_analysis/{{ contrast }}/results_dexseq.tsv | Full DEXSeq result table for one contrast. | aggregated |
| Differential transcript usage plots | de_analysis/{{ contrast }}/results_dtu.pdf | PDF plots generated during DEXSeq analysis for one contrast. | aggregated |
| Differential analysis QC summary | de_analysis/de_qc_stats.json | Structured DE/DTU QC summary. Use analysis_fallbacks for aggregate counts, and each contrast's deseq2_dispersion_fallback, dexseq_dispersion_method, and dexseq_covariates_dropped fields for interpretation. | aggregated |
| Differential analysis text summary | de_analysis/de_overall_summary.txt | Human-readable DE/DTU run summary across all contrasts. | aggregated |
| Per-contrast QC summary | de_analysis/{{ contrast }}/contrast_qc_summary.txt | Human-readable per-contrast DE/DTU QC summary including sample counts and key significance totals. | aggregated |
| DESeq2 fallback diagnostic | de_analysis/DESeq2_dispersion_fallback_{{ contrast }}.txt | Diagnostic details when DESeq2 falls back to gene-wise dispersion estimation. | aggregated |
| DGE failure diagnostic | de_analysis/{{ contrast }}/DGE_ANALYSIS_FAILED.txt | Diagnostic details when DESeq2 fails for a contrast. | aggregated |
| DTU failure diagnostic | de_analysis/{{ contrast }}/DTU_ANALYSIS_FAILED.txt | Diagnostic details when DEXSeq fails for a contrast. | aggregated |
| Multiple-testing warning | de_analysis/MULTIPLE_TESTING_WARNING.txt | Family-wise error-rate note generated when multiple contrasts are tested. | aggregated |
| IGV configuration | igv.json | JSON configuration for viewing the aligned BAMs in IGV. | aggregated |
| Reference FASTA index | reference/{{ ref_genome_file }}.fai | FAI index for the reference genome published for IGV. | aggregated |
| Reference GZI index | reference/{{ ref_genome_file }}.gzi | GZI index for a compressed reference genome published for IGV. | aggregated |




## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related RNA and cDNA sequencing protocols in the
[Nanopore community](https://community.nanoporetech.com/docs/).




## Troubleshooting

+ Check that the reference genome and reference annotation use overlapping
  sequence names. The workflow checks this early and will fail if there is no
  overlap at all.
+ If `--de_analysis` is enabled, ensure the sample sheet contains `alias`, the
  primary condition column, and any columns named in `--covariates`.
+ DE/DTU requires at least two condition levels and at least two samples per
  level.
+ See how to interpret common Nextflow exit codes
  [here](https://epi2me.nanoporetech.com/epi2me-docs/help/troubleshooting/).

### Common confusion when coming from the previous workflow version

#### I supplied `--ref_transcriptome`, but the workflow still built or used `bambu` outputs

The previous workflow version used `--ref_transcriptome` as a main driver for
transcript-level downstream analysis. In the current version the
`--ref_transcriptome` option has been removed and is not part of the current
`bambu` input setup.

Use `--transcriptome_mode fixed_annotation` together with
`--ref_genome` and `--ref_annotation` if you want annotation-driven quantification
without transcript discovery.

#### I omitted `--ref_annotation` because I expected the old input setup

The previous workflow version could be driven from a different combination of
transcriptome inputs. In the current version both `--ref_genome` and
`--ref_annotation` are required in `discover` and `fixed_annotation` modes.
Always provide a compatible genome FASTA and transcript annotation when
launching the workflow.

#### I used `--transcriptome_source` and got behaviour I did not expect

In the previous workflow version the `--transcriptome_source` was the main
setting that chose how the workflow behaved. In the current version,
`--transcriptome_mode` is the main setting that controls this. The
`--transcriptome_source` option has been removed from the workflow interface.

The options `--transcriptome_mode discover` or `--transcriptome_mode
fixed_annotation` should be used to choose between the modes of operation.

#### I cannot find the old flat DE output files

In the previous workflow version, DE files appeared directly under `de_analysis/`.
In the current version DE and DTU results are grouped by contrast under
`de_analysis/<contrast>/`. Look for outputs such as
`de_analysis/<contrast>/results_dge.tsv` and
`de_analysis/<contrast>/results_dtu_transcript.tsv`.

#### I expected the old output layout or transcriptome files

The previous workflow version emitted a single folder of all analysis outputs.
In the current version, the output folders are organised around
`cohort/`, `samples/<alias>/` and `de_analysis/<contrast>/`.

The previous workflow output a non-redundant transcriptome whereas
the current version of the workflow outputs a true joint transcriptome
in the `cohort/` folder, with individual transcriptome analyses
additionally output in the `samples/<alias>/` folders.

#### My DE/DTU run fails because of `--sample_sheet`

The sample sheet rules have been amended compared to the previous version
to allow for multi-way comparisons.

The sample sheet must contain `alias`, the primary condition
column, and any columns named in `--covariates`. At least two condition levels
are required, and each level must contain at least two samples.

What to change: verify the sample sheet columns first, then check
`--condition_column`, `--covariates`, and `--reference_level`.

#### I hit genome/annotation validation or strand-related annotation warnings

The workflow validates that the annotation and genome share
sequence names, and it excludes unstranded annotation entries from the
differential analysis path.

Confirm the genome and annotation come from a compatible source,
and ensure the annotation uses only `+` or `-` strand values where required for
DE/DTU and `SQANTI3`.




## FAQs

### Does the workflow support both cDNA and direct RNA?

Yes. Use `--direct_rna` for direct RNA data. cDNA is the default mode.

### Do I need both `--ref_genome` and `--ref_annotation` in fixed-annotation mode?

Yes. `bambu` still uses the genome together with the imported annotation.

### Does the workflow create one shared transcriptome or one per sample?

Both. The joint cohort model is the primary result for reporting and DE/DTU, and
the workflow also emits independent per-sample transcriptomes.

### Can I run DE/DTU without transcript discovery?

Yes. Use `--transcriptome_mode fixed_annotation` together with `--de_analysis`.

### What changed from the previous workflow version?

The current workflow uses `bambu`, `SQANTI3`, `DESeq2`, and `DEXSeq`.
The main transcriptome result is now one shared `bambu` model built from all
samples together, with separate per-sample `bambu` outputs published alongside
it. The most important differences are summarised in this FAQ section.

#### Why do the results now mention `bambu` and `SQANTI3` instead of StringTie or GffCompare?

The workflow now uses a different set of transcript analysis tools. It builds
its transcript models with `bambu`, classifies them with `SQANTI3`,
and performs DE/DTU from the cohort `bambu` outputs.

#### Why is `--ref_annotation` now required even in fixed-annotation mode?

`bambu` still requires the annotation together with the genome in both
`discover` and `fixed_annotation` modes. Fixed-annotation mode means
annotation-driven quantification, not annotation-only execution.

#### Can I still use `--ref_transcriptome`?

No. `--ref_transcriptome` has been removed from the workflow interface.

Short old-to-new example:

```text
Previous workflow version:
  --transcriptome_source precomputed --ref_transcriptome transcripts.fa

Current workflow:
  --transcriptome_mode fixed_annotation \
  --ref_genome genome.fa \
  --ref_annotation annotation.gtf
```

#### Why do I now get both cohort and per-sample transcriptomes?

The workflow now treats both as important outputs. The shared cohort model is
the main transcriptome used for reporting and optional DE/DTU, while the
per-sample transcriptomes are provided for looking at each sample separately.

#### Why are DE results under `de_analysis/<contrast>/`?

The workflow now writes one subdirectory per comparison instead of publishing
one flat DE result set. This makes the output folder clearer when more than one
comparison is present.

Short old-to-new example:

```text
Previous workflow version:
  de_analysis/results_dge.tsv

Current workflow:
  de_analysis/<contrast>/results_dge.tsv
```

#### Can I still run fixed-annotation quantification without transcript discovery?

Yes. Use `--transcriptome_mode fixed_annotation` together with `--ref_genome`,
`--ref_annotation`, and any optional DE/DTU settings.

#### What should I expect to differ in report contents and output filenames?

Expect the report and output folder to emphasise:

+ the joint cohort `bambu` transcriptome under `cohort/`
+ the per-sample `bambu` transcriptomes under `samples/<alias>/`
+ `SQANTI3` results under cohort and per-sample directories
+ contrast-specific DE/DTU outputs under `de_analysis/<contrast>/`

If your question is not answered here, please start a discussion on the
[community](https://community.nanoporetech.com/).




## Related blog posts



+ [Transcriptomic reference sets](https://epi2me.nanoporetech.com/transcriptomic-reference/).

+ See the [EPI2ME website](https://epi2me.nanoporetech.com/) for more workflow resources and blog posts.

### References

[Systematic assessment of long-read RNA-seq methods for transcript identification and quantification
](https://www.nature.com/articles/s41592-024-02298-3)



