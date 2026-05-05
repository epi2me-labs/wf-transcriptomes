+ Check that the reference genome and reference annotation use overlapping
  sequence names. The workflow checks this early and will fail if there is no
  overlap at all.
+ If `--de_analysis` is enabled, ensure the sample sheet contains `alias`, the
  primary condition column, and any columns named in `--covariates`.
+ DE/DTU requires at least two condition levels and at least two samples per
  level.
+ `--cdna_preprocess` must not be combined with `--direct_rna`.
+ See how to interpret common Nextflow exit codes
  [here](https://labs.epi2me.io/trouble-shooting/).

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

#### I expected pychopper-style preprocessing for cDNA libraries

cDNA preprocessing is available through
`--cdna_preprocess`, which runs `pychopper` before alignment. The
parameters for that stage are `--cdna_kit`, `--pychopper_backend`, and
`--pychopper_opts`.

Enable `--cdna_preprocess` for cDNA libraries when you want
`pychopper` preprocessing, or use `--direct_rna` without
`--cdna_preprocess` for direct RNA data.

#### I cannot find the old flat DE output files

In the previous workflow version, DE files appeared directly under `de_analysis/`.
In the current version DE and DTU results are grouped by contrast under
`de_analysis/<contrast>/`. Look for outputs such as
`de_analysis/<contrast>/results_dge.tsv` and
`de_analysis/<contrast>/results_dtu_transcript.tsv`.

#### I expected the old output layout or transcriptome files

The previous workflow version emitted one flat set of transcriptome.
In the current versions, the output folder is organised around `ingress_results/`,
`cohort/`, `samples/<alias>/`, `de_analysis/<contrast>/`, and
`igv_reference/`.

Primary shared transcriptome results are under `cohort/`, and
`samples/<alias>/` for sample-specific models. If `--cdna_preprocess` is enabled,
look under `ingress_results/<alias>/` for the `pychopper` outputs associated
with each sample.

#### My DE/DTU run fails because of `--sample_sheet`

The sample sheet rules have been amended compared to the previous version
to allow for multi-way comparisons.

The sample sheet must contain `alias`, the primary condition
column, and any columns named in `--covariates`. At least two condition levels
are required, and each level must contain at least two samples.

What to change: verify the sample sheet columns first, then check
`--condition_column`, `--covariates`, and `--reference_level`.

#### I combined `--direct_rna` with `--cdna_preprocess`

This is not a valid combination of parameters. Use `--direct_rna` for direct RNA
libraries and leave `--cdna_preprocess` unset.

#### I hit genome/annotation validation or strand-related annotation warnings

The workflow validates that the annotation and genome share
sequence names, and it excludes unstranded annotation entries from the
differential analysis path.

Confirm the genome and annotation come from a compatible source,
and ensure the annotation uses only `+` or `-` strand values where required for
DE/DTU and `SQANTI3`.
