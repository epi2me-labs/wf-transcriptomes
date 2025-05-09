# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]
### Changed
- Updated to wf-template v5.6.1, changing:
    - Reduce verbosity of debug logging from fastcat which can occasionally occlude errors found in FASTQ files during ingress.
    - Log banner art to say "EPI2ME" instead of "EPI2ME Labs" to match current branding. This has no effect on the workflow outputs.
    - pre-commit configuration to resolve an internal dependency problem with flake8. This has no effect on the workflow.

### Fixed
- Updated to wf-template v5.6.1, fixing:
    - dacite.exceptions.WrongTypeError during report generation when barcode is null.
    - Sequence summary read length N50 incorrectly displayed minimum read length, it now correctly shows the N50.
    - Sequence summary component alignment and coverage plots failed to plot under some conditions.

## [v1.7.0]
### Changed
- `split_bam` and `build_minimap_index_transcriptome` process memory allocation increased.
- Updated recommended memory requirement.
- Updated project description.
- A common user issue is providing a ref_annotation and ref_genome parameter that have mismatched reference IDs, which causes the DE_analysis to fail. The workflow will now do an upfront check and give an error message if no overlap is found or a warning if some IDs are present in one file but not in the other.
- Reconciled workflow with wf-template v5.5.0.
- Sort the columns and rows of the gene and transcript count files.
- DE_analysis alignment summary stats table no longer includes MAPQ or quality scores. MAPQ is not relevant for transcript alignment and quality scores are already available in the read summary section of the report. 
### Fixed
- `all_gene_counts.tsv` contained the DE counts results.
- Reduced memory usage of the report workflow process.
- Output BAM alignments in all cases unless the workflow is run with `transcriptome_source` set to `precomputed`.
- Corrected the demo command in the `README.md`.
- The merged transcriptome generated for differential expression analysis now only contains the exons and not the full genomic sequence.
- Output the gene name annotated differential expression analysis count files only.
- Only use full length reads in the differential expression analysis.

## [v1.6.1]
### Fixed
- `merge_gff_compare` failing with empty GFF files.

## [v1.6.0]
### Fixed
- v1.5.0 bug; access to undefined channel output bug when using precomputed transcriptome.
- Bug where incorrect gene_id assigned in the DE tables.

## [v1.5.0]
### Updated
- Workflow report updated to use `ezcharts`.
### Fixed
- Exons per isoforms histogram reporting incorrect numbers.
- Output the `results_dexseq.tsv` file when `--de_analysis` enabled.
### Removed
- per-class gffcompare tracking files as there exists a combine tracking file. 

## [v1.4.0]
## Added
- `--igv` parameter (default: false) for outputting IGV config allowing visualisation of read alignments in the EPI2ME App.
- If required for IGV, reference indexes are output in to a `igv_reference` directory
### Changed
- BAMS are output in to a BAMS directory.
- Reconcile with template 5.2.6.

## [v1.3.0]
### Removed
- Fusion detection subworkflow, as the functionality is not robust enough for general use at this time.
### Changed
- Updated pychopper to 2.7.10
## Added 
- new `cdna_kit` options: PCS114 and PCB111/114

## [v1.2.1]
### Changed
- Increase some memory and CPU allocations.

## [v1.2.0]
### Added
- Workflow now accepts BAM or FASTQ files as input (using the --bam or --fastq parameters, respectively).
### Changed
- MA plot in the `results_dge.pdf` has been updated to match the MA plot in the report.
### Added
- Error message when running in `de_analysis` mode and `ref_annotation` input file contains unstranded annotations.

## [v1.1.1]
### Changed
- Improved handling of different annotation file types (eg. `.gtf/.gff/.gff3`) in `de_analysis` mode.
- Improved handling of annotation files that do not contain version numbers in transcript_id (such as gtf's from Ensembl).
### Fixed
- Differential expression failing with 10 or more samples.
- Regression causing the DE analysis numeric parameters to not be evaluated correctly.

## [v1.1.0]
### Changed
- Improve documentation around filtering of transcripts done before DTU analysis.
- Renamed files:
  -  `de_analysis/all_counts_filtered.tsv` to `de_analysis/filtered_transcript_counts_with_genes.tsv`
  -  `de_analysis/de_tpm_transcript_counts.tsv` to `de_analysis/unfiltered_tpm_transcript_counts.tsv`
- Minimum memory requirements to `32 GB`.
### Added
- Published isoforms table to output directory.
- Output additional `de_analysis/cpm_gene_counts.tsv` with counts per million gene counts.
- Output additional `de_analysis/unfiltered_transcript_counts_with_genes.tsv` with unfiltered transcript counts with associated gene IDs.
- Add gene name column to the de_analysis counts TSV files.
### Fixed
- Mapping stage using a single thread only.
### Changed
- More memory assigned to the fusion detection process.
- When no `--ref_annotation` is provided the workflow will still run but the output transcripts will not be annotated. However `--de_analysis` mode still requires a `--ref_annotation`.

## [v1.0.0]
### Added
- Published minimap2 and pychopper results to output directory.
- Two extra pychopper parameters `--cdna_kit` and `--pychopper_backend`. `--pychopper_options` is still available to define any other options.
- Memory requirements for each process.
### Changed
- Documentation.
### Fixed
- When Jaffa is run only output one report.

## [v0.4.2]
### Changed
- Sample sheet must include a `control` type to indicate which samples are the reference for the differential expression pipeline.
### Removed
- Default local executor CPU and RAM limits.

## [v0.4.1]
### Changed
- Updated docker container with Pychopper to support LSK114.

## [v0.4.0]
### Fixed
- Remove dead links from README
### Removed
- Denovo `--transcriptome_source` option.

## [v0.3.1]
### Added
- Handling for input reference transcriptome headers that contain `|`

## [v0.3.0]
### Changed
- Improve differential expression outputs.
- Include transcript and gene count tables in DE_final folder.
- If differential expression subworkflow is used a non redundant transcriptome will be output which includes novel transcripts.
- Added wording to the report about how to identify novel transcripts in the DE tables.
- Nextflow minimum required version to 23.04.2
- `--minimap_index_opts` parameter has been changed to `minimap2_index_opts` for consistency.

### Added
- An additional gene name column to the differential gene expression results. This is especially handy for transcriptomes where the gene ID is not the same as gene name (e.g. Ensembl).
- Wording to the report about how to identify novel transcripts in the DE tables.

## [v0.2.1]
### Changed
- Any sample aliases that contain spaces will be replaced with underscores.
- Updated documentation to explain we only support Ensembl, NCBI and ENCODE annotation file types. 

### Fixed
- Documentation parameter examples corrected.
- Handling for annotation files that use gene as gene_id attribute.
- Handling for Ensembl annotation files.

## [v0.2.0]
### Changed
- GitHub issue templates
- Condition sheet is no longer required. The sample sheet is now used to indicate condition instead.
    - For differential expression, the sample sheet must have a `condition` column to indicate which condition group each sample in the sample sheet belongs to.
    - Values for the condition may be any two distinct strings, for example: treated/untreated; sample/control etc.

### Fixed
- Remove default of null for `--ref_transcriptome`.
- Read mapping summary table in the report has correct sample_ids.

## [v0.1.13]
### Added
- Handling for GFF3 reference_annotation file type.
- Warning for the `--transcriptome_source` denovo pipeline option.

### Changed
- Enum choices are enumerated in the `--help` output
- Enum choices are enumerated as part of the error message when a user has selected an invalid choice
- Bumped minimum required Nextflow version to 22.10.8

### Fixed
- Replaced `--threads` option in fastqingress with hardcoded values to remove warning about undefined `param.threads`
- Fix for the `--transcriptome_source` denovo pipeline option.

## [v0.1.12]
### Added
- Handling for GFF3 reference_annotation file type.
- Handling gzip input reference and annotation parameters.
- Handling for NCBI gtfs that contain some empty transcript ID fields.

## [v0.1.11]
### Changed
- LICENSE to Oxford Nanopore Technologies PLC. Public License Version 1.0.

### Added
- Configuration for running demo data in AWS

## [v0.1.10]
### Changed
- Condition sheet parameter description fixed to CSV
- Update fastqingress

## [v0.1.9]
### Changed
- Simplify JAFFAL docs

## [v0.1.8]
### Changed
- Description in manifest

## [v0.1.7]
### Changed
- `-profile conda` is no longer supported, users should use `-profile standard` (Docker) or `-profile singularity` instead
- `nextflow run epi2me-labs/wf-transcriptomes --version` will now print the workflow version number and exit
- Use parameter `--transcriptome-source` to define precalculated, reference-based or denovo

## [v0.1.6]
### Changed
- Removed sanitize option
- Reduce size of differential expression data.

### Added
- Improved DE explanation in docs
- Option to turn off transcript assembly steps with param transcript_assembly

### Fixed
- Fix JAFFAL terminating workflow when no fusions found.
- Error if condition sheet and sample sheet don't match.
- Failed to plot DE graphs when one of data sets is 0 length.

## [v0.1.5]
### Added
- Differential transcript and gene expression subworkflow

## [v0.1.4]
### Added
- JAFFAL fusion detection subworkflow

### Changed
- Args parser for fastqingress
- Set out_dir option type to ensure output is written to correct directory on Windows
- Skip unnecessary conversion to fasta from fastq
- Fastqingress metadata map
- Changed workflow name to wf-transcriptomes

## [v0.1.3]
### Changed
- Better help text on cli
- Use EPI2ME Labs-maintained version of pychopper

## [v0.1.2]
### Added
- direct_rna option
- Some extra error handling
- Minor report display improvements

## [v0.1.1]
### Fixed
- Incorrect numbers and of transcripts caused by merging gff files with same gene and transcript ids
- Error handling in de novo pipeline. Skip clusters in build_backbones that cause an isONclust2 error
- Several small fixes in report plotting

## [v0.1.0]
### Added
- Added the denovo pipeline

### Changed
- Updates to the report plots

## [v0.0.1]
### Added
- First release
- Initial port of Snakemake WF from https://github.com/nanoporetech/pipeline-nanopore-ref-isoforms

