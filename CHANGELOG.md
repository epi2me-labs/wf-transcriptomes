# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Added
- Published minimap2 and pychopper results to output directory.
- Two extra pychopper parameters `--cdna_kit` and `--pychopper_backend`. `--pychopper_options` is still available to define any other options.

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

