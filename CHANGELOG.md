# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.1.5]
### Added
- Differential transcript and gene expression subworkflow

## [v0.1.4]
### Added
- JAFFAL fusion detectoion subworkflow
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
