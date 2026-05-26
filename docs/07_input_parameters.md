### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ reads to analyse. | You can provide a single FASTQ, a folder of FASTQs, or a multiplexed folder containing one sub-folder per sample or barcode. |  |
| bam | string | BAM or uBAM reads to analyse. | You can provide a single BAM or uBAM, a folder of BAMs, or a multiplexed folder containing one sub-folder per sample or barcode. |  |
| analyse_unclassified | boolean | Include unclassified reads from multiplexed input directories. |  | False |
| analyse_fail | boolean | Include fail reads from multiplexed input directories. |  | False |
| fastq_chunk | integer | Maximum number of reads per ingress chunk. | Useful mainly for testing or for splitting very large inputs into smaller pieces. |  |


### Reference Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| ref_genome | string | Reference genome FASTA. | Required in both discover and fixed_annotation modes. |  |
| ref_annotation | string | Reference transcript annotation in GTF or GFF format. | Required in both discover and fixed_annotation modes. |  |
| transcriptome_mode | string | How bambu should prepare the transcriptome model. | Use discover for reference-guided transcript discovery and quantification, or fixed_annotation for quantification only against the supplied annotation. | discover |
| direct_rna | boolean | Set this for direct RNA sequencing libraries. |  | False |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | CSV file describing barcodes, aliases, and optional experimental design columns. | For multiplexed runs, the sample sheet should contain both barcode and alias. For differential analysis it must also contain alias, the condition column, and any extra columns named in `--covariates`. |  |
| sample | string | Single sample name for singleplexed input or to restrict multiplexed analysis to one sample. |  |  |


### Analysis Options

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
| minimap2_opts | string | Extra command-line options to pass to minimap2. |  |  |
| force_alignment | boolean | Force re-alignment of input BAM files. | Read alignment is skipped if the existing sequence names in the aligned BAM match the provided reference. Enable this if the existing alignments used incorrect minimap2 presets (e.g. missing --splice or direct RNA settings). | False |
| ndr | number | Optional bambu novel discovery rate override. |  |  |
| skip_sqanti | boolean | Skip SQANTI3 transcript classification and QC. |  | False |
| sqanti_skip_orf | boolean | Skip ORF prediction during SQANTI3 QC. |  | True |
| sqanti_extra_args | string | Extra command-line options to pass to SQANTI3. |  |  |


