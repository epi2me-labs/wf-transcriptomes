Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | wf-transcriptomes-report.html | a HTML report document detailing the primary findings of the workflow | aggregated |
| Per file read stats | fastq_ingress_results/reads/fastcat_stats/per-file-stats.tsv | A TSV with per file read stats, including all samples. | aggregated |
| Read stats | fastq_ingress_results/reads/fastcat_stats/per-read-stats.tsv | A TSV with per read stats, including all samples. | aggregated |
| Run ID's | fastq_ingress_results/reads/fastcat_stats/run_ids | List of run IDs present in reads. | aggregated |
| Meta map json | fastq_ingress_results/reads/metamap.json | Metadata used in workflow presented in a JSON. | aggregated |
| Concatenated sequence data | fastq_ingress_results/reads/{{ alias }}.fastq.gz | Per sample reads concatenated in to one FASTQ file. | per-sample |
| Assembled transcriptome | {{ alias }}_transcriptome.fas | Per sample assembled transcriptome. | per-sample |
| Annotated assembled transcriptome | {{ alias }}_merged_transcriptome.fas | Per sample annotated assembled transcriptome. | per-sample |
| Alignment summary statistics | {{ alias }}_read_aln_stats.tsv | Per sample alignment summary statistics. | per-sample |
| GFF compare results. | {{ alias }}_gffcompare | All GFF compare output files. | per-sample |
| Differential gene expression results | de_analysis/results_dge.tsv | This is a gene-level result file that describes genes and their probability of showing differential expression between experimental conditions. | aggregated |
| Differential gene expression report | de_analysis/results_dge.pdf | Summary report of differential gene expression analysis as a PDF. | aggregated |
| Differential transcript usage gene TSV | de_analysis/results_dtu_gene.tsv | This is a gene-level result file from DEXSeq that lists annotated genes and their probabilities of differential expression. | aggregated |
| Differential transcript usage report | de_analysis/results_dtu.pdf | Summary report of differential transcript usage results as a PDF. | aggregated |
| Differential transcript usage TSV | de_analysis/results_dtu_transcript.tsv | This is a transcript-level result file from DEXSeq that lists annotated genes and their probabilities of differential expression. | aggregated |
| Differential transcript usage stageR TSV | de_analysis/results_dtu_stageR.tsv  | This is the output from StageR and it shows both gene and transcript probabilities of differential expression | aggregated |
| Differential transcript usage DEXSeq TSV | de_analysis/results_dexseq.tsv | The complete output from the DEXSeq-analysis, shows both gene and transcript probabilities of differential expression. | aggregated |
| Gene counts | de_analysis/all_gene_counts.tsv | Raw gene counts created by the Salmon tool, before filtering. | aggregated |
| Gene counts per million | de_analysis/cpm_gene_counts.tsv | This file shows counts per million (CPM) of the raw gene counts to facilitate comparisons across samples. | aggregated |
| Transcript counts | de_analysis/unfiltered_transcript_counts_with_genes.tsv | Raw transcript counts created by the Salmon tool, before filtering. Includes reference to the associated gene ID. | aggregated |
| Transcript per million counts | de_analysis/unfiltered_tpm_transcript_counts.tsv | This file shows transcripts per million (TPM) of the raw counts to facilitate comparisons across samples. | aggregated |
| Transcript counts filtered | de_analysis/filtered_transcript_counts_with_genes.tsv | Filtered transcript counts, used for differential transcript usage analysis. Includes a reference to the associated gene ID. | aggregated |
| Transcript info table | {{ alias }}_transcripts_table.tsv | This file details each isoform that was reconstructed from the input reads. It contains a subset of columns from the .tmap output from [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) | per-sample |
| Final non redundant transcriptome | de_analysis/final_non_redundant_transcriptome.fasta | Transcripts that were used for differential expression analysis including novel transcripts with the identifiers used for DE analysis. | aggregated |
| Index of reference FASTA file | igv_reference/{{ ref_genome file }}.fai | Reference genome index of the FASTA file required for IGV config. | aggregated |
| GZI index of the reference FASTA file | igv_reference/{{ ref_genome file }}.gzi | GZI Index of the reference FASTA file. | aggregated |
| JSON configuration file for IGV browser | igv.json | JSON configuration file to be loaded in IGV for visualising alignments against the reference. | aggregated |
| BAM file (minimap2) | BAMS/{{ alias }}.reads_aln_sorted.bam | BAM file generated from mapping input reads to the reference. | per-sample |
| BAM index file (minimap2) | BAMS/{{ alias }}.reads_aln_sort.bam.bai | Index file generated from mapping input reads to the reference. | per-sample |
