Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | wf-transcriptomes-report.html | HTML report summarising transcript discovery, quantification, optional SQANTI3 classification, and optional differential analysis results. | aggregated |
| Per-file read stats | ingress_results/{{ alias }}/fastcat_stats/per-file-stats.tsv | Read statistics for each input FASTQ file in a sample, when FASTQ read stats are available. | per-sample |
| Per-read stats | ingress_results/{{ alias }}/fastcat_stats/per-read-stats.tsv.gz | Read statistics for individual reads in a sample, when this output is enabled. | per-sample |
| Ingress reads | ingress_results/{{ alias }}/seqs.fastq.gz | Reads prepared from the input data for downstream analysis. | per-sample |
| Ingress metadata | ingress_results/{{ alias }}/metamap.json | Per-sample metadata used by the workflow. | per-sample |
| Aligned BAM | cohort/alignments/{{ alias }}.aligned.sorted.bam | Genome-aligned BAM used for bambu, optional SQANTI3 QC, and IGV. | per-sample |
| Aligned BAM index | cohort/alignments/{{ alias }}.aligned.sorted.bam.bai | Index for the aligned BAM. | per-sample |
| Alignment summary | cohort/alignments/{{ alias }}.flagstat.txt | samtools flagstat output for the aligned BAM. | per-sample |
| Cohort transcriptome GTF | cohort/transcripts.gtf | Joint bambu transcript model used as the primary cohort transcriptome. | aggregated |
| Cohort transcriptome FASTA | cohort/cohort.transcriptome.fa | Transcript sequences derived from the joint cohort GTF. | aggregated |
| Cohort transcript counts | cohort/transcript_counts.tsv | Transcript-level count matrix produced by bambu. | aggregated |
| Cohort gene counts | cohort/gene_counts.tsv | Gene-level count matrix derived from bambu output. | aggregated |
| Cohort transcript metadata | cohort/transcript_metadata.tsv | Transcript annotations and bambu transcript classes for the cohort model. | aggregated |
| Cohort SQANTI3 summary | cohort/sqanti_cohort/classification_summary.tsv | SQANTI3 classification summary for the cohort transcriptome when SQANTI3 QC is enabled. | aggregated |
| Per-sample transcriptome GTF | samples/{{ alias }}/transcripts.gtf | Independent bambu transcript model for an individual sample. | per-sample |
| Per-sample transcriptome FASTA | samples/{{ alias }}/{{ alias }}.transcriptome.fa | Transcript sequences derived from the per-sample GTF. | per-sample |
| Per-sample transcript counts | samples/{{ alias }}/transcript_counts.tsv | Transcript-level abundance estimates for the per-sample bambu model. | per-sample |
| Per-sample gene counts | samples/{{ alias }}/gene_counts.tsv | Gene-level abundance estimates for the per-sample bambu model. | per-sample |
| Per-sample transcript metadata | samples/{{ alias }}/transcript_metadata.tsv | Transcript annotations and bambu transcript classes for the per-sample model. | per-sample |
| Per-sample SQANTI3 summary | samples/{{ alias }}/{{ alias }}_sqanti/classification_summary.tsv | SQANTI3 classification summary for the per-sample transcriptome when SQANTI3 QC is enabled. | per-sample |
| Differential gene expression results | de_analysis/{{ contrast }}/results_dge.tsv | DESeq2 gene-level differential expression results for one contrast. | aggregated |
| Differential gene expression plots | de_analysis/{{ contrast }}/results_dge.pdf | PDF plots generated during DESeq2 analysis for one contrast. | aggregated |
| Differential transcript usage results | de_analysis/{{ contrast }}/results_dtu_transcript.tsv | Transcript-level DTU results for one contrast. | aggregated |
| Differential transcript usage gene summary | de_analysis/{{ contrast }}/results_dtu_gene.tsv | Gene-level DTU summary for one contrast. | aggregated |
| DEXSeq results | de_analysis/{{ contrast }}/results_dexseq.tsv | Full DEXSeq result table for one contrast. | aggregated |
| Differential transcript usage plots | de_analysis/{{ contrast }}/results_dtu.pdf | PDF plots generated during DEXSeq analysis for one contrast. | aggregated |
| IGV configuration | igv.json | JSON configuration for viewing the aligned BAMs in IGV. | aggregated |
| Reference FASTA index | igv_reference/{{ ref_genome_file }}.fai | FAI index for the reference genome published for IGV. | aggregated |
| Reference GZI index | igv_reference/{{ ref_genome_file }}.gzi | GZI index for a compressed reference genome published for IGV. | aggregated |
