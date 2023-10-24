#!/usr/bin/env Rscript

# script to install Bioconductor R packages required for generation of reports by CRISPR-pooled-RSL pipeline
# based on https://github.com/casbap/ncRNA/blob/main/docker/rpkgs.R
# used in images based on rocker/tidyverse:4.2.3

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(version = "3.16")

# BiocManager::install(c('MAGeCKFlute','edgeR','fgsea','clusterProfiler','enrichplot','DOSE','biomaRt','ReactomePA','reactome.db','org.Hs.eg.db','org.Mm.eg.db'),
# 	version = "3.16", ask=FALSE, update=FALSE)

BiocManager::install(c('Rsamtools','GenomicAlignments','BiocParallel'),
	ask=FALSE, update=FALSE)

print("Install Bioconductor packages, done!")
