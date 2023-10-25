#!/usr/bin/env Rscript

suppressMessages(library("DRIMSeq"))
suppressMessages(library("GenomicFeatures"))
args <- commandArgs(trailingOnly=TRUE)
ref_annotation <- args[1]
min_samps_gene_expr <- args[2]
min_samps_feature_expr <- args[3]
min_gene_expr <- args[4] 
min_feature_expr <- args[5]
annotation_type <- args[6]
strip_version <- args[7]

cat("Loading counts, conditions and parameters.\n")
cts <- as.matrix(read.csv("merged/all_counts.tsv", sep="\t", row.names="Reference", stringsAsFactors=FALSE))

# Set up sample data frame:
#changed this to sample_id
coldata <- read.csv("de_analysis/coldata.tsv", row.names="alias", sep=",", stringsAsFactors=TRUE)

coldata$sample_id <- rownames(coldata)
# check if control condition exists, sets as reference 
if(!"control" %in% coldata$condition)
  stop("sample_sheet.csv does not contain 'control' 
       condition - unable to set reference.")
coldata$condition <- relevel(coldata$condition, ref = "control")

cat("Loading annotation database.\n")

txdb <- makeTxDbFromGFF(ref_annotation,  format = annotation_type)
txdf <- select(txdb, keys(txdb,"GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx<- tab[match(txdf$GENEID, names(tab))]


if (strip_version == "true"){
  rownames(cts) <- lapply(rownames(cts),  sub, pattern = "\\.\\d+$", replacement = "")
}

cts <- cts[rownames(cts) %in% txdf$TXNAME, ] # FIXME: filter for transcripts which are in the annotation. Why they are not all there? 

# Reorder transcript/gene database to match input counts:
txdf <- txdf[match(rownames(cts), txdf$TXNAME), ]
rownames(txdf) <- NULL

# Create counts data frame:
counts<-data.frame(gene_id=txdf$GENEID, feature_id=txdf$TXNAME, cts)

cat("Filtering counts using DRIMSeq.\n")

d <- dmDSdata(counts=counts, samples=coldata)
trs_cts_unfiltered <- counts(d)

d <- dmFilter(d, min_samps_gene_expr = min_samps_gene_expr, min_samps_feature_expr = min_samps_feature_expr,
        min_gene_expr = min_gene_expr, min_feature_expr = min_feature_expr)

cat("Building model matrix.\n")
design <- model.matrix(~condition, data=DRIMSeq::samples(d))



suppressMessages(library("dplyr"))

# Sum transcript counts into gene counts:
cat("Sum transcript counts into gene counts.\n")
trs_cts <- counts(d)
write.table(trs_cts, file="merged/all_counts_filtered.tsv",sep="\t")

gene_cts <- trs_cts_unfiltered %>% dplyr::select(c(1, 3:ncol(trs_cts)))  %>% group_by(gene_id) %>% summarise_all(tibble::lst(sum)) %>% data.frame()
rownames(gene_cts) <- gene_cts$gene_id
gene_cts$gene_id <- NULL
write.table(gene_cts, file="merged/all_gene_counts.tsv",sep="\t")

# Differential gene expression using edgeR:
suppressMessages(library("edgeR"))
cat("Running differential gene expression analysis using edgeR.\n")

y <- DGEList(gene_cts)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit)
edger_res <- topTags(qlf, n=nrow(y), sort.by="PValue")[[1]]

pdf("de_analysis/results_dge.pdf")
plotMD(qlf)
abline(h=c(-1,1), col="blue")
plotQLDisp(fit)

write.table(as.data.frame(edger_res), file="de_analysis/results_dge.tsv", sep="\t")

# Differential transcript usage using DEXSeq:
suppressMessages(library("DEXSeq"))
cat("Running differential transcript usage analysis using DEXSeq.\n")

sample.data<-DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data, sampleData=sample.data, design=~sample + exon + condition:exon, featureID=trs_cts$feature_id, groupID=trs_cts$gene_id)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd, reducedModel=~sample + exon)
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="condition")
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)

dev.off()
pdf("de_analysis/results_dtu.pdf")
plotMA(dxr, cex=0.8, alpha=0.05) 
plotDispEsts(dxd)

qval <- perGeneQValue(dxr) 
dxr.g<-data.frame(gene=names(qval), qval)
dxr.g <- dxr.g[order(dxr.g$qval),]

dxr_out <- as.data.frame(dxr[,c("featureID", "groupID", "pvalue")])
dxr_out <- dxr_out[order(dxr$pvalue),]

write.table(dxr.g, file="de_analysis/results_dtu_gene.tsv", sep="\t")
write.table(dxr_out, file="de_analysis/results_dtu_transcript.tsv", sep="\t")

# and writing out some of the DEXSeq metrics to accompany EPI2ME Labs tutorial
colnames(dxr)[grep("log2fold", colnames(dxr))] <- "log2fold"
MADTUdata <- data.frame(dxr)[order(dxr$padj),c("exonBaseMean", "log2fold", "pvalue", "padj")]
MADTUdata$exonBaseMean <- log2(MADTUdata$exonBaseMean)
colnames(MADTUdata)[which(colnames(MADTUdata)=="exonBaseMean")] <- "Log2MeanExon"
colnames(MADTUdata)[which(colnames(MADTUdata)=="log2fold")] <- "Log2FC"
write.table(MADTUdata, file="de_analysis/results_dexseq.tsv", sep="\t")

# stageR analysis of DEXSeq results:
cat("stageR analysis\n")
library(stageR)

cat("Running stageR analysis on the differential transcript usage results.\n")
pConfirmation <- matrix(dxr$pvalue, ncol=1)

dimnames(pConfirmation) <- list(dxr$featureID, "transcript")
pScreen <- qval
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
# note: the choice of 0.05 here means you can *only* threshold at 5% OFDR later
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.10)
suppressWarnings({dex.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE)})

# dex.padj <- dex.padj[,-1]
write.table(dex.padj, file="de_analysis/results_dtu_stageR.tsv", sep="\t")
