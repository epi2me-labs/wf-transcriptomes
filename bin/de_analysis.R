#!/usr/bin/env Rscript

suppressMessages(library("DRIMSeq"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("edgeR"))
args <- commandArgs(trailingOnly=TRUE)
ref_annotation <- args[1]
min_samps_gene_expr <- as.numeric(args[2])
min_samps_feature_expr <- as.numeric(args[3])
min_gene_expr <- as.numeric(args[4]) 
min_feature_expr <- as.numeric(args[5])

cat("Loading counts, conditions and parameters.\n")
cts <- as.matrix(read.csv("all_counts.tsv", sep="\t", row.names="Reference", stringsAsFactors=FALSE))

# Set up sample data frame:
#changed this to sample_id
coldata <- read.csv("sample_sheet.csv", row.names="alias", sep=",", stringsAsFactors=TRUE)

coldata$sample_id <- rownames(coldata)
# check if control condition exists, sets as reference 
if(!"control" %in% coldata$condition)
  stop("sample_sheet.csv does not contain 'control' 
       condition - unable to set reference.")
coldata$condition <- relevel(coldata$condition, ref = "control")

# a .gff annotation file extension may be gff2(gtf) or gff3 so check in files for use of = in the attribute field
# if '=' present it is gff3 if not it is gtf.
# see https://www.ensembl.org/info/website/upload/gff.html
# and http://gmod.org/wiki/GFF2#Converting_GFF2_to_GFF3
cat("Checking annotation file type.\n")
lines <- readLines(file(ref_annotation), n=10000)
# If transcript_id containing '=' (format eg. transcript_id=xxx)
# annotation type is gff3
check_file_type <- sum(grepl("transcript_id=", lines))
if (check_file_type != 0){
    cat("Annotation file type is gff3.\n")
    annotation_type <- "gff3"
} else {
    # otherwise gtf
    cat("Annotation file type is gtf.\n")
    annotation_type <- "gtf"
}

# Transcript_id versions (eg. ENTXXX.1, eg. ENTXXX.2) represent how many times that transcript reference has been changed 
# during its time in the database.
# Not all annotation files include it as part of the transcript_id - notably Ensembl
# The following handles this.
cat("Checking annotation file for presence of transcript_id versions.\n")
# Get the first transcript_id from the annotation file by parsing
lines <- readLines(file(ref_annotation), n=100000)
# Find transcript_ids in first 1000 lines and check if they contain dot (format eg. ENTXXX.1)
check_version <- sum(grepl("transcript_id[^;]+\\.", lines))
if (check_version != 0){
        # we do not need to strip the count file rows if ref_annotation includes versions
        cat("Annotation file transcript_ids include versions.\n")
    } else {
       # otherwise remove the versions
        rownames(cts) <- lapply(rownames(cts),  sub, pattern = "\\.\\d+$", replacement = "")
        cat("Annotation file transcript_ids do not include versions so also strip versions from the counts df.\n")
    }

cat("Loading annotation database.\n")
txdb <- makeTxDbFromGFF(ref_annotation,  format = annotation_type)
txdf <- select(txdb, keys(txdb,"GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx<- tab[match(txdf$GENEID, names(tab))]


cts <- cts[rownames(cts) %in% txdf$TXNAME, ] # FIXME: filter for transcripts which are in the annotation. Why they are not all there? 

# Reorder transcript/gene database to match input counts:
txdf <- txdf[match(rownames(cts), txdf$TXNAME), ]
rownames(txdf) <- NULL

# Create counts data frame:
counts<-data.frame(gene_id=txdf$GENEID, feature_id=txdf$TXNAME, cts)

# output unfiltered version of the counts table now we have paired transcripts with gene ids
write.table(counts, file="de_analysis/unfiltered_transcript_counts_with_genes.tsv", sep="\t", row.names = FALSE, quote=FALSE)

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
write.table(trs_cts, file="merged/filtered_transcript_counts_with_genes.tsv", sep="\t", row.names = FALSE, quote=FALSE)

gene_cts <- trs_cts_unfiltered %>% dplyr::select(c(1, 3:ncol(trs_cts)))  %>% group_by(gene_id) %>% summarise_all(tibble::lst(sum)) %>% data.frame()
rownames(gene_cts) <- gene_cts$gene_id
gene_cts$gene_id <- NULL
write.table(gene_cts, file="merged/all_gene_counts.tsv", sep="\t", quote=FALSE)

# Output count per million of the gene counts using edgeR CPM
cpm_gene_counts <- cpm(gene_cts)
# Add gene_id as index column header
cpm_gene_counts <- cbind(var_name = rownames(cpm_gene_counts), cpm_gene_counts)
rownames(cpm_gene_counts) <- NULL
colnames(cpm_gene_counts)[1] <- "gene_id"
write.table(cpm_gene_counts, file="de_analysis/cpm_gene_counts.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Differential gene expression using edgeR:
cat("Running differential gene expression analysis using edgeR.\n")

y <- DGEList(gene_cts)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit)
edger_res <- topTags(qlf, n=nrow(y), sort.by="PValue")[[1]]

pdf("de_analysis/results_dge.pdf")

# create status vector
status <- ifelse(
  qlf$PValue<0.01 & qlf$logFC>0, 
  'up', 
  ifelse(
    qlf$PValue<0.01 & qlf$logFC<=0,
    'down',
    'notsig'
  )
)
plotMD(qlf, status=status,  values=c("up","down","notsig"), hl.col=c("red","blue","black"))
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
