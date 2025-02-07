#!/usr/bin/env Rscript

suppressMessages(library(argparser))

parser <- arg_parser("Plot results")
parser <- add_argument(parser, "--counts", help="Filtered transcript counts with genes.")
parser <- add_argument(parser, "--results_dtu", help="stageR results.")
parser <- add_argument(parser, "--sample_sheet", help="Sample sheet.")
parser <- add_argument(parser, "--pdf_out", help="PDF file name.")
argv <- parse_args(parser)

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))

# Set up sample data frame:
coldata <- read.csv(argv$sample_sheet, row.names="alias", sep=",")
coldata$condition <- factor(coldata$condition, levels=rev(levels(coldata$condition)))
coldata$type <-NULL
coldata$patient <-NULL

# Read stageR results:
stageR <- read.csv(argv$results_dtu, sep="\t")
names(stageR) <- c("gene_id", "transcript_id", "p_gene", "p_transcript");

# Read filtered counts:
counts <- read.csv(argv$counts, sep="\t");
names(counts)[2]<-"transcript_id"

# Join counts and stageR results:
df <- counts %>% left_join(stageR, by = c("gene_id", "transcript_id"))
df <- df[order(df$p_gene),]

scols <- setdiff(names(df),c("gene_id", "transcript_id", "p_gene", "p_transcript"))

# Normalise counts:
for(sc in scols){
    df[sc] <- df[sc] / sum(df[sc])
}

# Melt data frame:
tdf <- df %>% gather(key='sample', value='norm_count',-gene_id, -transcript_id, -p_gene, -p_transcript)

# Add sample group column:
sampleToGroup<-function(x){
    return(coldata[x,]$condition)
}

tdf$group <- sampleToGroup(tdf$sample)

# Filter for significant genes:
sig_level <- 0.05
genes <- as.character(tdf[which(tdf$p_gene < sig_level),]$gene_id)
genes <- unique(genes)

pdf(argv$pdf_out)

for(gene in genes){
    gdf<-tdf[which(tdf$gene_id==gene),]
    p_gene <- unique(gdf$p_gene)
    dtu_plot <- ggplot(gdf, aes(x=transcript_id, y=norm_count)) + geom_bar(stat="identity", aes(fill=sample), position="dodge")
    dtu_plot <- dtu_plot + facet_wrap(~ group) + coord_flip()
    dtu_plot <- dtu_plot + ggtitle(paste(gene," : p_value=",p_gene,sep=""))
    print(dtu_plot)
}
