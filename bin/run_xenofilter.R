#!/usr/bin/env Rscript


args <- commandArgs(TRUE)

if (length(args) < 3) stop("Not all args are set; required: smpl_name bam_graft bam_host ")

smpl_name <- args[1]
bam_graft <- args[2]
bam_host <- args[3]


library(XenofilteR)


sample.list=as.data.frame(cbind(bam_graft,bam_host))
colnames(sample.list)=c("graft","host")

destination.folder="."

MM_threshold_ont=200

output.names=c(smpl_name)

bp.param=SnowParam(workers = 1, type = "SOCK")


XenofilteR(sample.list, destination.folder, bp.param=bp.param, output.names, MM_threshold=MM_threshold_ont)


