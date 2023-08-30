#!/usr/bin/env Rscript

rm(list=ls())

library(knitr)
args <- commandArgs(TRUE)

if (length(args) < 4) stop("Not all args are set; required: projdir proj.name.prefix sample.info xenofilt")

proj.dir <- args[1]
proj.name.pref <- args[2]
sample.info <- args[3]
xenofilt <- args[4] # TRUE or FALSE

#user.run=Sys.getenv("USER")

#wrk.dir=file.path(proj.dir,"results",data.type,"report")

dir.create(proj.dir, recursive = TRUE)

wrk.dir=file.path(proj.dir,paste(proj.name.pref,"report",sep="."))
dir.create(wrk.dir, recursive = TRUE)

print(proj.dir)
print(proj.name.pref)
print(sample.info)
print(wrk.dir)

print(getwd())

rmarkdown::render('6556_QC_report_v0.1.Rmd', output_file = file.path(wrk.dir,paste("QC_report",proj.name.pref,'html', sep=".")))


