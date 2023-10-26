# 6556
# preliminary qc of cdna seq data from nanopore
# collected from output of adapted wf-transcriptomes
# read mapping stats
# exploration of salmon quantified cdna seq data


## ---- local_tests
#proj.dir="/Users/agata.smialowska/NBISproj/6556_isoforms_nanopore/results/pilot_30v2023/"
#proj.name.pref = "nanopore_pilot_30v2023"
#sample.info = "/Users/agata.smialowska/NBISproj/6556_isoforms_nanopore/scripts/qc_report/sample_info.txt"
#xenofilt = "TRUE" ## or FALSE

resdir=file.path(proj.dir,"results")

## ---- libraries

library(AnnotationHub)
library(ensembldb)

library(tidyverse)
library(dplyr)
library(tidyr)
library(kableExtra)

#library(pheatmap)
library(ggplot2)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)

library(DESeq2)
library(tximport)
#library(tximeta)
# library(edgeR)

#library(htmltools)
library(matrixStats)
library(org.Hs.eg.db)

is.xenofilt=isTRUE(xenofilt=="TRUE")


## ---- annot

ah <- AnnotationHub(localHub = FALSE)
ah_107 <- query(ah, "EnsDb.Hsapiens.v107")
edb_107 <- ah[[names(ah_107)]]


# loading from cache but just in case

#fname=file.path(datadir, "edb_107.rds")
#saveRDS(edb_107, file = fname)
#edb_107=readRDS(file = fname)


tx2gene_107 <- transcripts(edb_107, columns = c("tx_id_version", "gene_id", "tx_id"), return.type = "data.frame")
colnames(tx2gene_107) <- c("TXNAME", "GENEID", "tx_id")

tx2gene.2=tx2gene_107
tx2gene.2$TXNAME=tx2gene.2$tx_id


## ---- metadata

smpls_pths=read.table(sample.info ,sep="\t", header=FALSE, quote = "\"", dec = ".", fill = TRUE, row.names=NULL,blank.lines.skip=TRUE)
colnames(smpls_pths)=c("sample","smpl_rootdir")

n_smpls=nrow(smpls_pths)

samples=smpls_pths$smpl_rootdir
names(samples)=smpls_pths$sample

## ---- outdirs

#in the launcher
#wrk.dir=file.path(projdir,paste(proj.name.prefix,"report",sep="."))
reportdir=file.path(proj.dir,paste(proj.name.pref,"report","plots",sep="."))
dir.create(reportdir, recursive = TRUE)



## ---- read-stats-data


rlen_list=list()
xlim_readlen=2000

figurecapreadlen=paste("Histograms of read length distribution after preprocessing. X axis limited to", xlim_readlen,"to visualise most of data. Please note the Y axes range may differ between plots.",sep=" ")


# 
for (i in samples ){
    rootdir.i=i
    name.i=names(samples)[samples==i]

    print(rootdir.i)
    print(name.i)

    readstats.i=file.path(rootdir.i,"fastq_ingress_results",name.i,"fastcat_stats","per-read-stats.tsv")

    file_stats.i=read.table(readstats.i,sep="\t", header=TRUE, quote = "\"", dec = ".", fill = TRUE, row.names=NULL,blank.lines.skip=TRUE)


    hist_rlen=ggplot(file_stats.i, aes(x=read_length)) + 
      geom_histogram(binwidth=1, color="grey40", fill="grey70") +
      theme_bw(base_size = 18) + 
      theme(panel.grid.minor = element_blank() )+
      ggtitle(name.i) + 
      theme(aspect.ratio = 1/1.618) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    hist_rlen2=hist_rlen + 
     coord_cartesian(xlim=c(0,xlim_readlen) )

   rlen_list[[name.i]]=hist_rlen2

}

outfile_readlen=file.path(reportdir,"read_length_histograms.pdf")

  # cowplot::plot_grid(plotlist=rlen_list,
  #   ncol=2, labels=LETTERS[1:n_smpls])



pdf(outfile_readlen)

  cowplot::plot_grid(plotlist=rlen_list,
    ncol=2, labels=LETTERS[1:n_smpls])

dev.off()

## ---- summary-host-filt

#only in XenofilteR processed runs

xenofilt_stats=data.frame(sample=as.character(),graft=as.numeric())

xenofilt_stats_list=list()

if(is.xenofilt){

  for (i in samples ){
    rootdir.i=i
    name.i=names(samples)[samples==i]

    print(rootdir.i)
    print(name.i)

    host_filt_stats.i=file.path(rootdir.i,"fastq_filtered_graft",paste0(name.i,".host_filtering_stats.txt"))

    graft_reads=read.table(host_filt_stats.i,sep=" ", header=FALSE, quote = "\"", dec = ".", fill = TRUE, row.names=NULL,blank.lines.skip=TRUE)

    graft_read_cnt=graft_reads[1,1]/4


    x_stats=c(name.i, graft_read_cnt)

    xenofilt_stats=rbind(xenofilt_stats,x_stats)

    xenofilt_stats_list[[name.i]]=graft_read_cnt

  }

  colnames(xenofilt_stats)=c("name","graft_reads")
  # this will be later added to the table in the next chunk


}


## ---- read-map-stats-data

#if no host read filtering employed

if(!is.xenofilt){

  rlen_map=data.frame(sample=as.character(),reads_preproc=as.numeric(),mapped=as.numeric(),mapped_filt=as.numeric(),unmapped=as.numeric(),fraction_mapped=as.numeric(),fraction_mapped_filt=as.numeric(),
    avg_len_mapped=as.numeric(),max_len_mapped=as.numeric(),avg_len_unmapped=as.numeric(),max_len_unmapped=as.numeric() )

  for (i in samples ){
      rootdir.i=i
      name.i=names(samples)[samples==i]

      print(rootdir.i)
      print(name.i)

      read_map_stats.i=file.path(rootdir.i,"bam_minimap_genome_all",paste0(name.i,"_all_alns.minimap2.mapped.txt"))
      read_unmap_stats.i=file.path(rootdir.i,"bam_minimap_genome_all",paste0(name.i,"_all_alns.minimap2.unmapped.txt"))
      read_aln_stats.i=file.path(rootdir.i,"bam_minimap_genome_all",paste0(name.i,"_all_alns.minimap2.mapstats.txt"))

      pychopper_stats.i=file.path(rootdir.i,"fastq_pychopper",paste0(name.i,"_pychopper.tsv"))
      fastq_stats=read.table(pychopper_stats.i,sep="\t", header=FALSE, quote = "\"", dec = ".", fill = TRUE, row.names=NULL,blank.lines.skip=TRUE)
      reads_preproc=as.numeric(fastq_stats$V3[fastq_stats$V2=="Primers_found"])

      mapeed_filt.i=file.path(rootdir.i,"bam_minimap_genome_mapped",paste0(name.i,"_read_aln_stats.tsv"))
      mapped_filt_wftrx=read.table(mapeed_filt.i,sep="\t", header=TRUE, quote = "\"", dec = ".", fill = TRUE, row.names=NULL,blank.lines.skip=TRUE)
      reads_mapped_filt=mapped_filt_wftrx$PrimAln[1]


      read_map_stats=read.table(read_map_stats.i,sep="\t", header=FALSE, quote = "\"", dec = ".", fill = TRUE, row.names=NULL,blank.lines.skip=TRUE)

      reads_mapped=as.numeric(read_map_stats$V3[read_map_stats$V2=="reads mapped:"])
      total_bases_mapped=read_map_stats$V3[read_map_stats$V2=="total length:"]
      avg_len_mapped=read_map_stats$V3[read_map_stats$V2=="average length:"]
      max_len_mapped=read_map_stats$V3[read_map_stats$V2=="maximum length:"]

      read_unmap_stats=read.table(read_unmap_stats.i,sep="\t", header=FALSE, quote = "\"", dec = ".", fill = TRUE, row.names=NULL,blank.lines.skip=TRUE)
      unmapped_reads=as.numeric(read_unmap_stats$V3[read_unmap_stats$V2=="reads unmapped:"])
      total_bases_unmapped=read_unmap_stats$V3[read_unmap_stats$V2=="total length:"]
      avg_len_unmapped=read_unmap_stats$V3[read_unmap_stats$V2=="average length:"]
      max_len_unmapped=read_unmap_stats$V3[read_unmap_stats$V2=="maximum length:"]

      fraction_mapped=reads_mapped/reads_preproc
      fraction_mapped_filt=reads_mapped_filt/reads_preproc

      mapstats.i=c(name.i, reads_preproc, reads_mapped,reads_mapped_filt, unmapped_reads,fraction_mapped,fraction_mapped_filt,avg_len_mapped,max_len_mapped,avg_len_unmapped,max_len_unmapped)

      rlen_map=rbind(rlen_map,mapstats.i)

  }

  colnames(rlen_map)=c("name", "reads_preproc", "reads_mapped", "reads_mapped_filtered", "reads_unmapped", "fraction_mapped","fraction_mapped_filt","avg_len_mapped","max_len_mapped","avg_len_unmapped","max_len_unmapped")
  tab.title.1="Summary of read mapping statistics."
  tab.title.2="Read length of mapped and unmapped reads."


  cols.num=c("reads_preproc", "reads_mapped", "reads_mapped_filtered", "reads_unmapped", "fraction_mapped","fraction_mapped_filt","avg_len_mapped","max_len_mapped","avg_len_unmapped","max_len_unmapped")
  rlen_map[cols.num]=sapply(rlen_map[cols.num],as.numeric)

  cols.bigmark=c("reads_preproc", "reads_mapped", "reads_mapped_filtered", "reads_unmapped","avg_len_mapped","max_len_mapped","avg_len_unmapped","max_len_unmapped")
  rlen_map[cols.bigmark]=sapply(rlen_map[cols.bigmark], \(x) format(x,big.mark = ","))

  cols.frac=c("fraction_mapped","fraction_mapped_filt")
  rlen_map[cols.frac]=sapply(rlen_map[cols.frac], \(x) format(x,digits=3))

  cols_tab1=c("name", "reads_preproc", "reads_mapped", "reads_mapped_filtered", "reads_unmapped", "fraction_mapped","fraction_mapped_filt")
  cols_tab2=c("name","avg_len_mapped","max_len_mapped","avg_len_unmapped","max_len_unmapped")


  newnames.1=c("sample name", "reads preprocessed", "reads mapped", "reads mapped filtered", "reads unmapped", "fraction mapped","fraction mapped filtered")
  newnames.2=c("sample name","average length, mapped reads","max length, mapped reads","average length, unmapped reads","max length, unmapped reads")

  # rlen_map[cols_tab1] %>% 
  #   kable(booktabs = TRUE, row.names = FALSE, caption = tab.title.1, col.names=newnames.1) %>% 
  #   kable_minimal(full_width = TRUE) %>% 
  #   footnote(general = "Alignments were filtered using MAPQ cutoff 40, where indicated; otherwise no MAPQ filter was applied.")



  # rlen_map[cols_tab2] %>% 
  #   kable(booktabs = TRUE, row.names = FALSE, caption = tab.title.2, col.names=newnames.2) %>% 
  #   kable_minimal(full_width = TRUE)

}


if(is.xenofilt){  ## ---- read-map-stats-xenofilt


  rlen_map=data.frame(sample=as.character(),reads_preproc=as.numeric(),
    reads_graft=as.numeric(),fraction_graft_filt=as.numeric(),
    mapped=as.numeric(),mapped_filt=as.numeric(),unmapped=as.numeric(),fraction_mapped=as.numeric(),fraction_mapped_filt=as.numeric(),
    avg_len_mapped=as.numeric(),max_len_mapped=as.numeric(),avg_len_unmapped=as.numeric(),max_len_unmapped=as.numeric() )

  for (i in samples ){
      rootdir.i=i
      name.i=names(samples)[samples==i]

      print(rootdir.i)
      print(name.i)

      read_map_stats.i=file.path(rootdir.i,"bam_minimap_genome_all",paste0(name.i,"_all_alns.minimap2.mapped.txt"))
      read_unmap_stats.i=file.path(rootdir.i,"bam_minimap_genome_all",paste0(name.i,"_all_alns.minimap2.unmapped.txt"))
      read_aln_stats.i=file.path(rootdir.i,"bam_minimap_genome_all",paste0(name.i,"_all_alns.minimap2.mapstats.txt"))

      pychopper_stats.i=file.path(rootdir.i,"fastq_pychopper",paste0(name.i,"_pychopper.tsv"))
      fastq_stats=read.table(pychopper_stats.i,sep="\t", header=FALSE, quote = "\"", dec = ".", fill = TRUE, row.names=NULL,blank.lines.skip=TRUE)
      reads_preproc=as.numeric(fastq_stats$V3[fastq_stats$V2=="Primers_found"])

      mapeed_filt.i=file.path(rootdir.i,"bam_minimap_genome_mapped",paste0(name.i,"_read_aln_stats.tsv"))
      mapped_filt_wftrx=read.table(mapeed_filt.i,sep="\t", header=TRUE, quote = "\"", dec = ".", fill = TRUE, row.names=NULL,blank.lines.skip=TRUE)
      reads_mapped_filt=mapped_filt_wftrx$PrimAln[1]


      read_map_stats=read.table(read_map_stats.i,sep="\t", header=FALSE, quote = "\"", dec = ".", fill = TRUE, row.names=NULL,blank.lines.skip=TRUE)

      reads_mapped=as.numeric(read_map_stats$V3[read_map_stats$V2=="reads mapped:"])
      total_bases_mapped=read_map_stats$V3[read_map_stats$V2=="total length:"]
      avg_len_mapped=read_map_stats$V3[read_map_stats$V2=="average length:"]
      max_len_mapped=read_map_stats$V3[read_map_stats$V2=="maximum length:"]

      read_unmap_stats=read.table(read_unmap_stats.i,sep="\t", header=FALSE, quote = "\"", dec = ".", fill = TRUE, row.names=NULL,blank.lines.skip=TRUE)
      unmapped_reads=as.numeric(read_unmap_stats$V3[read_unmap_stats$V2=="reads unmapped:"])
      total_bases_unmapped=read_unmap_stats$V3[read_unmap_stats$V2=="total length:"]
      avg_len_unmapped=read_unmap_stats$V3[read_unmap_stats$V2=="average length:"]
      max_len_unmapped=read_unmap_stats$V3[read_unmap_stats$V2=="maximum length:"]

      reads_graft=xenofilt_stats_list[[name.i]]

      fraction_graft_filt=reads_graft/reads_preproc

      fraction_mapped=reads_mapped/reads_preproc
      fraction_mapped_filt=reads_mapped_filt/reads_preproc

      mapstats.i=c(name.i, reads_preproc, reads_graft, fraction_graft_filt,  reads_mapped,reads_mapped_filt, unmapped_reads,fraction_mapped,fraction_mapped_filt,avg_len_mapped,max_len_mapped,avg_len_unmapped,max_len_unmapped)

      rlen_map=rbind(rlen_map,mapstats.i)

  }

  colnames(rlen_map)=c("name", "reads_preproc", "reads_graft", "fraction_graft_filt","reads_mapped", "reads_mapped_filtered", "reads_unmapped", "fraction_mapped","fraction_mapped_filt","avg_len_mapped","max_len_mapped","avg_len_unmapped","max_len_unmapped")
  tab.title.1="Summary of read mapping statistics."
  tab.title.2="Read length of mapped and unmapped reads."


  cols.num=c("reads_preproc","reads_graft", "fraction_graft_filt", "reads_mapped", "reads_mapped_filtered", "reads_unmapped", "fraction_mapped","fraction_mapped_filt","avg_len_mapped","max_len_mapped","avg_len_unmapped","max_len_unmapped")
  rlen_map[cols.num]=sapply(rlen_map[cols.num],as.numeric)

  cols.bigmark=c("reads_preproc","reads_graft", "reads_mapped", "reads_mapped_filtered", "reads_unmapped","avg_len_mapped","max_len_mapped","avg_len_unmapped","max_len_unmapped")
  rlen_map[cols.bigmark]=sapply(rlen_map[cols.bigmark], \(x) format(x,big.mark = ","))

  cols.frac=c("fraction_graft_filt","fraction_mapped","fraction_mapped_filt")
  rlen_map[cols.frac]=sapply(rlen_map[cols.frac], \(x) format(x,digits=3))

  cols_tab1=c("name", "reads_preproc", "reads_graft", "fraction_graft_filt","reads_mapped", "reads_mapped_filtered", "reads_unmapped", "fraction_mapped","fraction_mapped_filt")
  cols_tab2=c("name","avg_len_mapped","max_len_mapped","avg_len_unmapped","max_len_unmapped")


  newnames.1=c("sample name", "reads preprocessed", "reads graft" , "fraction reads graft" , "reads mapped", "reads mapped filtered", "reads unmapped", "fraction mapped","fraction mapped filtered")
  newnames.2=c("sample name","average length, mapped reads","max length, mapped reads","average length, unmapped reads","max length, unmapped reads")

  # rlen_map[cols_tab1] %>% 
  #   kable(booktabs = TRUE, row.names = FALSE, caption = tab.title.1, col.names=newnames.1) %>% 
  #   kable_minimal(full_width = TRUE) %>% 
  #   footnote(general = "Alignments were filtered using MAPQ cutoff 40, where indicated; otherwise no MAPQ filter was applied.")%>% 
  #   footnote(general = "Filtering host reads also removed reads unmapped to neither graft nor host references.")


  # rlen_map[cols_tab2] %>% 
  #   kable(booktabs = TRUE, row.names = FALSE, caption = tab.title.2, col.names=newnames.2) %>% 
  #   kable_minimal(full_width = TRUE)

}



## ---- jaffal-data

fusions_list=list()
xlim_fus_spanning=50

fusions_stats=data.frame(sample=as.character(),n=as.numeric(),n_CO5=as.numeric(),
  max=as.numeric(),mean=as.numeric(),mean_CO5=as.numeric(),
  median=as.numeric(),median_CO5=as.numeric() )

top_fusion_genes=data.frame(sample=as.character(),top_fusions=as.character())

fus.all=data.frame()


for (i in samples ){
    rootdir.i=i
    name.i=names(samples)[samples==i]

    print(rootdir.i)
    print(name.i)

    jaffal_out.i=file.path(rootdir.i,paste0("jaffal_output_",name.i),paste0(name.i,"_jaffa_results.csv"))

    fusions=read.csv(jaffal_out.i, header=TRUE)


    # histogram
    hist_fus=ggplot(fusions, aes(x=spanning.reads)) + geom_histogram(binwidth=1, color="grey40", fill="grey70") +
      theme_bw(base_size = 18) + 
      theme(panel.grid.minor = element_blank() )+
      ggtitle(name.i) + 
      theme(aspect.ratio = 1/1.618) + 
      coord_cartesian(xlim=c(0,xlim_fus_spanning) )

    # tables
    fusions_list[[name.i]]=hist_fus

    # count fusions with at least 5 spanning reads
    # fusions.5=fusions[fusions$spanning.reads>=5,]

    # fusions_stats.i=c(name.i,nrow(fusions),nrow(fusions.5),
    #     max(fusions$spanning.reads),mean(fusions$spanning.reads),mean(fusions.5$spanning.reads),
    #     median(fusions$spanning.reads), median(fusions.5$spanning.reads))



    #count HighConfidence fusions
    fusions.HC=fusions[fusions$classification=="HighConfidence",]

    fusions_stats.i=c(name.i,nrow(fusions),nrow(fusions.HC),
        max(fusions$spanning.reads),mean(fusions$spanning.reads),mean(fusions.HC$spanning.reads),
        median(fusions$spanning.reads), median(fusions.HC$spanning.reads))


    fusions_stats=rbind(fusions_stats,fusions_stats.i)

    top_fusions=fusions %>% slice_max(fusions$spanning.reads,n=5)

    top_fusion_genes=rbind(top_fusion_genes,c(name.i,paste(top_fusions$fusion.genes, collapse=", ")))

    # for boxplot
    fus.i=cbind(fusions$spanning.reads,name.i)
    colnames(fus.i)=c("spanning.reads","sample")
    fus.all=rbind(fus.all,fus.i)

}

#colnames(fusions_stats)=c("name", "n_fusions", "n_fusions_CO5", "max_spanning_reads", "mean_spanning_reads", "mean_spanning_reads_CO5","median_spanning_reads", "median_spanning_reads_CO5")
colnames(fusions_stats)=c("name", "n_fusions", "n_fusions_HighConfidence", "max_spanning_reads", "mean_spanning_reads", "mean_spanning_reads_HighConfidence","median_spanning_reads", "median_spanning_reads_HighConfidence")

tab.title.3="Summary of detected transcript fusions."

#cols.num=c("n_fusions", "n_fusions_CO5", "max_spanning_reads", "mean_spanning_reads", "mean_spanning_reads_CO5", "median_spanning_reads", "median_spanning_reads_CO5")
cols.num=c("n_fusions", "n_fusions_HighConfidence", "max_spanning_reads", "mean_spanning_reads", "mean_spanning_reads_HighConfidence", "median_spanning_reads", "median_spanning_reads_HighConfidence")

fusions_stats[cols.num]=sapply(fusions_stats[cols.num],as.numeric)

#cols.frac=c("mean_spanning_reads", "mean_spanning_reads_CO5", "median_spanning_reads", "median_spanning_reads_CO5")
cols.frac=c("mean_spanning_reads", "mean_spanning_reads_HighConfidence", "median_spanning_reads", "median_spanning_reads_HighConfidence")

fusions_stats[cols.frac]=sapply(fusions_stats[cols.frac], \(x) format(x,digits=3))

newnames.3=c("name", "n fusions", "n fusions, HighConfidence", "max spanning reads", "mean spanning reads", "mean spanning reads, HighConfidence", "median spanning reads", "median spanning reads, HighConfidence")

# fusions_stats %>% 
#   kable(booktabs = TRUE, row.names = FALSE, caption = tab.title.3, col.names=newnames.3) %>% 
#   kable_minimal(full_width = TRUE) %>% 
#   footnote(general = "CO 5 - fusions detected with at least 5 spanning reads.")


tab.title.4="Top 5 transcript fusions by number of spanning reads."

# top_fusions_allsmpl=stack(top_fusion_genes)

# top_fusions_allsmpl=enframe(top_fusion_genes) %>% # creates the 'value' as a `list` column
#    mutate(value = map(value, as.character)) %>% # change to single type
#    unnest(cols=c("value"))

colnames(top_fusion_genes)=c("sample","top_fusions")

newnames.4=c("sample","top transcript fusions")

# top_fusion_genes %>% 
#   kable(booktabs = TRUE, row.names = FALSE, caption = tab.title.4, col.names=newnames.4) %>% 
#   kable_minimal(full_width = TRUE)







  # cowplot::plot_grid(plotlist=fusions_list,
  #   ncol=2, labels=LETTERS[1:n_smpls])


outfile_fusions1=file.path(reportdir,"reads_spanning_fusions_histograms.pdf")
pdf(outfile_fusions1)

  cowplot::plot_grid(plotlist=fusions_list,
    ncol=2, labels=LETTERS[1:n_smpls])

dev.off()


figurecapfusions1=paste("Histograms of fusion spanning read counts. X axis limited to", xlim_fus_spanning,"to visualise most of data. Please note the Y axes range may differ between plots.",sep=" ")



fus_spanning_reads_boxplot =  ggplot(fus.all, aes(x = factor(sample), y = as.numeric(spanning.reads), fill = factor(sample))) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log10') +
  theme_bw(base_size = 18) + 
  theme(panel.grid.minor = element_blank(),axis.text.x = element_text(margin=margin(5, b = 10), angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_viridis_d(alpha=0.6)+
  ggtitle("Number of reads spanning fusions") + 
  xlab("sample") + ylab("spanning reads") + labs(fill = "Sample") +
  theme(aspect.ratio = 1/1.618)


outfile_fusions2=file.path(reportdir,"reads_spanning_fusions_boxplot.pdf")
pdf(outfile_fusions2)

  fus_spanning_reads_boxplot

dev.off()

figurecapfusions2="Boxplot of fusion spanning read counts. Y axis is log10-transformed."


## ---- salmon-data

salmon_paths=list()

for (i in samples ){
    rootdir.i=i
    name.i=names(samples)[samples==i]

    salmon_dir.i=file.path(rootdir.i,"salmon",paste0(name.i,"_quant_salmon"),"quant.sf")
    salmon_paths[[name.i]]=salmon_dir.i

}

salmon_paths_list=unlist(salmon_paths)

txi.salmon_all <- tximport(salmon_paths_list, type = "salmon", tx2gene = tx2gene.2)

## add info on n expr, avg counts etc





## ---- gene_ids

# get gene ids and descriptions for genes in the count table

gene_names=data.frame(mapIds(org.Hs.eg.db, keys=rownames(txi.salmon_all$counts), keytype="ENSEMBL", column="GENENAME", multiVals="first"))
colnames(gene_names)="gene_name"
gene_names$ensembl_gene_id=rownames(gene_names)

gene_ids=data.frame(mapIds(org.Hs.eg.db, keys=rownames(txi.salmon_all$counts), keytype="ENSEMBL", column="SYMBOL", multiVals="first"))
colnames(gene_ids)="symbol"
gene_ids$ensembl_gene_id=rownames(gene_ids)

gene_annot=left_join(gene_ids,gene_names,by="ensembl_gene_id")



## ---- deseq2

### need the second sample for this

sampleTable_all = data.frame(condition = factor(colnames(txi.salmon_all$counts)))
rownames(sampleTable_all) = colnames(txi.salmon_all$counts)


dds_all = DESeqDataSetFromTximport(txi.salmon_all,sampleTable_all, ~1)
dds_all = DESeq(dds_all)
dds_vnorm=vst(dds_all,blind = TRUE, fitType = "parametric")


pca_data=plotPCA(dds_vnorm, intgroup = "condition", ntop = 500, returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot=ggplot(pca_data, aes(PC1, PC2, color=name)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw(base_size = 18) + 
  theme(panel.grid.minor = element_blank() )+
  scale_color_viridis_d(alpha=0.6)+
  coord_fixed(ratio = 1) + theme(aspect.ratio = 1) 

figurecapPCA="PCA plot. VST-transformed count estimates scaled to library size were used for calculations."




fpm_all=as.data.frame(fpm(dds_all))

# fpm.3r=fpm(dds, robust=TRUE) ## the same as above
fpm.all.long=fpm_all %>% 
  gather(smpls_all)

colnames(fpm.all.long)=c("sample","fpm")

figurecapTPM="Scaled counts per million summarised at gene level. FPM - fragments / reads per million."

#boxplot_gene_counts=ggplot(fpm.all.long, aes(x = factor(sample), y = as.numeric(fpm), fill = factor(sample))) + 

boxplot_gene_counts=ggplot(fpm.all.long, aes(x =sample, y =fpm, fill = sample)) +
  geom_boxplot() + 
  scale_y_continuous(trans='log10') +
  theme_bw(base_size = 18) + 
  theme(panel.grid.minor = element_blank(),axis.text.x = element_text(margin=margin(5, b = 10), angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_viridis_d(alpha=0.6)+
  ggtitle("Scaled counts per million") + 
  theme(aspect.ratio = 1/1.618)



## ---- heatmaps-data

ntop=60

# highest expression
select_exprs <- order(rowMeans(counts(dds_all,normalized=TRUE)),decreasing=TRUE)[1:ntop]
#standard deviation
select_sd <- order(apply(as.data.frame(counts(dds_all,normalized=TRUE)), 1, sd, na.rm=TRUE),decreasing=TRUE)[1:ntop]
#variance in data
select_var <- order(rowVars(counts(dds_all,normalized=TRUE)),decreasing=TRUE)[1:ntop]

dat_hm_exprs=assay(dds_vnorm)[select_exprs,]
dat_hm_var=assay(dds_vnorm)[select_var,]


# change rownames to symbol

g_annot=gene_annot
g_annot=g_annot %>%
 mutate(symbol = coalesce(symbol,ensembl_gene_id ))

g_names=g_annot$symbol
names(g_names)=g_annot$ensembl_gene_id

rownames(dat_hm_exprs)<-g_names[rownames(dat_hm_exprs)]
rownames(dat_hm_var)<-g_names[rownames(dat_hm_var)]




#for report
#heatmaps without scaling
hm1=pheatmap(dat_hm_exprs)
#with scaling - when more samples
#hm2=pheatmap(dat_hm_exprs, cluster_rows=TRUE, scale="row")


#highest standard deviation
#row_stdev <- apply(as.data.frame(counts(dds,normalized=TRUE)), 1, sd, na.rm=TRUE)
#hm3=pheatmap(assay(dds_vnorm)[select_sd,])
#hm4=pheatmap(assay(dds_vnorm)[select_sd,], cluster_rows=TRUE,show_rownames=FALSE, scale="row")

# highest variance
hm5=pheatmap(dat_hm_var)
#hm6=pheatmap(assay(dds_vnorm)[select_var,], cluster_rows=TRUE,show_rownames=FALSE, scale="row")

#figurecap3="Heatmap of rlog-transformed and library size-scaled log-counts for top 60 genes, selected by expression."
#figurecap4="Heatmap of rlog-transformed and library size-scaled log-counts for top 60 genes, selected by variance."

figurecap_hm1="Heatmaps of vst-transformed and library size-scaled log-counts for top 60 genes, selected by expression. A - z-scored log-counts (scaled by row to show deviation from the row average); A - unscaled log-counts."
figurecap_hm2="Heatmaps of vst-transformed and library size-scaled log-counts for top 60 genes, selected by variance. A - z-scored log-counts (scaled by row to show deviation from the row average); A - unscaled log-counts."


hm2=pheatmap(dat_hm_exprs, cluster_rows=TRUE, column_title="B", heatmap_legend_param=list(title="Logcounts"),fontsize_row = 7, fontsize_col = 7,color=colorRampPalette(brewer.pal(6,name="PuOr"))(12))
hm1=pheatmap(dat_hm_exprs, cluster_rows=TRUE,  scale="row",column_title="A", heatmap_legend_param=list(title="Z score"),fontsize_row = 7,fontsize_col = 7,color=colorRampPalette(rev(brewer.pal(6,name="PRGn")))(12))

hms_exprs = hm1 + hm2

draw(hms_exprs, 
    column_title = "Highest expression", column_title_gp = gpar(fontsize = 16))

hm4=pheatmap(dat_hm_var, cluster_rows=TRUE, column_title="B", heatmap_legend_param=list(title="Logcounts"),fontsize_row = 7, fontsize_col = 7,color=colorRampPalette(brewer.pal(6,name="PuOr"))(12))
hm3=pheatmap(dat_hm_var, cluster_rows=TRUE,  scale="row",column_title="A", heatmap_legend_param=list(title="Z score"),fontsize_row = 7,fontsize_col = 7,color=colorRampPalette(rev(brewer.pal(6,name="PRGn")))(12))

hms_var = hm3 + hm4



