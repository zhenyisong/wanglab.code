library(gplots)
library(xlsx)
library(Rsubread)
library(edgeR)
library(limma)
library(DESeq2)
library(genefilter)
library(RColorBrewer)
library(gridExtra)
library(Rsamtools)
library(clusterProfiler)
library(pathview)
library(tidyverse)
library(biomaRt)
library(cluster)
library(factoextra)


# data source 
# download from ENCODE database
# https://www.encodeproject.org/
# using filter to select targeted data set
# Assay category Transcription
# Assay total RNA-seq
# Organism Homo Sapiens
# Biosample type tissues
# Organ heart
# Available data fastq
# Using batch download
# nohup xargs -n 1 curl -O -C - -L < download.txt &
#---


genome_ref.path      <- "/mnt/date/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
unlink('rsubread', force = TRUE, recursive = TRUE)
dir.create('rsubread')
output.path          <- '/home/zhenyisong/biodata/cardiodata/encodeHeart/humanHeart/rsubread'
setwd(output.path)
yaoyan.qc.rawdata    <- list.files( path = "../",pattern = "\\.fastq$", all.files = FALSE,
                                    full.names = TRUE, recursive = TRUE, 
                                    ignore.case = FALSE, include.dirs = TRUE)
base.string          <-  'hg19_index'
#buildindex( basename = base.string, reference = genome_ref.path )
reads.paths.1         <- yaoyan.qc.rawdata[c(6,4,2,8,11,12)]
reads.paths.2         <- yaoyan.qc.rawdata[c(5,3,1,7,9,10)]
outputs.files         <- c('GSE87943.bam', 'GSE78567.bam','GSE78567.2.bam','GSE88367.bam','GSE88247.bam','GSE88271.bam')
align( index          = base.string, 
       readfile1      = reads.paths.1, 
       readfile2      = reads.paths.2, 
       input_format   = "FASTQ", 
       type           = 'rna', 
       output_file    = outputs.files, 
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 25, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )


yaoyan.qc.mRNA      <- featureCounts( outputs.files, useMetaFeatures = TRUE,
                                      countMultiMappingReads = FALSE,
                                      strandSpecific         = 0, 
                                      isPairedEnd            = TRUE,
                                      requireBothEndsMapped  = TRUE,
                                      autosort               = TRUE,
                                      nthreads               = 25,
                                      annot.inbuilt = "hg19", allowMultiOverlap = TRUE)


gene         <- yaoyan.qc.mRNA
gene.counts  <- gene$counts
gene.ids     <- gene$annotation$GeneID

columns      <- c("ENTREZID","SYMBOL", "GENENAME");
GeneInfo     <- AnnotationDbi::select( org.Hs.eg.db, keys= as.character(gene.ids), 
                       keytype = "ENTREZID", columns = columns);
m            <- match(gene$annotation$GeneID, GeneInfo$ENTREZID);
Ann          <- cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                              GeneInfo[m, c("SYMBOL", "GENENAME")]);
             
Ann$Chr      <- unlist( lapply(strsplit(Ann$Chr, ";"), 
                        function(x) paste(unique(x), collapse = "|")))
Ann$Chr      <- gsub("chr", "", Ann$Chr)

yaoyan.qc.gene.exprs.log  <- DGEList(counts = gene.counts, genes = Ann) %>%
                                     calcNormFactors()  %>% rpkm(log = T)

save.image(file = 'yaoyan.qc.Rdata')
quit('no')