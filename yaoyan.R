# @data source
#     deposited at Fuwai SuperComputer
#     /home/zhenyisong/data/wanglilab/projects/yaoyan/rawdata
# @author Yisong Zhen
# @since 2017-03-08
# cd /home/zhenyisong/data/wanglilab/wangcode
# nohup R CMD BATCH yaoyan.R &
#---
library(gplots)
library(xlsx)
library(Rsubread)
library(edgeR)
library(limma)
library(DESeq2)
library(genefilter)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(affy)
library(annotate)
library(org.Hs.eg.db)
library(VennDiagram)
library(Rsamtools)


setwd('/wa/zhenyisong/wanglilab/projects/yaoyan/rawdata/rsubread')
load("yaoyan.Rdata")

genome_ref.path      <- "/home/zhenyisong/data/bringback/igenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
setwd('/wa/zhenyisong/wanglilab/projects/yaoyan/rawdata')
unlink('rsubread', force = TRUE)
dir.create('rsubread')
output.path          <- '/wa/zhenyisong/wanglilab/projects/yaoyan/rawdata/rsubread'
setwd('/wa/zhenyisong/wanglilab/projects/yaoyan/rawdata/rsubread')
yaoyan.rawdata       <- list.files( path = "../",pattern = "\\.fq.\\gz$", all.files = FALSE,
                                    full.names = TRUE, recursive = TRUE, ignore.case = FALSE, include.dirs = TRUE)
base.string          <-  'hg19_index'
buildindex( basename = base.string, reference = genome_ref.path )
mRNA.yaoyan.rawdata  <- yaoyan.rawdata[1:60]
reads.paths.1        <- mRNA.yaoyan.rawdata[grep("\\.1\\.clean\\.fq\\.gz",mRNA.yaoyan.rawdata)]
reads.paths.2        <- mRNA.yaoyan.rawdata[grep("\\.2\\.clean\\.fq\\.gz",mRNA.yaoyan.rawdata)]
read.file.names      <- basename(reads.paths.1)
read.file.names      <- sub("\\.\\d\\D+$", "\\1", read.file.names, perl = TRUE)
outputs.files        <- paste0(output.path,'/', read.file.names,'.bam')
align( index          = base.string, 
       readfile1      = reads.paths.1, 
       readfile2      = reads.paths.2, 
       input_format   = "gzFASTQ", 
       type           = 'rna', 
       output_file    = outputs.files, 
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 4, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

# merge 2 bam file into one
# file needed to sort
# 
file.needed.sort <- outputs.files[7:30]
file.odd.bam     <- file.needed.sort[1:length(file.needed.sort) %% 2 == 1]
file.even.bam    <- file.needed.sort[1:length(file.needed.sort) %% 2 == 0]
merged.bams      <- c()
for( i in c(1:length(file.odd.bam)) ) {
   odd.name  <- basename(file.odd.bam[i])
   even.name <- basename(file.even.bam[i])
   file.name <- sub("_lane.*$", "\\1", odd.name, perl = TRUE)
   file.odd.temp  <- paste0(file.name,"_temp_1")
   file.even.temp <- paste0(file.name,"_temp_2")
   bam1 <- sortBam(odd.name, file.odd.temp);
   bam2 <- sortBam(even.name, file.even.temp);
   file.name <- paste0(file.name,'.bam')
   bam3 <- mergeBam(c(bam1,bam2), file.name, overwrite = TRUE)
   merged.bams <- c(merged.bams,bam3)
}

unlink("*_temp_*\\.bam")


yaoyan.mRNA.genes      <- featureCounts( outputs.files, useMetaFeatures = TRUE,
                                         countMultiMappingReads = FALSE,
                                         strandSpecific         = 0, 
                                         isPairedEnd            = TRUE,
                                         requireBothEndsMapped  = TRUE,
                                         autosort               = TRUE,
                                         nthreads               = 4,
                                         annot.inbuilt = "hg19", allowMultiOverlap = TRUE)

save.image(file = "yaoyan.Rdata")
quit("no")




