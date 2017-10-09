library(Rsubread)
library(tidyverse)
library(org.Hs.eg.db)

##################################################################
###先进行LP_hh的数据处理
##################################################################


genome_ref.path   <- "/mnt/date/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
reads.path        <- '/mnt/date/Project/yaofang/20170605_LF_RNAseq_LPBD/hh/'
output.path       <- '/home/zhenyisong/biodata/wanglab/wangdata/yaofang/'

setwd(output.path)

rawdata.path      <- file.path(reads.path)
yaofang.files     <- list.files( path = reads.path, pattern = "\\.fastq.\\gz$", all.files = FALSE,
                                 full.names = TRUE, recursive = TRUE, 
                                 ignore.case = FALSE, include.dirs = TRUE)

read.path.1       <- yaofang.files[grep("_001",yaofang.files)]
read.path.2       <- yaofang.files[grep("_002",yaofang.files)]

base.string       <-  'hg19_index'
buildindex( basename = base.string, reference = genome_ref.path )

read.file.names      <- basename(read.path.1)
read.file.names      <- sub("\\.fastq\\.gz$", "", read.file.names, perl = TRUE)

outputs.files        <- paste0(output.path,'/', read.file.names,'.bam')
align( index          = base.string, 
       readfile1      = read.path.1, 
       readfile2      = read.path.2, 
       input_format   = "gzFASTQ", 
       type           = 'rna', 
       output_file    = outputs.files, 
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 15, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

yaofang.mRNA.genes      <- featureCounts( outputs.files, useMetaFeatures = TRUE,
                                         countMultiMappingReads = FALSE,
                                         strandSpecific         = 0, 
                                         isPairedEnd            = TRUE,
                                         requireBothEndsMapped  = TRUE,
                                         autosort               = TRUE,
                                         nthreads               = 25,
                                         annot.inbuilt = "hg19", allowMultiOverlap = TRUE)

