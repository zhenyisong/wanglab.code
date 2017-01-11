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
#library(DiagrammeR)
#library(magrittr)
#library(cowplot)
#library(clusterProfiler)
#library(png)
library(grid)
library(GGally)
library(reshape)

genome_ref.path      <- "/home/zhenyisong/data/bringback/igenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
setwd('/home/zhenyisong/data/wanglilab/projects/VZR20160421/Rawdata')
reads.files.names    <- list.files(pattern = '*.fastq')
read.path.1          <- reads.files.names[grep("R1",reads.files.names)]
read.path.2          <- reads.files.names[grep("R2",reads.files.names)]
reads.paths.1        <- paste0(getwd(),'/',read.path.1)
reads.paths.2        <- paste0(getwd(),'/',read.path.2)
setwd('/home/zhenyisong/data/wanglilab/projects/VZR20160421')
unlink('rsubread')
dir.create('rsubread')
output.path <- '/home/zhenyisong/data/wanglilab/projects/VZR20160421/rsubread'
setwd(output.path)

base.string       <-  'hg19_index'
buildindex( basename = base.string, reference = genome_ref.path )
align( index          = base.string, 
       readfile1      = reads.paths.1, 
       readfile2      = reads.paths.2, 
       input_format   = "FASTQ", 
       type           = 'rna',
       output_file    = outputs.files, 
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 4, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

liangping.genes      <- featureCounts( outputs.files, useMetaFeatures = TRUE,
                                       countMultiMappingReads = FALSE,
                                       strandSpecific         = 0, 
                                       isPairedEnd            = TRUE,
                                       requireBothEndsMapped  = TRUE,
                                       autosort               = TRUE,
                                       nthreads               = 8,
                                       annot.inbuilt = "hg19", allowMultiOverlap = TRUE)

setwd("/home/zhenyisong/data/cardiodata")
save.image(file = 'liangping.Rdata')
quit("no")