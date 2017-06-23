# @author Yisong Zhen
# @ rawdata
# contact with Dr. Liang, Ping, his paired-end project
# the preprocessing was completed by yaofang
# nohup R CMD BATCH /home/zhenyisong/biodata/wanglab/wangcode/liangping.restart.R &

# the rawdata is not completed
# find ./ -type f -name "*LP*"
# @since 2017-06-20
#---

library(tidyverse)
library(Rsubread)
library(org.Hs.eg.db)
library(annotate)
library(edgeR)
library(limma)
library(DESeq2)

# data transfer 
# in window 10
#---

#x1.runing.path <- file.path('D:\\wangli_data\\Rdata')
#setwd(x1.runing.path)
#load('liangp.Rdata')
#

#--

human.genome.path      <- file.path('/mnt/date/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa')
liangp.output.path     <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/lianglab/nextdata')

human.rawdata.path.1   <- file.path('/mnt/date/Sequencing/FastQ/20170605_LF_RNAseq_LPBD')
liangp.files.1         <- list.files( path = human.rawdata.path.1, pattern = "h\\_.*\\.gz$", all.files = FALSE,
                                      full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = F) %>%
                                      `[`(-c(13:24))

human.rawdata.path.2   <- file.path('/mnt/date/Sequencing/FastQ/170616_LF_RNAseq')
liangp.files.2         <- list.files( path = human.rawdata.path.2, pattern = "LP-2.*$", all.files = FALSE,
                                      full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = F)

human.rawdata.path.3   <- file.path('/mnt/date/Sequencing/FastQ/20170621_LF_RNAseq')
liangp.files.3        <- list.files( path = human.rawdata.path.3, pattern = "LP-1.*$", all.files = FALSE,
                                      full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = F)

liangp.files           <- c(liangp.files.1, liangp.files.2, liangp.files.3)
read.1.files           <- liangp.files[grep("R1",liangp.files)]
read.2.files           <- liangp.files[grep("R2",liangp.files)]

human.output.filenames <- basename(read.1.files) %>% sub(pattern = '_R1_001.fastq.gz', replacement = '') %>%
                          paste0(liangp.output.path,'/', . ,'.bam')

setwd(liangp.output.path)
base.string          <-  'hg38'
buildindex( basename = base.string, reference = human.genome.path)
align( index          = base.string, 
       readfile1      = read.1.files, 
       readfile2      = read.2.files, 
       input_format   = "gzFASTQ", 
       type           = 'rna',
       output_file    = human.output.filenames, 
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 15, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )
liangp.genes      <- featureCounts( human.output.filenames, useMetaFeatures = TRUE,
                                       countMultiMappingReads = FALSE,
                                       strandSpecific         = 0, 
                                       isPairedEnd            = TRUE,
                                       requireBothEndsMapped  = TRUE,
                                       autosort               = TRUE,
                                       nthreads               = 15,
                                       annot.inbuilt = "hg38", allowMultiOverlap = TRUE)


gene         <- liangp.genes
gene.counts  <- gene$counts
gene.ids     <- gene$annotation$GeneID

columns      <- c("ENTREZID","SYMBOL", "GENENAME");
GeneInfo     <- select( org.Hs.eg.db, keys= as.character(gene.ids), 
                        keytype = "ENTREZID", columns = columns);
m            <- match(gene$annotation$GeneID, GeneInfo$ENTREZID);
Ann          <- cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                       GeneInfo[m, c("SYMBOL", "GENENAME")]);

Ann$Chr      <- unlist( lapply(strsplit(Ann$Chr, ";"), 
                        function(x) paste(unique(x), collapse = "|")))
Ann$Chr      <- gsub("chr", "", Ann$Chr)

gene.exprs   <- DGEList(counts = gene.counts, genes = Ann) %>% calcNormFactors()

save.image('liangp.Rdata')
quit('no')