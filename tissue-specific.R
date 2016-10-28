# see similiar program
# cardioFibroblast.R

library(gplots)
library(xlsx)
library(Rsubread)
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(DESeq2)
library(gplots)
library(genefilter)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(annotate)
library(limma)


setwd('/home/zhenyisong/data/cardiodata/SRP012040')
targets.file      <- '/home/zhenyisong/data/cardiodata/SRP012040/targets.txt'
reads.files       <- read.table(targets.file,header = F)
reads.path        <- '/home/zhenyisong/data/cardiodata/SRP012040/'
output.path       <- '/home/zhenyisong/data/cardiodata/SRP012040/results/'

reads.files.names <- reads.files$V1
read.path.1       <- reads.files.names[grep("_1",reads.files.names)]
read.path.2       <- reads.files.names[grep("_2",reads.files.names)]
reads.paths.1     <- paste0(reads.path,read.path.1)
reads.paths.2     <- paste0(reads.path,read.path.2)
outputs.files     <- paste0(output.path,read.path.1,'.bam')
base.string       <- 'mm10_index'
setwd('/home/zhenyisong/data/cardiodata/SRP051406/')
align( index         = base.string, 
       readfile1     = reads.paths.1, 
       readfile2     = reads.paths.2, 
       input_format  = "FASTQ", 
       type          = 'rna',
       output_file   = outputs.files, 
       output_format = "BAM", 
       nthreads      = 4, 
       indels        = 1,
       maxMismatches = 3,
       phredOffset   = 33,
       unique        = T )
tissue.list      <- featureCounts( outputs.files, useMetaFeatures = TRUE, 
                                   annot.inbuilt = "mm10", allowMultiOverlap = TRUE)

fibro.path       <- '/home/zhenyisong/data/cardiodata/SRP043191/results/'
fibro.path.files <- c('SRR1390714.fastq.sam','SRR1390715.fastq.sam','SRR1390716.fastq.sam')
outputs.files    <- paste0(fibro.path ,fibro.path.files)

fibro.list       <- featureCounts( outputs.files, useMetaFeatures = TRUE, 
                                   annot.inbuilt = "mm10", allowMultiOverlap = TRUE)
gene.counts      <- cbind(tissue.list$counts, fibro.list$counts)
gene.ids         <- tissue.list$annotation$GeneID
columns  <- c("ENTREZID","SYMBOL", "MGI", "GENENAME");
GeneInfo <- select( org.Mm.eg.db, keys= as.character(gene.ids), 
                   keytype="ENTREZID", columns = columns);
m        <- match(gene.ids, GeneInfo$ENTREZID);
Ann      <- cbind( tissue.list$annotation[, c("GeneID", "Chr","Length")],
                          GeneInfo[m, c("SYMBOL", "MGI", "GENENAME")]);
Ann$Chr  <-  unlist( lapply(strsplit(Ann$Chr, ";"), 
                    function(x) paste(unique(x), collapse = "|")))
Ann$Chr    <- gsub("chr", "", Ann$Chr)
gene.exprs <- DGEList(counts = gene.counts, genes = Ann)
gene.exprs <- calcNormFactors(gene.exprs)
dge.tmm    <- t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
dge.tmm    <- log(dge.tmm + 1)
rownames(dge.tmm) <- GeneInfo[m,'SYMBOL'];

#------------------------------------------------------------------------
# PMID: 15388519 
# Bioinformatics. 2005 Mar 1;21(5):650-9. Epub 2004 Sep 23.
# Genome-wide midrange transcription profiles reveal expression 
# level relationships in human tissue specification.
#------------------------------------------------------------------------
tao.func <- function(x) {
    max.val  <- max(x)
    size     <- length(x)
    norm.val <- x/max.val
    norm.val <- 1 - norm.val
    norm.val <- sum(norm.val)
    norm.val <- norm.val/(size - 1)
    return(norm.val)
}

tao.index <- apply(dge.tmm, 1, tao.func)
setwd('/home/zhenyisong/data/cardiodata')
save.image('tissue-specific.Rdata')
q("no")