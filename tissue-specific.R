# see similiar program
# cardioFibroblast.R
# @parent
#      hisat2.R
# @raw data
#    GSE36025/SRP012040
#    PMID: 25582902
# I fetch the raw data from NCBI usng aspera
# aspera.script

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


#---
# raw data was processed by Hisat2
#---
setwd('/home/zhenyisong/data/cardiodata/SRP012040/hisat2')
file.name <- 'gene_count_matrix.csv'
SRP012040.hisat2.data <- read.csv( file.name, header = TRUE, sep = ",",
                         row.names = 1, stringsAsFactors = FALSE)

setwd('/home/zhenyisong/data/cardiodata/SRP047225/hisat2')
file.name         <- 'gene_count_matrix.csv'
SRP047225.hisat2.data <- read.csv( file.name, header = TRUE, sep = ",",
                                   row.names = 1, stringsAsFactors = FALSE)

hisat2.data <- merge( SRP012040.hisat2.data, SRP047225.hisat2.data, 
                      by = "row.names", all.x = FALSE)


#---
# this protocol see the limma manual
#   PDF version: First edition 2 December 2002
#   Last revised 1 March 2016
#   pp.117-120
#---

gene.counts     <- as.matrix(sapply(hisat2.data[,-1], as.integer))
gene.exprs      <- DGEList(counts = gene.counts)
gene.exprs      <- calcNormFactors(gene.exprs)
dge.tmm         <- t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
dge.tmm.counts  <- apply(dge.tmm,2, as.integer)
log.trans       <- log(dge.tmm.counts + 1)

adult.Ovary          <- apply(log.trans[,c(1:10)],1, median)
adlut.MammaryGland   <- apply(log.trans[,c(11:16)],1, median)
adult.Stomach        <- apply(log.trans[,c(17:22)],1, median)
adlut.SmIntestine    <- apply(log.trans[,c(23:32)],1, median)
adult.Duodenum       <- apply(log.trans[,c(33:39)],1, median)
adult.Adrenal        <- apply(log.trans[,c(40:45)],1, median)
adult.LgIntestine    <- apply(log.trans[,c(46:49)],1, median)
adult.GenitalFatPad  <- apply(log.trans[,c(50:53)],1, median)
adult.SubcFatPad     <- apply(log.trans[,c(54:57)],1, median)
adult.Thymus         <- apply(log.trans[,c(58:63)],1, median)
adult.Testis         <- apply(log.trans[,c(64:67)],1, median)
adult.Kidney         <- apply(log.trans[,c(68:73)],1, median)
adult.Liver          <- apply(log.trans[,c(74:79)],1, median)
adult.Lung           <- apply(log.trans[,c(80:83)],1, median)
adult.Spleen         <- apply(log.trans[,c(84:89)],1, median)
adult.Colon          <- apply(log.trans[,c(90:95)],1, median)
adult.Heart          <- apply(log.trans[,c(96:99)],1, median)
adult.FrontalLobe    <- apply(log.trans[,c(100:101)],1, median)
adult.Cortex         <- apply(log.trans[,c(102:103)],1, median)
adult.Bladder        <- apply(log.trans[,c(104:105)],1, median)
adult.Placenta       <- apply(log.trans[,c(106:107)],1, median)
fetal.Liver          <- apply(log.trans[,c(108:109)],1, median)
adult.Cerebellum     <- apply(log.trans[,c(110:111)],1, median)
fetal.Limb           <- apply(log.trans[,c(112:113)],1, median)
fetal.CNS.E14        <- apply(log.trans[,c(114:115)],1, median)
fetal.CNS.E18        <- apply(log.trans[,c(116:117)],1, median)
fetal.Liver.E14.5    <- apply(log.trans[,c(118:119)],1, median)
fetal.WholeBrain_E14.5    <- apply(log.trans[,c(120:121)],1, median)
fetal.CNS.E11.5           <- apply(log.trans[,c(122:123)],1, median)
fetal.Liver.E14           <- apply(log.trans[,c(124:125)],1, median)
adult.fibroblast          <- apply(log.trans[,c(126:128)],1, median)

tisse.matrix  <- matrix()
tissue.matrix <- cbind(adult.Ovary,adlut.MammaryGland)
tissue.matrix <- cbind(tissue.matrix,adult.Stomach)
tissue.matrix <- cbind(tissue.matrix,adlut.SmIntestine)
tissue.matrix <- cbind(tissue.matrix,adult.Duodenum)
tissue.matrix <- cbind(tissue.matrix,adult.Adrenal)
tissue.matrix <- cbind(tissue.matrix,adult.LgIntestine)
tissue.matrix <- cbind(tissue.matrix,adult.GenitalFatPad)
tissue.matrix <- cbind(tissue.matrix,adult.SubcFatPad)
tissue.matrix <- cbind(tissue.matrix,adult.Thymus)
tissue.matrix <- cbind(tissue.matrix,adult.Testis)
tissue.matrix <- cbind(tissue.matrix,adult.Kidney)
tissue.matrix <- cbind(tissue.matrix,adult.Liver)
tissue.matrix <- cbind(tissue.matrix,adult.Lung)
tissue.matrix <- cbind(tissue.matrix,adult.Colon)
tissue.matrix <- cbind(tissue.matrix,adult.Spleen)
tissue.matrix <- cbind(tissue.matrix,adult.Heart)
tissue.matrix <- cbind(tissue.matrix,adult.FrontalLobe)
tissue.matrix <- cbind(tissue.matrix,adult.Cortex)
tissue.matrix <- cbind(tissue.matrix,adult.Bladder)
tissue.matrix <- cbind(tissue.matrix,adult.Placenta)
tissue.matrix <- cbind(tissue.matrix,fetal.Liver)
tissue.matrix <- cbind(tissue.matrix,adult.Cerebellum)
tissue.matrix <- cbind(tissue.matrix,fetal.Limb)
tissue.matrix <- cbind(tissue.matrix,fetal.CNS.E14)
tissue.matrix <- cbind(tissue.matrix,fetal.CNS.E18)
tissue.matrix <- cbind(tissue.matrix,fetal.Liver.E14.5)
tissue.matrix <- cbind(tissue.matrix,fetal.WholeBrain_E14.5)
tissue.matrix <- cbind(tissue.matrix,fetal.CNS.E11.5)
tissue.matrix <- cbind(tissue.matrix,fetal.Liver.E14)
tissue.matrix <- cbind(tissue.matrix,adult.fibroblast)

dim(tissue.matrix)

#---
# index of heart tissue is 17
#---

tissue.max <- apply(tissue.matrix, 1, which.max)
tissue.max[tissue.max != 17] <- 0
tissue.max[tissue.max == 17] <- 1


#---
# index of fibroblast is 31
#---
tissue.max <- apply(tissue.matrix, 1, which.max)
tissue.max[tissue.max != 31] <- 0
tissue.max[tissue.max == 31] <- 1

#------------------------------------------------------------------------
# PMID: 15388519 
# Bioinformatics. 2005 Mar 1;21(5):650-9. Epub 2004 Sep 23.
# Genome-wide midrange transcription profiles reveal expression 
# level relationships in human tissue specification.
# see also:
# PMID: 26891983
# see also:
# https://www.biostars.org/p/209984/
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

tao.index    <- apply(tissue.matrix, 1, tao.func)
tissue.score <- tao.index * tissue.max

hisat2.data[,1][tissue.score != 0 & !is.na(tissue.score)]
heatmap.matrix <- tissue.matrix[tissue.score != 0 & !is.na(tissue.score),]
heatmap.result <- heatmap.2( heatmap.matrix, col = greenred(75),scale  = 'row', 
						     Rowv = TRUE,Colv = FALSE, density.info = 'none',key = TRUE, trace = 'none', 
						     cexCol = 0.6,distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						     hclustfun = function(d) hclust(d, method = 'complete'),
						     dendrogram = 'row',margins = c(12,9),labRow = NA, srtCol = 30,
						     lmat = rbind(c(4,0), c(2,1),c(0,3)), lhei = c(1,3, 0.5), lwid = c(1,4));
match('Nppa',rownames(hisat2.data))
setwd('/home/zhenyisong/data/cardiodata')
save.image('tissue-specific.Rdata')
q("no")