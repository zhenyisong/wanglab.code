#
# author Yisong Zhen
# since 2016-08-23
# update
# version 1.001
# aim
#    to extract genes of maturation and immaturation
#    of cardiomyocytes from public RNA-seq data
#
# raw data is from Boyer's lab: GSE64403
#  and GSE47948
# PubMed ID:
#     25477501
#     22981692
# cd /home/zhenyisong/data/cardiodata/SRP051406
# cp -r S*/*.sra ./
# https://www.biostars.org/p/156909/
# fastq-dump.2.4.5  --split-3 *.sra
# mkdir results
# ls *.fastq>targets.txt

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
library(pd.mogene.1.0.st.v1)
library(mogene10sttranscriptcluster.db)
library(oligo)
library(VennDiagram)

#setwd("C:\\Users\\Yisong\\Desktop")
#setwd("/home/zhenyisong/data/cardiodata")
#load("fibroblast.Rdata")
"
read GSE58453/SRP043191
"
setwd("/home/zhenyisong/data/cardiodata/SRP043191")
# ln -s /home/zhenyisong/data/cardiodata/SRP051406/mm10_index* ./
# ls *.fastq > targets.txt
targets.file      <- '/home/zhenyisong/data/cardiodata/SRP043191/targets.txt'
reads.files       <- read.table(targets.file,header = F)


"
output path, where the resuls are saved
"
reads.path        <- '/home/zhenyisong/data/cardiodata/SRP043191/'
output.path       <- '/home/zhenyisong/data/cardiodata/SRP043191/results/'

reads.files.names <- reads.files$V1
reads.paths       <- paste0(reads.path,reads.files$V1)
outputs.files     <- paste0(output.path,reads.files$V1,'.sam')



base.string       <- 'mm10_index'

align( index         = base.string, 
       readfile1     = reads.paths, 
       input_format  = "FASTQ", 
       type          = 'rna',
       output_file   = outputs.files, 
       output_format = "SAM", 
       nthreads      = 8, 
       indels        = 1,
       maxMismatches = 3,
       phredOffset   = 33,
       unique        = T )


polya.gene <- featureCounts( outputs.files, useMetaFeatures = TRUE, 
                               annot.inbuilt = "mm10", allowMultiOverlap = TRUE)



"
GSE49906/SRP029464
"
targets.file      <- '/home/zhenyisong/data/cardiodata/SRP029464/targets.txt'
reads.files       <- read.table(targets.file,header = F)


"
output path, where the resuls are saved
"
reads.path        <- '/home/zhenyisong/data/cardiodata/SRP029464/'
output.path       <- '/home/zhenyisong/data/cardiodata/SRP029464/results/'

reads.files.names <- reads.files$V1
read.path.1       <- reads.files.names[grep("_1",reads.files.names)]
read.path.2       <- reads.files.names[grep("_2",reads.files.names)]

"
generate the path vectors
"
reads.paths.1       <- paste0(reads.path,read.path.1)
reads.paths.2       <- paste0(reads.path,read.path.2)
outputs.files       <- paste0(output.path,read.path.1,'.bam')

"
the base index name
"
base.string       = 'mm10_index'

"
use the Rsubread command to generate index file
this index file will be generated and saved at getwd()
you do not need to generate the script
"
setwd("/home/zhenyisong/data/cardiodata/SRP029464")
#ln -s /home/zhenyisong/data/cardiodata/SRP051406/mm10_index* ./

duncan.gene <- featureCounts( outputs.files, useMetaFeatures = TRUE, 
                               annot.inbuilt = "mm10", allowMultiOverlap = TRUE)


setwd("D:\\wangli_data\\GSE14414")
#
# I extract those information from here
# https://support.bioconductor.org/p/59396/
#
celPath   <- list.files( getwd(), pattern = "\\.CEL|\\.CEL\\.gz", 
                      full.names = TRUE, ignore.case = TRUE)
raw.data  <- read.celfiles(celPath, pkgname = "pd.mogene.1.0.st.v1")
rma.data  <- rma(raw.data,target = 'core')
annotateGene <- function ( db , what , missing ) {
    tab <- toTable(db[intersect(featureNames(rma.data),mappedkeys(db)) ])
    mt  <- match ( featureNames (rma.data) , tab$probe_id )
    ifelse ( is.na(mt), missing , tab[[ what ]][ mt ])
}
fData(rma.data)$symbol <- annotateGene( mogene10sttranscriptclusterSYMBOL,"symbol" , missing = NA )
gene.symbol    <- annotateGene( mogene10sttranscriptclusterSYMBOL,"symbol" , missing = NA )
exprs.data <- exprs(rma.data)

adult1.heart <- apply(exprs.data[,1:3],1,median)
heart1.ord   <- order(adult1.heart, decreasing = TRUE)
adult1.fibro <- apply(exprs.data[,13:15],1,median)
fibro1.ord   <- order(adult1.fibro, decreasing = TRUE)

N <- as.integer(length(gene.symbol) * 0.05)
heart1.top.name <- gene.symbol[heart1.ord[1:N]]
fibro1.top.name <- gene.symbol[fibro1.ord[1:N]]



setwd("D:\\wangli_data\\GSE49192")
celPath   <- list.files( getwd(), pattern = "\\.CEL|\\.CEL\\.gz", 
                         full.names = TRUE, ignore.case = TRUE)
raw.data  <- read.celfiles(celPath[c(30:32,35:37)], pkgname = "pd.mogene.1.0.st.v1")
rma.data  <- rma(raw.data,target = 'core')
annotateGene <- function ( db , what , missing ) {
    tab <- toTable(db[intersect(featureNames(rma.data),mappedkeys(db)) ])
    mt  <- match ( featureNames (rma.data) , tab$probe_id )
    ifelse ( is.na(mt), missing , tab[[ what ]][ mt ])
}
fData(rma.data)$symbol <- annotateGene( mogene10sttranscriptclusterSYMBOL,"symbol" , missing = NA )
gene.symbol    <- annotateGene( mogene10sttranscriptclusterSYMBOL,"symbol" , missing = NA )
exprs.data     <- exprs(rma.data)

adult2.heart <- apply(exprs.data[,1:3],1,median)
heart2.ord   <- order(adult2.heart, decreasing = TRUE)
adult2.fibro <- apply(exprs.data[,4:6],1,median)
fibro2.ord   <- order(adult2.fibro, decreasing = TRUE)

N <- as.integer(length(gene.symbol) * 0.05)
heart2.top.name <- gene.symbol[heart2.ord[1:N]]
fibro2.top.name <- gene.symbol[fibro2.ord[1:N]]

#--------------------------------------------------------
#
#  extract ploya.gene
#
#--------------------------------------------------------

gene         <- polya.gene
gene.counts  <- gene$counts
gene.ids     <- gene$annotation$GeneID
colnames(gene.counts)
colnames(gene.counts) <- c('FVB-CM_polyA_1','FVB-CM_polyA_2','FVB-CM_polyA_3',
                           'FVB-fibro_polyA_1','FVB-fibro_polyA_2','FVB-fibro_polyA_3',
                           'FVB-CM_smallRNA_1','FVB-CM_smallRNA_2','FVB-CM_smallRNA_3',
                           'FVB-fibro_smallRNA_1','FVB-fibro_smallRNA_2','FVB-fibro_smallRNA_3')


keytypes(org.Mm.eg.db)

columns  <- c("ENTREZID","SYMBOL", "MGI", "GENENAME");
GeneInfo <- select( org.Mm.eg.db, keys= as.character(gene.ids), 
                   keytype="ENTREZID", columns = columns);
m        <- match(gene$annotation$GeneID, GeneInfo$ENTREZID);
Ann      <- cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                          GeneInfo[m, c("SYMBOL", "MGI", "GENENAME")]);

rownames(gene.counts) <- GeneInfo[m,'SYMBOL'];

Ann$Chr  <-  unlist( lapply(strsplit(Ann$Chr, ";"), 
                    function(x) paste(unique(x), collapse = "|")))
Ann$Chr  <- gsub("chr", "", Ann$Chr)
gene.exprs <- DGEList(counts = gene.counts, genes = Ann)
gene.exprs <- calcNormFactors(gene.exprs)
dge.tmm                  = t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
#dge.tmm.counts <- round(dge.tmm, digits = 0)
dge.tmm.counts           <- apply(dge.tmm,2, as.integer)

sample.info              <- data.frame( treat  = c('CM','CM','CM','fibro','fibro','fibro',
                                                   'CM_smallRNA','CM_smallRNA','CM_smallRNA',
                                                   'fibro_smallRNA','fibro_smallRNA','fibro_smallRNA') )
                                                
dds                      <- DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                                    colData   = sample.info,
                                                    design    = ~ treat)
vsd                      <- varianceStabilizingTransformation(dds, blind = FALSE);
vsd.expr                 <- assay(vsd)
vsd.expr                 <- vsd.expr[,1:6]
gene.symbol              <- gene.exprs$genes$SYMBOL
colnames(vsd.expr)       <- colnames(gene.counts)[1:6]
exprs.data               <- vsd.expr

adult3.heart <- apply(exprs.data[,1:3],1,median)
heart3.ord   <- order(adult3.heart, decreasing = TRUE)
adult3.fibro <- apply(exprs.data[,4:6],1,median)
fibro3.ord   <- order(adult3.fibro, decreasing = TRUE)

N <- as.integer(length(gene.symbol) * 0.05)
heart3.top.name <- gene.symbol[heart3.ord[1:N]]
fibro3.top.name <- gene.symbol[fibro3.ord[1:N]]


#--------------------------------------------------------
#
#  extract duncan.gene
#
#--------------------------------------------------------

gene         <- duncan.gene
gene.counts  <- gene$counts
gene.ids     <- gene$annotation$GeneID
colnames(gene.counts)
colnames(gene.counts) <- c('ventricle-PN90','ventricle-PN28','ventricle-PN10',
                           'ventricle-PN1','ventricle-E17',' fibroblasts-PN60',
                           'fibroblasts-PN28','fibroblasts-PN1-3','fibroblasts-PN1-2',
                           'cardiomyocytes-PN67','cardiomyocytes-PN30','cardiomyocytes-PN1-2',
                           'cardiomyocytes-PN1')

gene.counts <- gene.counts[,c(1,2,6,7)]
keytypes(org.Mm.eg.db)

columns  <- c("ENTREZID","SYMBOL", "MGI", "GENENAME");
GeneInfo <- select( org.Mm.eg.db, keys= as.character(gene.ids), 
                   keytype="ENTREZID", columns = columns);
m        <- match(gene$annotation$GeneID, GeneInfo$ENTREZID);
Ann      <- cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                          GeneInfo[m, c("SYMBOL", "MGI", "GENENAME")]);

rownames(gene.counts) <- GeneInfo[m,'SYMBOL'];

Ann$Chr  <-  unlist( lapply(strsplit(Ann$Chr, ";"), 
                    function(x) paste(unique(x), collapse = "|")))
Ann$Chr  <- gsub("chr", "", Ann$Chr)
gene.exprs <- DGEList(counts = gene.counts, genes = Ann)
gene.exprs <- calcNormFactors(gene.exprs)
dge.tmm                  = t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
#dge.tmm.counts <- round(dge.tmm, digits = 0)
dge.tmm.counts           <- apply(dge.tmm,2, as.integer)

sample.info              <- data.frame( treat  = c('heart','heart','fibro','fibro') )
                                                
dds                      <- DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                                    colData   = sample.info,
                                                    design    = ~ treat)
vsd                      <- varianceStabilizingTransformation(dds, blind = FALSE);
vsd.expr                 <- assay(vsd)
gene.symbol              <- gene.exprs$genes$SYMBOL
colnames(vsd.expr)       <- colnames(gene.counts)[1:4]
exprs.data               <- vsd.expr

adult4.heart <- apply(exprs.data[,1:2],1,median)
heart4.ord   <- order(adult4.heart, decreasing = TRUE)
adult4.fibro <- apply(exprs.data[,3:4],1,median)
fibro4.ord   <- order(adult4.fibro, decreasing = TRUE)

N <- as.integer(length(gene.symbol) * 0.05)
heart4.top.name <- gene.symbol[heart4.ord[1:N]]
fibro4.top.name <- gene.symbol[fibro4.ord[1:N]]


heart.name <- intersect(heart1.top.name,heart2.top.name)
heart.name <- intersect(heart.name, heart3.top.name)
heart.name <- intersect(heart.name, heart4.top.name)


fibro.name <- intersect(fibro1.top.name,fibro2.top.name)
fibro.name <- intersect(fibro.name, fibro3.top.name)
fibro.name <- intersect(fibro.name, fibro4.top.name)

area1 <- length(heart1.top.name)
area2 <- length(heart2.top.name)
area3 <- length(heart3.top.name)
area4 <- length(heart4.top.name)
n12   <- length(intersect(heart1.top.name,heart2.top.name))
n13   <- length(intersect(heart1.top.name,heart3.top.name))
n14   <- length(intersect(heart1.top.name,heart4.top.name))
n23   <- length(intersect(heart2.top.name,heart3.top.name))
n24   <- length(intersect(heart2.top.name,heart4.top.name))
n34   <- length(intersect(heart3.top.name,heart4.top.name))
n123  <- intersect(heart1.top.name,heart2.top.name)
n123  <- intersect(n123,heart3.top.name)
n123  <- length(n123)
n124  <- intersect(heart1.top.name,heart2.top.name)
n124  <- intersect(n124,heart4.top.name)
n124  <- length(n124)
n134  <- intersect(heart1.top.name,heart3.top.name)
n134  <- intersect(n134,heart4.top.name)
n134  <- length(n134)
n234  <- intersect(heart2.top.name,heart3.top.name)
n234  <- intersect(n234,heart4.top.name)
n234  <- length(n234)
n1234  <- intersect(heart1.top.name,heart2.top.name)
n1234  <- intersect(n1234,heart3.top.name)
n1234  <- intersect(n1234,heart4.top.name)
n1234  <- length(n1234)
draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24,
n34, n123, n124, n134, n234, n1234,alpha = rep(0.5, 4),fill = c('red','blue','green','purple'))cat.col = rep("black", 4) )  


category = rep("",
4), lwd = rep(2, 4), lty = rep("solid", 4), col =
rep("black", 4), fill = NULL, ,
label.col = rep("black", 15), cex = rep(1, 15),
fontface = rep("plain", 15), fontfamily = rep("serif",
15), cat.pos = c(-15, 15, 0, 0), cat.dist = c(0.22,
0.22, 0.11, 0.11), , cat.cex
= rep(1, 4), cat.fontface = rep("plain", 4),
cat.fontfamily = rep("serif", 4), cat.just =
rep(list(c(0.5, 0.5)), 4), rotation.degree = 0,
rotation.centre = c(0.5, 0.5), ind = TRUE, cex.prop =
NULL, print.mode = "raw", sigdigs = 3, direct.area =
FALSE, area.vector = 0, ...)
setwd("/home/zhenyisong/data/cardiodata")
save.image(file = 'fibroblast.Rdata')
quit("no")