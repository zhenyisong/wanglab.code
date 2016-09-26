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
# @reference
# https://www.biostars.org/p/86563/

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
library(affy)
library(annotate)
library(limma)
library(org.Hs.eg.db)
library(mgu74a.db)
library(RFLPtools)


setwd("C:\\Users\\Yisong\\Desktop")
#setwd("/home/zhenyisong/data/cardiodata")
load("maturation.Rdata")
#load("GSE1479.Rdata")
#save.image("maturation.Rdata")

"
affy GSE75
"

#setwd('/home/zhenyisong/data/cardiodata/GSE75')
setwd("E:\\FuWai\\PaperPublished\\CardiacAgeing\\GSE75_RAW")
raw.data       <- ReadAffy();
rma.data       <- rma(raw.data);
exprs.data     <- exprs(rma.data)
"
define the path where the reference genome sequence is deposited
"
genome_ref.path     <- "/home/zhenyisong/data/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
gtf_annotation.file <- "/home/zhenyisong/data/genome/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"

#"
#this is the targets file which indicate the fastq file path and other experiment
#inforamtion regarding the sequence
#ls *.fastq>targets.txt
#"
#targets.file      <- '/home/zhenyisong/data/cardiodata/SRP051406/targets.txt'
#reads.files       <- read.table(targets.file,header = F)
#
#
#"
#output path, where the resuls are saved
#"
#reads.path        <- '/home/zhenyisong/data/cardiodata/SRP051406/'
#output.path       <- '/home/zhenyisong/data/cardiodata/SRP051406/results/'
#
#reads.files.names <- reads.files$V1
#read.path.1       <- reads.files.names[grep("_1",reads.files.names)]
#read.path.2       <- reads.files.names[grep("_2",reads.files.names)]
#
#"
#generate the path vectors
#"
#reads.paths.1       <- paste0(reads.path,read.path.1)
#reads.paths.2       <- paste0(reads.path,read.path.2)
#outputs.files       <- paste0(output.path,read.path.1,'.bam')
#
#"
#the base index name
#"
#base.string       = 'mm10_index'
#
#"
#use the Rsubread command to generate index file
#this index file will be generated and saved at getwd()
#you do not need to generate the script
#"
#setwd("/home/zhenyisong/data/cardiodata/SRP051406")
##buildindex( basename = base.string, reference = genome_ref.path )
#
#"
#this is the function which is called to align the genome
#sequence
#"
#
#align( index         = base.string, 
#       readfile1     = reads.paths.1, 
#       readfile2     = reads.paths.2, 
#       input_format  = "FASTQ", 
#       type          = 'rna',
#       output_file   = outputs.files, 
#       output_format = "BAM", 
#       nthreads      = 8, 
#       indels        = 1,
#       maxMismatches = 3,
#       phredOffset   = 33,
#       unique        = T )
#
#
## get gene's counts
#boyer.gene       <- featureCounts( outputs.files, useMetaFeatures = TRUE, 
#                                   annot.inbuilt = "mm10", allowMultiOverlap = TRUE)
#
#"
#read GSE59970/SRP045149
#"
#setwd("/home/zhenyisong/data/cardiodata/SRP045149")
## ln -s /home/zhenyisong/data/cardiodata/SRP051406/mm10_index* ./
#targets.file      <- '/home/zhenyisong/data/cardiodata/SRP045149/targets.txt'
#reads.files       <- read.table(targets.file,header = F)
#
#
#"
#output path, where the resuls are saved
#"
#reads.path        <- '/home/zhenyisong/data/cardiodata/SRP045149/'
#output.path       <- '/home/zhenyisong/data/cardiodata/SRP045149/results/'
#
#reads.files.names <- reads.files$V1
#reads.paths       <- paste0(reads.path,reads.files$V1)
#outputs.files     <- paste0(output.path,reads.files$V1,'.sam')
#
#
#
#base.string       <- 'mm10_index'
#
#align( index         = base.string, 
#       readfile1     = reads.paths, 
#       input_format  = "FASTQ", 
#       type          = 'rna',
#       output_file   = outputs.files, 
#       output_format = "SAM", 
#       nthreads      = 8, 
#       indels        = 1,
#       maxMismatches = 3,
#       phredOffset   = 33,
#       unique        = T )
#
#
#perrelo.gene <- featureCounts( outputs.files, useMetaFeatures = TRUE, 
#                               annot.inbuilt = "mm10", allowMultiOverlap = TRUE)
#
#
#
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

"
this is the function which is called to align the genome
sequence
"

align( index         = base.string, 
       readfile1     = reads.paths.1, 
       readfile2     = reads.paths.2, 
       input_format  = "FASTQ", 
       type          = 'rna',
       output_file   = outputs.files, 
       output_format = "BAM", 
       nthreads      = 8, 
       indels        = 1,
       maxMismatches = 3,
       phredOffset   = 33,
       unique        = T )

duncan.gene <- featureCounts( outputs.files, useMetaFeatures = TRUE, 
                               annot.inbuilt = "mm10", allowMultiOverlap = TRUE)

#setwd("/home/zhenyisong/data/cardiodata")
#save.image(file = 'maturation.Rdata')
#quit("no")
#

#
#  Now, extract data from above output
#  GSE64403/SRP051406


gene         <- boyer.gene
gene.counts  <- gene$counts
gene.ids     <- gene$annotation$GeneID
colnames(gene.counts)
colnames(gene.counts) <- c( 'iP0_1','iP0_2','iP0_3','iP4_1','iP4_2','iP4_3',
                            'iD7S_1','iD7S_2','iD7S_3','iD7R_1','iD7R_2','iD7R_3',
                            'ex0hr_1','ex0hr_2','ex24hr_1','ex24hr_2','ex48hr_1','ex48hr_2',
                            'ex72hr_1','ex72hr_2','vP0_1','vP0_2','vP4_1','vP4_2',
                            'vP7_1','vP7_2','vAd_1','vAd_2','vD1S_1','vD1S_2','vD1R_1','vD1R_2',
                            'vD7R_1','vD7R_2','vD7S_1','vD7S_2')



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
gene.exprs <- DGEList(counts = boyer.gene, genes = Ann)
gene.exprs <- calcNormFactors(gene.exprs)
dge.tmm                  = t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
#dge.tmm.counts <- round(dge.tmm, digits = 0)
dge.tmm.counts           <- apply(dge.tmm,2, as.integer)

sample.info              <- data.frame( treat  = c('iP0','iP0','iP0','iP4','iP4','iP4',
                                                   'iD7S','iD7S','iD7S','iD7R','iD7R','iD7R',
                                                   'ex0hr','ex0hr','ex24hr','ex24hr','ex48hr','ex48hr',
                                                   'ex72hr','ex72hr','vP0','vP0','vP4','vP4','vP7','vP7',
                                                   'vAd','vAd','vD1S','vD1S','vD1R','vD1R','vD7R','vD7R',
                                                   'vD7S','vD7S') )
dds                      <- DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                                    colData   = sample.info,
                                                    design    = ~ treat)
vsd                      <- varianceStabilizingTransformation(dds, blind = FALSE);
vsd.expr                 <- assay(vsd)
rownames(vsd.expr)       <- gene.exprs$genes$SYMBOL
colnames(vsd.expr)       <- colnames(gene.counts)
vsd.subset.exprs         <- vsd.expr[,c('vP0_1','vP0_2','vP4_1','vP4_2','vP7_1','vP7_2','vAd_1','vAd_2')]

f1 <- function(x) (IQR(x) > 0.5)
f2 <- pOverA(0.25, log2(100))
f3 <- function(x) (median(2^x) > 300)
f4 <- function(x) (shapiro.test(x)$p.value > 0.05)
f5 <- function(x) (sd(x)/abs(mean(x)) < 0.1 )
f6 <- function(x) (sqrt(10) * abs(mean(x))/sd(x) > qt(0.975,9))
ff <- filterfun(f1,f2,f3,f4,f5,f6)
mature.boyer.markers <- genefilter(vsd.subset.exprs, ff)
mature.boyer.markers <- vsd.subset.exprs[mature.boyer.markers,]
"
sds <- rowSds(vsd.subset.exprs)
hist(sds)
sh  <- shorth(sds)
"


heatmap.result <- heatmap.2( mature.boyer.markers, col = greenred(75),scale  = 'row', 
						     Rowv = TRUE,Colv = FALSE, density.info = 'none',key = TRUE, trace = 'none', 
						     cexCol = 1.5,distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						     hclustfun = function(d) hclust(d, method = 'complete'),
						     dendrogram = 'row',margins = c(12,9),labRow = NA, srtCol = 30,
						     lmat = rbind(c(4,0), c(2,1),c(0,3)), lhei = c(1,3, 0.5), lwid = c(1,4));
mature.boyer.tree    <- hclust(as.dist(1 - cor(t(mature.boyer.markers), method='pearson')), method='complete');
mature.boyer.result  <- cutree(mature.boyer.tree, k = 2)
mature.boyer.markers <- names(mature.boyer.result)[mature.boyer.result == 1]

#-- end of Boyer


"
GSE49906,SRP029464
"
gene         <- duncan.gene
gene.counts  <- gene$counts
gene.ids     <- gene$annotation$GeneID
colnames(gene.counts)
colnames(gene.counts) <- c('ventricle-PN90','ventricle-PN28','ventricle-PN10',
                           'ventricle-PN1','ventricle-E17',' fibroblasts-PN60',
                           'fibroblasts-PN28','fibroblasts-PN1-3','fibroblasts-PN1-2',
                           'cardiomyocytes-PN67','cardiomyocytes-PN30','cardiomyocytes-PN1-2',
                           'cardiomyocytes-PN1')

gene.counts <- gene.counts[,1:5]


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

sample.info              <- data.frame( treat  = c('PN90','PN28','PN10','PN1','E17') )
                                                
dds                      <- DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                                    colData   = sample.info,
                                                    design    = ~ treat)
vsd                      <- varianceStabilizingTransformation(dds);
vsd.expr                 <- assay(vsd)
rownames(vsd.expr)       <- gene.exprs$genes$SYMBOL
colnames(vsd.expr)       <- colnames(gene.counts)[1:5]


f1 <- function(x) (IQR(x) > 0.5)
f2 <- pOverA(0.25, log2(100))
f3 <- function(x) (median(2^x) > 300)
f4 <- function(x) (shapiro.test(x)$p.value > 0.05)
f5 <- function(x) (sd(x)/abs(mean(x)) < 0.1 )
f6 <- function(x) (sqrt(10) * abs(mean(x))/sd(x) > qt(0.975,9))
ff <- filterfun(f1,f2,f3,f4,f5,f6)
mature.duncan.markers <- genefilter(vsd.expr, ff)
mature.duncan.markers <- vsd.expr[mature.duncan.markers,]

heatmap.result <- heatmap.2( mature.duncan.markers, col = greenred(75),scale  = 'row', 
						     Rowv = TRUE,Colv = FALSE, density.info = 'none',key = TRUE, trace = 'none', 
						     cexCol = 1.5,distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						     hclustfun = function(d) hclust(d, method = 'complete'),
						     dendrogram = 'row',margins = c(12,9),labRow = NA, srtCol = 30,
						     lmat = rbind(c(4,0), c(2,1),c(0,3)), lhei = c(1,3, 0.5), lwid = c(1,4));
mature.duncan.tree    <- hclust(as.dist(1 - cor(t(mature.duncan.markers), method='pearson')), method='complete');
mature.duncan.result  <- cutree(mature.duncan.tree, k = 2)
mature.duncan.markers <- names(mature.duncan.result)[mature.duncan.result == 1]

#-- end of duncan

#intersect(perrelo.markers,boyer.markers)

"
GSE59970/SRP045149
"
gene         <- perrelo.gene
gene.counts  <- gene$counts
gene.ids     <- gene$annotation$GeneID
colnames(gene.counts)
colnames(gene.counts) <- c('P1-1','P1-2','P1-3',
                           'P1-14','P1-14',' P1-14')


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

sample.info              <- data.frame( treat  = c('P1','P1','P1','P14','P14','P14'))
                                                
dds                      <- DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                                    colData   = sample.info,
                                                    design    = ~ treat)
vsd                      <- varianceStabilizingTransformation(dds, blind = FALSE);
vsd.expr                 <- assay(vsd)
rownames(vsd.expr)       <- gene.exprs$genes$SYMBOL
colnames(vsd.expr)       <- colnames(gene.counts)


f1 <- function(x) (IQR(x) > 0.5)
f2 <- pOverA(0.25, log2(100))
f3 <- function(x) (median(2^x) > 300)
f4 <- function(x) (shapiro.test(x)$p.value > 0.05)
f5 <- function(x) (sd(x)/abs(mean(x)) < 0.1 )
f6 <- function(x) (sqrt(10) * abs(mean(x))/sd(x) > qt(0.975,9))
ff <- filterfun(f1,f2,f3,f4,f5,f6)
mature.perrelo.markers <- genefilter(vsd.expr, ff)
mature.perrelo.markers <- vsd.expr[mature.perrelo.markers,]

heatmap.result <- heatmap.2( mature.perrelo.markers, col = greenred(75),scale  = 'row', 
						     Rowv = TRUE,Colv = FALSE, density.info = 'none',key = TRUE, trace = 'none', 
						     cexCol = 1.5,distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						     hclustfun = function(d) hclust(d, method = 'complete'),
						     dendrogram = 'row',margins = c(12,9),labRow = NA, srtCol = 30,
						     lmat = rbind(c(4,0), c(2,1),c(0,3)), lhei = c(1,3, 0.5), lwid = c(1,4));
mature.perrelo.tree    <- hclust(as.dist(1 - cor(t(mature.perrelo.markers), method='pearson')), method='complete');
mature.perrelo.result  <- cutree(mature.perrelo.tree, k = 2)
mature.perrelo.markers <- names(mature.perrelo.result)[mature.perrelo.result == 1]

"
GSE75
"
colnames(exprs.data) <- c('1w-1','1w-2','1w-3','1y-1','1y-2','1y-3','4w-1','4w-2','4w-3',
                          '5m-1','5m-2','5m-3','E12.5-1','E12.5-2','E12.5-3',
                          'NN-1','NN-2','NN-3','3m-1f','3m-2f','3m-3f','1y-1f','1y-2f','1y-3f')
exprs.data           <- exprs.data[,c('E12.5-1','E12.5-2','E12.5-3',
                                      'NN-1','NN-2','NN-3',
                                      '1w-1','1w-2','1w-3',
                                      '4w-1','4w-2','4w-3',
                                      '3m-1f','3m-2f','3m-3f',
                                      '5m-1','5m-2','5m-3',
                                      '1y-1f','1y-2f','1y-3f',
                                      '1y-1','1y-2','1y-3' )]
probes               <- rownames(exprs.data)
gene.symbols         <- unlist(mget(probes, mgu74aSYMBOL, ifnotfound = NA))

f1 <- function(x) (IQR(x) > 0.5)
f2 <- pOverA(0.25, log2(100))
f3 <- function(x) (median(2^x) > 300)
f4 <- function(x) (shapiro.test(x)$p.value > 0.05)
f5 <- function(x) (sd(x)/abs(mean(x)) < 0.1 )
f6 <- function(x) (sqrt(10) * abs(mean(x))/sd(x) > qt(0.975,9))
ff <- filterfun(f1,f2,f3,f4,f5,f6)

mature.harvard.markers <- genefilter(exprs.data, ff)
mature.harvard.name    <- gene.symbols[mature.harvard.markers]
mature.harvard.markers <- exprs.data[mature.harvard.markers,]


heatmap.result <- heatmap.2( mature.harvard.markers, col = greenred(75),scale  = 'row', 
						     Rowv = TRUE,Colv = FALSE, density.info = 'none',key = TRUE, trace = 'none', 
						     cexCol = 1.5,distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						     hclustfun = function(d) hclust(d, method = 'complete'),
						     dendrogram = 'row',margins = c(12,9),labRow = NA, srtCol = 30,
						     lmat = rbind(c(4,0), c(2,1),c(0,3)), lhei = c(1,3, 0.5), lwid = c(1,4));
mature.harvard.tree      <- hclust(as.dist(1 - cor(t(mature.harvard.markers), method='pearson')), method='complete');
mature.harvard.result    <- cutree(mature.harvard.tree, k = 2)
mature.harvard.markers   <- mature.harvard.name[mature.harvard.result == 1]


mature.genes <- intersect(mature.harvard.markers,mature.perrelo.markers)
mature.genes <- intersect(mature.perrelo.markers,mature.duncan.markers)
mature.genes <- intersect(mature.genes,mature.boyer.markers)

setwd('/home/zhenyisong/data/cardiodata')
save.image(file = 'maturation.Rdata')
quit("no")