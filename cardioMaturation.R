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


##---
## strandness check with 
##    http://rseqc.sourceforge.net/
##    the folowing data was downloaded from above site
##    cooper.data
##    GSE49906
##    SRP029464
##    PMID: 24752171
##    Title:
##    Alternative splicing regulates vesicular trafficking genes in 
##    cardiomyocytes during postnatal heart development. 
##    Nat Commun 2014 Apr 22;5:3603. 
##    mm10_bed='/home/zhenyisong/data/genome/mm10_RefSeq.bed'
##    cd /home/zhenyisong/data/cardiodata/SRP029464/results/
##    infer_experiment.py -i SRR964790_1.fastq.bam -r $mm10_bed
##    output:
##    This is PairEnd Data
##    Fraction of reads explained by "1++,1--,2+-,2-+": 0.9405
##    Fraction of reads explained by "1+-,1-+,2++,2--": 0.0243
##    Fraction of reads explained by other combinations: 0.0352
##    dataset SRP029464 use the parameter: FR
##    perrelo.data
##    cd /home/zhenyisong/data/cardiodata/SRP045149/results/
##    infer_experiment.py -i SRR1533761.fastq.sam  -r $mm10_bed
##    output
##    This is SingleEnd Data
##    Fraction of reads explained by "++,--": 0.4747
##    Fraction of reads explained by "+-,-+": 0.4986
##    Fraction of reads explained by other combinations: 0.0267
##    naked eye check
##        try to extract paired reads from raw data
##        sed 's/@//g;s/ /_/g' SRR4044044_1.fastq  | awk '{ if(NR%4 == 1) print ">"$0;if (NR%4 == 2) print $0; }' | head
##   and use UCSC Genome Browser to view the direction.

##
# http://stackoverflow.com/questions/9917049/inserting-an-image-to-ggplot2
# layout in the final output
##

library(gplots)
library(xlsx)
library(Rsubread)
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(DESeq2)
library(genefilter)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(affy)
library(annotate)
library(org.Hs.eg.db)
library(mgu74a.db)
library(RFLPtools)
library(VennDiagram)
library(DiagrammeR)
library(magrittr)
library(cowplot)
library(clusterProfiler)
library(png)
library(grid)
library(GGally)
library(reshape)


setwd("C:\\Users\\Yisong\\Desktop")
#setwd("/home/zhenyisong/data/cardiodata")
load("maturation.update.Rdata")
#load("GSE1479.Rdata")
#save.image("maturation.Rdata")

##"
##affy GSE75
##"
##
##setwd('/home/zhenyisong/data/cardiodata/GSE75')
##setwd("E:\\FuWai\\PaperPublished\\CardiacAgeing\\GSE75_RAW")
##raw.data       <- ReadAffy();
##rma.data       <- rma(raw.data);
##exprs.data     <- exprs(rma.data)


"
define the path where the reference genome sequence is deposited
"
genome_ref.path     <- "/home/zhenyisong/data/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
gtf_annotation.file <- "/home/zhenyisong/data/genome/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"

"
this is the targets file which indicate the fastq file path and other experiment
inforamtion regarding the sequence
ls *.fastq>targets.txt
"
targets.file      <- '/home/zhenyisong/data/cardiodata/SRP051406/targets.txt'
reads.files       <- read.table(targets.file,header = F)


"
output path, where the resuls are saved
"
reads.path        <- '/home/zhenyisong/data/cardiodata/SRP051406/'
output.path       <- '/home/zhenyisong/data/cardiodata/SRP051406/results/'

reads.files.names <- reads.files$V1
read.path.1       <- reads.files.names[grep("_1",reads.files.names)]
read.path.2       <- reads.files.names[grep("_2",reads.files.names)]

"
generate the path vectors
"
reads.paths.1     <- paste0(reads.path,read.path.1)
reads.paths.2     <- paste0(reads.path,read.path.2)
outputs.files     <- paste0(output.path,read.path.1,'.bam')

"
the base index name
"
base.string       = 'mm10_index'

"
use the Rsubread command to generate index file
this index file will be generated and saved at getwd()
you do not need to generate the script
"
setwd("/home/zhenyisong/data/cardiodata/SRP051406")
buildindex( basename = base.string, reference = genome_ref.path )

"
this is the function which is called to align the genome
sequence
"

align( index          = base.string, 
       readfile1      = reads.paths.1, 
       readfile2      = reads.paths.2, 
       input_format   = "FASTQ", 
       type           = 'rna',
       output_file    = outputs.files, 
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 8, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )


# get gene's counts
boyer.gene       <- featureCounts( outputs.files, useMetaFeatures = TRUE,
                                   countMultiMappingReads = FALSE,
                                   strandSpecific         = 0, 
                                   isPairedEnd            = TRUE,
                                   requireBothEndsMapped  = TRUE,
                                   autosort               = TRUE,
                                   nthreads               = 8,
                                   annot.inbuilt = "mm10", allowMultiOverlap = TRUE)

"
read GSE59970/SRP045149
PMID:?
"
setwd("/home/zhenyisong/data/cardiodata/SRP045149")
# ln -s /home/zhenyisong/data/cardiodata/SRP051406/mm10_index* ./
targets.file      <- '/home/zhenyisong/data/cardiodata/SRP045149/targets.txt'
reads.files       <- read.table(targets.file,header = F)


"
output path, where the resuls are saved
"
reads.path        <- '/home/zhenyisong/data/cardiodata/SRP045149/'
output.path       <- '/home/zhenyisong/data/cardiodata/SRP045149/results/'

reads.files.names <- reads.files$V1
reads.paths       <- paste0(reads.path,reads.files$V1)
outputs.files     <- paste0(output.path,reads.files$V1,'.sam')



base.string       <- 'mm10_index'
setwd("/home/zhenyisong/data/cardiodata/SRP045149/")
buildindex( basename = base.string, reference = genome_ref.path )

align( index          = base.string, 
       readfile1      = reads.paths, 
       input_format   = "FASTQ", 
       type           = 'rna',
       output_file    = outputs.files, 
       output_format  = "SAM", 
       nthreads       = 8, 
       indels         = 1, 
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )


perrelo.gene    <- featureCounts( outputs.files, useMetaFeatures = TRUE,
                                   countMultiMappingReads = FALSE,
                                   strandSpecific         = 0, 
                                   nthreads               = 8,
                                   annot.inbuilt = "mm10", allowMultiOverlap = TRUE)



"
GSE49906/SRP029464
PMID:24752171
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
base.string         <- 'mm10_index'

"
use the Rsubread command to generate index file
this index file will be generated and saved at getwd()
you do not need to generate the script
"
setwd("/home/zhenyisong/data/cardiodata/SRP029464")
buildindex( basename = base.string, reference = genome_ref.path )
"
this is the function which is called to align the genome
sequence
"

align( index          = base.string, 
       readfile1      = reads.paths.1, 
       readfile2      = reads.paths.2, 
       input_format   = "FASTQ", 
       type           = 'rna',
       output_file    = outputs.files, 
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 8, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

cooper.gene    <-   featureCounts( outputs.files, useMetaFeatures = TRUE,
                                   countMultiMappingReads = FALSE,
                                   strandSpecific         = 1, 
                                   isPairedEnd            = TRUE,
                                   requireBothEndsMapped  = TRUE,
                                   autosort               = TRUE,
                                   nthreads               = 8,
                                   annot.inbuilt = "mm10", allowMultiOverlap = TRUE)

setwd("/home/zhenyisong/data/cardiodata")
save.image(file = 'maturation.update.Rdata')
quit("no")


#---
#  Now, extract data from above output
#  GSE64403/SRP051406
#---


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
gene.exprs <- DGEList(counts = gene.counts, genes = Ann)
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
vsd.lineage.exprs        <- vsd.expr[,c('iP0_1','iP0_2','iP0_3','iP4_1','iP4_2','iP4_3',
                                        'iD7S_1','iD7S_2','iD7S_3','ex0hr_1','ex0hr_2')]


#---
# the dicarded filtering strategy is much more vigour than expect
# thus leading to no genes in intersection result
#---

# please see the interpretation of the prcomp
# in 'An Introdcution to Statistical Learning'
# Page: 407-409


#---
# seperate the line size and point size
# Can ggplot2 control point size and line size 
# (lineweight) separately in one legend?
# How to scale the size of line and point separately in ggplot2
#---
vsd.no.na.subset.exprs           <- vsd.subset.exprs[!is.na(rownames(vsd.subset.exprs)),]
rownames(vsd.no.na.subset.exprs) <- rownames(vsd.subset.exprs)[!is.na(rownames(vsd.subset.exprs))]
boyer.PCA  <- prcomp(t(vsd.no.na.subset.exprs))
names(boyer.PCA)
boyer.cord <- as.data.frame(boyer.PCA$x)
#biplot(boyer.PCA, scale = 0)
pc.no <- dim(boyer.cord)[2]

pr.var <- boyer.PCA$sdev^2
pve    <- pr.var/sum(pr.var)
pve.df <- data.frame(variance = pve, pca = c(1:8))
ggplot(pve.df) +
        xlab('Principle Component') +
        ylab('Proportion of Variance Explained') +
        scale_x_continuous( breaks = c(1:8), labels = as.character(c(1:8), 
                            limits = as.character(c(1:8)))) +
        geom_point(aes(x = pca, y = variance), size = 3) +
        geom_line(aes(x = pca, y = variance), size = 0.8) +
        scale_linetype_discrete() +
        theme(legend.position="none") +
        theme_gray()
pvecum.df <- data.frame(variance = cumsum(pve), pca = c(1:8))
ggplot(pvecum.df) +
        xlab('Principle Component') +
        ylab('Culmulative Proportion of Variance Explained') +
        scale_x_continuous( breaks = c(1:8), labels = as.character(c(1:8), 
                            limits = as.character(c(1:8)))) +
        geom_point(aes(x = pca, y = variance), size = 3) +
        geom_line(aes(x = pca, y = variance),size = 0.8) +
        scale_linetype_discrete() +
        theme(legend.position="none")



boyer.cord$color.type <- factor(c(1,1,2,2,3,3,4,4))

ggplot(boyer.cord) + 
           geom_point(aes(x = PC1, y = PC2, color = color.type), size = 3) + 
           scale_colour_manual( name   = 'developement stages',
                                values = c("black", "blue", "red","green"),
                                labels = c( '0 day postnatal ventricle, P0',
                                            '4 day postnatal ventricle, P4',
                                             '7 day postnatal ventricle, P7','adult')) +
           theme_gray() + 
           theme(legend.position = c(0.65, 0.2),legend.title.align = 0.5)
           


rotations <- order(boyer.PCA$rotation[,1], decreasing = FALSE)
gene.len  <- 200
total.len <- length(rownames(vsd.no.na.subset.exprs))
mature.boyer.markers <- unique(c( rownames(vsd.no.na.subset.exprs)[rotations[1:gene.len]],
                                  rownames(vsd.no.na.subset.exprs)[rotations[(total.len - gene.len):total.len]] ) )

mature.markers.pca.exprs <- vsd.subset.exprs[mature.markers,]

heatmap.result <- heatmap.2( mature.markers.pca.exprs, col = greenred(75),scale  = 'row', 
						     Rowv = TRUE,Colv = FALSE, density.info = 'none',key = TRUE, trace = 'none', 
						     cexCol = 1.5,distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						     hclustfun = function(d) hclust(d, method = 'complete'),
						     dendrogram = 'row',margins = c(12,9),labRow = NA, srtCol = 30,
						     lmat = rbind(c(4,0), c(2,1),c(0,3)), lhei = c(1,3, 0.5), lwid = c(1,4));

mature.markers.pca.exprs <- vsd.lineage.exprs[mature.markers,]

#---
# Boyer data 
# analyse the same data set from boyer lab and 
# find the lineage specifc genes using PCA
# we would like to test if the data result can be replicated
#---
vsd.no.na.lineage.exprs           <- vsd.lineage.exprs[!is.na(rownames(vsd.lineage.exprs)),]
rownames(vsd.no.na.lineage.exprs) <- rownames(vsd.lineage.exprs)[!is.na(rownames(vsd.lineage.exprs))]
boyer.lineage.PCA  <- prcomp(t(vsd.no.na.lineage.exprs))
names(boyer.lineage.PCA)
boyer.lineage.cord <- as.data.frame(boyer.lineage.PCA$x)
#biplot(boyer.PCA, scale = 0)
pc.no <- dim(boyer.lineage.cord)[2]

pr.var <- boyer.lineage.PCA$sdev^2
pve    <- pr.var/sum(pr.var)
pve.df <- data.frame(variance = pve, pca = c(1:pc.no))
#---
# Figure 1B
#---
PCA.curve <- ggplot(pve.df) +
             xlab('Principle Component') +
             ylab('Proportion of Variance Explained') +
             scale_x_continuous( breaks = c(1:pc.no), labels = as.character(c(1:pc.no), 
                                 limits = as.character(c(1:pc.no)))) +
             geom_point(aes(x = pca, y = variance), size = 3) +
             geom_line(aes(x = pca, y = variance), size = 0.8) +
             scale_linetype_discrete() +
             theme_gray() + 
             theme(legend.position="none")

pvecum.df <- data.frame(variance = cumsum(pve), pca = c(1:pc.no))
ggplot(pvecum.df) +
        xlab('Principle Component') +
        ylab('Culmulative Proportion of Variance Explained') +
        scale_x_continuous( breaks = c(1:pc.no), labels = as.character(c(1:pc.no), 
                            limits = as.character(c(1:pc.no)))) +
        geom_point(aes(x = pca, y = variance), size = 3) +
        geom_line(aes(x = pca, y = variance),size = 0.8) +
        scale_linetype_discrete() +
        theme(legend.position="none")

boyer.lineage.cord$color.type<- factor(c(1,1,1,2,2,2,3,3,3,4,4))

#---
# Figure 1C?
#---

PCA.cluster <- ggplot(boyer.lineage.cord) + 
               geom_point(aes(x = PC1, y = PC2, color = color.type), size = 3) + 
               scale_colour_manual( name   = 'isolated cardiomyocyte developement stages',
                                    values = c("black", "blue", "red","green"),
                                    labels = c( '0 day isolated cardiomyocytes, P0',
                                                '4 day isolated cardiomyocytes, P4',
                                                '7 day isolated cardiomyocytes, P7','adult')) +
               theme_gray() + 
               theme(legend.position = 'bottom',legend.direction = 'vertical') +              
               guides( color = guide_legend(title.position = 'top') )

ggplot(boyer.lineage.cord) + 
           geom_point(aes(x = PC2, y = PC3, color = color.type), size = 3) + 
           scale_colour_manual( name   = 'isolated cardiomyocyte developement stages',
                                values = c("black", "blue", "red","green"),
                                labels = c( '0 day isolated cardiomyocytes, P0',
                                            '4 day isolated cardiomyocytes, P4',
                                            '7 day isolated cardiomyocytes, P7','adult')) +
           theme_gray() + 
           theme(legend.position = 'bottom',legend.direction = 'vertical') + 
           guides( color = guide_legend(title.position = 'top') )

rotations <- order(boyer.lineage.PCA$rotation[,1], decreasing = FALSE)
gene.len  <- 400
total.len <- length(rownames(vsd.no.na.lineage.exprs))
mature.lineage.boyer.markers <- unique(c( rownames(vsd.no.na.lineage.exprs)[rotations[1:gene.len]],
                                           rownames(vsd.no.na.lineage.exprs)[rotations[(total.len - gene.len):total.len]] ) )

#---
# redefine fetal/mature genes
#---



# check if they are fetal genes by the definition

sds <- rowSds(vsd.no.na.lineage.exprs)
hist(sds)
sh  <- shorth(sds)
lineage.boyer.markers            <- vsd.no.na.lineage.exprs[sds > 0.3,]
lineage.boyer.tree               <- hclust(as.dist(1 - cor(t(lineage.boyer.markers), method='pearson')), method='complete');
mature.boyer.lineage.result      <- cutree(lineage.boyer.tree, k = 2)
mature.boyer.lineage.markers     <- names(mature.boyer.lineage.result)[mature.boyer.lineage.result == 2]
fetal.boyer.lineage.markers      <- names(mature.boyer.lineage.result)[mature.boyer.lineage.result == 1]

fetal.genes               <- intersect(mature.lineage.boyer.markers ,fetal.boyer.lineage.markers)
mature.genes              <- intersect(mature.lineage.boyer.markers ,mature.boyer.lineage.markers)
fetal.lineage.pca.exprs   <- vsd.no.na.lineage.exprs[fetal.genes,]
mature.lineage.pca.exprs  <- vsd.no.na.lineage.exprs[mature.genes,]
heatmap.result <- heatmap.2( fetal.lineage.pca.exprs , col = greenred(75), scale  = 'row', 
						         Rowv = TRUE,Colv = FALSE, density.info = 'none',
                                 key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
						         cexRow  = 1, cexCol = 1,srtCol = 30,
                                 distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						         hclustfun  = function(d) hclust(d, method = 'complete'),
						         dendrogram = 'no',margins = c(12,6),labRow = NA);
#---
# export the data to Desktop according to
# the rquirement by Wang Yin 2016-12-29
#---
setwd("C:\\Users\\Yisong\\Desktop")
write.xlsx(fetal.lineage.pca.exprs, file = 'fetal.genes.xlsx')
write.xlsx(mature.lineage.pca.exprs, file = 'mature.genes.xlsx')
##---
## this is the test Venn's graph
##---

area1 <- length(mature.boyer.markers)
area2 <- length(mature.lineage.boyer.markers) 
cross.area <- length(intersect(mature.boyer.markers,mature.lineage.boyer.markers))

setwd('C:/Users/Yisong/Desktop')
common.common <- intersect(mature.boyer.markers,mature.lineage.boyer.markers)
write.xlsx(file = 'common_set.xlsx', common.common)

draw.pairwise.venn( area1, area2, cross.area,
                    category = c('ventricle','cardiomyocyte'),
                    fill     = c('purple','blue'),
                    alpha    = 0.8,
                    ext.text = F,
                    cat.col  = c('green','green'),
                    cat.cex  = c(1.5,1.5),
                    cat.just = list(c(-0.5,2),c(1.2,1)),
                    label.col= c('yellow','yellow','yellow'),
                    lwd      = c(0.5,2),
                    cex      = c(2,3,2)
                   )  


#---
#
#---
"
read the curated mature/immature dataset
"
setwd("E:\\FuWai\\wangli.lab")
curated.filename  <- 'cardio_manual_maturation_markers.xlsx'
curated.file.df   <- read.xlsx(curated.filename, sheetIndex = 1, header = TRUE)
commom.commom <- intersect(curated.file.df$geneSymbol, mature.lineage.boyer.markers)
common.matrix <- vsd.no.na.lineage.exprs[commom.commom,]
P0            <- apply(common.matrix[,1:3],1,median)
P4            <- apply(common.matrix[,4:6],1,median)
P7            <- apply(common.matrix[,7:9],1,median)
P60           <- apply(common.matrix[,10:11],1,median)


common.df   <- matrix()
common.df   <- cbind(P0,P4,P7,P60)
common.melt <- melt(common.df)
ggplot(data = common.melt, aes( x = factor(X2, levels = c('P0','P4','P7','P60')), 
                                y = value, group = factor(X1), colour = factor(X1))) + 
    geom_line(size = 1) +
    geom_point(size = 2) + 
    theme_grey() +
    xlab('Cardiomyocyte lineage stages') +
    ylab('Gene Expression Level in log transformation') +
    scale_colour_discrete(  name  = "Gene names in Consistcy")

#---
# KEGG analysis
#--- START

common.names   <- intersect(mature.boyer.markers,mature.lineage.boyer.markers)
gene.entrez.id <- unlist( mget(x = common.names, envir = org.Mm.egALIAS2EG) )

kegg.table     <- enrichKEGG( gene.entrez.id, organism = "mouse", 
                              pvalueCutoff  = 0.05, 
                              pAdjustMethod = "BH", 
                              qvalueCutoff  = 0.1, readable = TRUE)
kegg.result       <- summary(kegg.table)
kegg.qvalue       <- -log(kegg.result$qvalue)
kegg.pathway.name <- kegg.result$Description

#---
# KEGG has negative report on pahtway enrichment analysis
#--- END
"
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

"
alternative filtering method which is 
from the book biocondcutor case
"
sds <- rowSds(vsd.subset.exprs)
hist(sds)
sh  <- shorth(sds)
mature.boyer.markers <- vsd.subset.exprs[sds > 0.3,]

"
test <- vsd.subset.exprs[fetal.genes,]

heatmap.result <- heatmap.2( test, col = greenred(75),scale  = 'row', 
						     Rowv = TRUE,Colv = FALSE, density.info = 'none',key = TRUE, trace = 'none', 
						     cexCol = 1.5,distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						     hclustfun = function(d) hclust(d, method = 'complete'),
						     dendrogram = 'row',margins = c(12,9),labRow = NA, srtCol = 30,
						     lmat = rbind(c(4,0), c(2,1),c(0,3)), lhei = c(1,3, 0.5), lwid = c(1,4));
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
fetal.boyer.markers  <- names(mature.boyer.result)[mature.boyer.result == 2]

#---
# PCA analysis
# please see PMID: 24739965
# see this article attachment: 
#    Ranalysis_scRNAseq_E14-16-18-AT2_paper_corr
#    line: 276
#---



#-- end of Boyer


"
duncan set
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
                          GeneInfo[m, c("SYMBOL", "MGI", "GENENAME")])

Ann$Chr    <-  unlist( lapply(strsplit(Ann$Chr, ";"), 
                            function(x) paste(unique(x), collapse = "|")))
Ann$Chr    <- gsub("chr", "", Ann$Chr)
gene.exprs <- DGEList(counts = gene.counts, genes = Ann)
gene.exprs <- calcNormFactors(gene.exprs)
dge.tmm                  = t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
#dge.tmm.counts <- round(dge.tmm, digits = 0)
dge.tmm.counts           <- apply(dge.tmm,2, as.integer)

sample.info              <- data.frame( treat  = c( 'PN90','PN28','PN10','PN1','E17','fPN60',
                                                    'fPN28','fPN1-3','fPN1-2','cPN67',
                                                    'cPN30','cPN1-2','cPN1') )
                                                
dds                      <- DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                                    colData   = sample.info,
                                                    design    = ~ treat)
vsd                      <- varianceStabilizingTransformation(dds);
vsd.expr                 <- assay(vsd)
rownames(vsd.expr)       <- gene.exprs$genes$SYMBOL
#colnames(vsd.expr)       <- colnames(gene.counts)[1:5]
colnames(vsd.expr)       <- colnames(gene.counts)
duncan.lineage.epxrs     <- vsd.expr[,

"
f1 <- function(x) (IQR(x) > 0.5)
f2 <- pOverA(0.25, log2(100))
f3 <- function(x) (median(2^x) > 300)
f4 <- function(x) (shapiro.test(x)$p.value > 0.05)
f5 <- function(x) (sd(x)/abs(mean(x)) < 0.1 )
f6 <- function(x) (sqrt(10) * abs(mean(x))/sd(x) > qt(0.975,9))
ff <- filterfun(f1,f2,f3,f4,f5,f6)
mature.duncan.markers <- genefilter(vsd.expr, ff)
mature.duncan.markers <- vsd.expr[mature.duncan.markers,]
"

sds <- rowSds(vsd.expr)
hist(sds)
sh  <- shorth(sds)
summary(sds)
mature.duncan.markers <- vsd.expr[sds > 0.5,]

heatmap.result <- heatmap.2( mature.duncan.markers, col = greenred(75),scale  = 'row', 
						     Rowv = TRUE,Colv = FALSE, density.info = 'none',key = TRUE, trace = 'none', 
						     cexCol = 1.5,distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						     hclustfun = function(d) hclust(d, method = 'complete'),
						     dendrogram = 'row',margins = c(12,9),labRow = NA, srtCol = 30,
						     lmat = rbind(c(4,0), c(2,1),c(0,3)), lhei = c(1,3, 0.5), lwid = c(1,4));


mature.duncan.tree    <- hclust(as.dist(1 - cor(t(mature.duncan.markers), method='pearson')), method='complete');
mature.duncan.result  <- cutree(mature.duncan.tree, k = 2)
mature.duncan.markers <- names(mature.duncan.result)[mature.duncan.result == 2]
fetal.duncan.markers  <- names(mature.duncan.result)[mature.duncan.result == 1]


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

"
f1 <- function(x) (IQR(x) > 0.5)
f2 <- pOverA(0.25, log2(100))
f3 <- function(x) (median(2^x) > 300)
f4 <- function(x) (shapiro.test(x)$p.value > 0.05)
f5 <- function(x) (sd(x)/abs(mean(x)) < 0.1 )
f6 <- function(x) (sqrt(10) * abs(mean(x))/sd(x) > qt(0.975,9))
ff <- filterfun(f1,f2,f3,f4,f5,f6)
mature.perrelo.markers <- genefilter(vsd.expr, ff)
mature.perrelo.markers <- vsd.expr[mature.perrelo.markers,]
"

sds <- rowSds(vsd.expr)
hist(sds)
sh  <- shorth(sds)
summary(sds)
mature.perrelo.markers <- vsd.expr[sds > 0.3,]

heatmap.result <- heatmap.2( mature.perrelo.markers, col = greenred(75),scale  = 'row', 
						     Rowv = TRUE,Colv = FALSE, density.info = 'none',key = TRUE, trace = 'none', 
						     cexCol = 1.5,distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						     hclustfun = function(d) hclust(d, method = 'complete'),
						     dendrogram = 'row',margins = c(12,9),labRow = NA, srtCol = 30,
						     lmat = rbind(c(4,0), c(2,1),c(0,3)), lhei = c(1,3, 0.5), lwid = c(1,4));



mature.perrelo.tree    <- hclust(as.dist(1 - cor(t(mature.perrelo.markers), method='pearson')), method='complete');
mature.perrelo.result  <- cutree(mature.perrelo.tree, k = 2)
mature.perrelo.markers <- names(mature.perrelo.result)[mature.perrelo.result == 1]
fetal.perrelo.markers  <- names(mature.perrelo.result)[mature.perrelo.result == 2]

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

"
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
"

sds <- rowSds(exprs.data)
hist(sds)
sh  <- shorth(sds)
summary(sds)
mature.harvard.markers <- exprs.data[sds > 0.3,]
mature.harvard.name    <- gene.symbols[sds > 0.3]


heatmap.result <- heatmap.2( mature.harvard.markers, col = greenred(75),scale  = 'row', 
						     Rowv = TRUE,Colv = FALSE, density.info = 'none',key = TRUE, trace = 'none', 
						     cexCol = 1.5,distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						     hclustfun = function(d) hclust(d, method = 'complete'),
						     dendrogram = 'row',margins = c(12,9),labRow = NA, srtCol = 30,
						     lmat = rbind(c(4,0), c(2,1),c(0,3)), lhei = c(1,3, 0.5), lwid = c(1,4));
mature.harvard.tree      <- hclust(as.dist(1 - cor(t(mature.harvard.markers), method='pearson')), method='complete');
mature.harvard.result    <- cutree(mature.harvard.tree, k = 2)
mature.harvard.markers   <- mature.harvard.name[mature.harvard.result == 1]
fetal.harvard.markers    <- mature.harvard.name[mature.harvard.result == 2]


mature.genes <- intersect(mature.perrelo.markers,mature.perrelo.markers)
mature.genes <- intersect(mature.genes,mature.duncan.markers)
mature.genes <- intersect(mature.genes,mature.boyer.markers)
mature.genes <- mature.genes[!is.na(mature.genes)]

fetal.genes  <-  intersect(fetal.perrelo.markers,fetal.perrelo.markers)
fetal.genes  <-  intersect(fetal.genes,fetal.duncan.markers)
fetal.genes  <-  intersect(fetal.genes,fetal.boyer.markers)
fetal.genes  <- fetal.genes[!is.na(fetal.genes)]

# mature genes ven's graph
area1 <- length(mature.perrelo.markers)
area2 <- length(mature.duncan.markers)
area3 <- length(mature.boyer.markers)
n12   <- length(intersect(mature.perrelo.markers,mature.duncan.markers))
n23   <- length(intersect(mature.boyer.markers,mature.duncan.markers))
n13   <- length(intersect(mature.boyer.markers,mature.perrelo.markers))
n123  <- intersect(mature.perrelo.markers,mature.duncan.markers)
n123  <- length( intersect(n123, mature.boyer.markers) )

draw.triple.venn(area1, area2, area3, n12, n23, n13, n123,
                 alpha = rep(0.5, 3),fill = c('red','blue','green'))

# fetal genes ven's graph
area1 <- length(fetal.perrelo.markers)
area2 <- length(fetal.duncan.markers)
area3 <- length(fetal.boyer.markers)
n12   <- length(intersect(fetal.perrelo.markers,fetal.duncan.markers))
n23   <- length(intersect(fetal.boyer.markers,fetal.duncan.markers))
n13   <- length(intersect(fetal.boyer.markers,fetal.perrelo.markers))
n123  <- intersect(fetal.perrelo.markers,fetal.duncan.markers)
n123  <- length( intersect(n123, fetal.boyer.markers) )

draw.triple.venn(area1, area2, area3, n12, n23, n13, n123,
                 alpha = rep(0.5, 3),fill = c('red','blue','green'))
setwd('/home/zhenyisong/data/cardiodata')
save.image(file = 'maturation.Rdata')
quit("no")


#--
# import hisat2 results for Cooper dataset
#--
setwd('/home/zhenyisong/data/cardiodata/SRP029464/hisat2')
file.name <- 'gene_count_matrix.csv'
cooper.hisat2.data <- read.csv( file.name, header = TRUE, sep = ",",
                                row.names = 1, stringsAsFactors = FALSE)
gene.counts     <- as.matrix(sapply(cooper.hisat2.data, as.integer))
gene.exprs      <- DGEList(counts = gene.counts)
gene.exprs      <- calcNormFactors(gene.exprs)
dge.tmm         <- t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
dge.tmm.counts  <- apply(dge.tmm,2, as.integer)

sample.info              <- data.frame( treat  = c( 'PN90','PN28','PN10','PN1','E17','fPN60',
                                                    'fPN28','fPN1-3','fPN1-2','cPN67',
                                                    'cPN30','cPN1-2','cPN1') )
                                                
dds                      <- DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                                    colData   = sample.info,
                                                    design    = ~ treat)
vsd                      <- varianceStabilizingTransformation(dds);
cooper.expr              <- assay(vsd)
colnames(cooper.expr)    <- c('ventricle-PN90','ventricle-PN28','ventricle-PN10',
                           'ventricle-PN1','ventricle-E17',' fibroblasts-PN60',
                           'fibroblasts-PN28','fibroblasts-PN1-3','fibroblasts-PN1-2',
                           'cardiomyocytes-PN67','cardiomyocytes-PN30','cardiomyocytes-PN1-2',
                           'cardiomyocytes-PN1')

cooper.no.na.lineage.exprs           <- cooper.expr[!is.na(rownames(cooper.hisat2.data)),c(10:13)]
rownames(cooper.no.na.lineage.exprs) <- rownames(cooper.hisat2.data)[!is.na(rownames(cooper.hisat2.data))]
cooper.lineage.PCA  <- prcomp(t(cooper.no.na.lineage.exprs))
names(cooper.lineage.PCA)
cooper.lineage.cord <- as.data.frame(cooper.lineage.PCA$x)
pc.no <- dim(cooper.lineage.cord)[2]

pr.var <- cooper.lineage.PCA$sdev^2
pve    <- pr.var/sum(pr.var)
pve.df <- data.frame(variance = pve, pca = c(1:pc.no))
cooper.lineage.cord$color.type<- factor(c(1,2,3,4))
ggplot(cooper.lineage.cord) + 
               geom_point(aes(x = PC1, y = PC2, color = color.type), size = 3) + 
               scale_colour_manual( name   = 'isolated cardiomyocyte developement stages',
                                    values = c("black", "blue", "red","green"),
                                    labels = c( 'cardiomyocytes-PN67',
                                                'cardiomyocytes-PN30',
                                                'cardiomyocytes-PN1-2',
                                                'cardiomyocytes-PN1')) +
               theme_gray() + 
               theme(legend.position = 'bottom',legend.direction = 'vertical') +              
               guides( color = guide_legend(title.position = 'top') )

#setwd('/home/zhenyisong/data/cardiodata')
#save.image(file = 'cooper.Rdata')
#quit("no")
setwd("C:\\Users\\Yisong\\Desktop")
load('cooper.Rdata')
rotations <- order(cooper.lineage.PCA$rotation[,1], decreasing = FALSE)
gene.len  <- 400
total.len <- length(rownames(cooper.no.na.lineage.exprs))
mature.lineage.cooper.markers <- unique(c( rownames(cooper.no.na.lineage.exprs)[rotations[1:gene.len]],
                                           rownames(cooper.no.na.lineage.exprs)[rotations[(total.len - gene.len):total.len]] ) )
# this trick is not to recyle
#  the common gene names in table manner
#-----
common.features <- sort(intersect(mature.lineage.boyer.markers,mature.lineage.cooper.markers))
length(common.features)       <- prod(dim(matrix(common.features, ncol = 8)))
common.matrix   <- matrix(common.features, ncol = 8, byrow = TRUE)
g <- tableGrob(common.matrix , rows = NULL)
grid.newpage()
grid.draw(g)
#---
# Venn's digram
#---
area1 <- length(mature.lineage.boyer.markers)
area2 <- length(mature.lineage.cooper.markers)
cross.area <- length(common.features)
draw.pairwise.venn( area1, area2, cross.area,
                    category = c('Boyer.data','Cooper.data'),
                    fill     = c('purple','blue'),
                    alpha    = 0.8,
                    ext.text = F,
                    cat.col  = c('green','green'),
                    cat.cex  = c(1.5,1.5),
                    cat.just = list(c(-0.5,2),c(1.2,1)),
                    label.col= c('yellow','yellow','yellow'),
                    lwd      = c(0.5,2),
                    cex      = c(2,3,2))
#
# flowchart creation
#


graph <-
    create_graph() %>%
    set_graph_name("Flow Chart of Mature gene selection") %>%
    set_global_graph_attrs("graph", "overlap", "true") %>%
    set_global_graph_attrs("node", "color", "grey75") %>%
    set_global_graph_attrs("node", "fontname", "Helvetica") %>%
    add_n_nodes(5) %>%
    select_nodes_by_id(c(1:5)) %>% 
    set_node_attrs_ws("shape", "retangle") %>%
    set_node_attrs_ws("style", "filled") %>%
    clear_selection %>%
    add_edges_w_string(
      "1->2 2->3 3->4 5->3", "black") %>%
    set_node_attrs("label",c( 'PCA analysis','loading selection',
                              'gene intersection','mature/immature\ngeneset',
                              'Hierarchical clustering')) %>%
    set_node_attrs("fontsize",12) %>%
    set_edge_attrs("arrowsize", 1)
render_graph(graph)
setwd("C:\\Users\\Yisong\\Desktop")
file.temp.name <- 'temp_zhen3.png'
export_graph(graph, file_name = file.temp.name, file_type = 'png')
img             <- readPNG('temp_zhen3.png')
PCA.img         <- rasterGrob(img, interpolate = TRUE)
#PCA.flowchart   <- qplot() + annotation_raster(img, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, interpolate = TRUE)
#PCA.flowchart   <- ggplot() + annotation_custom(PCA.flowchart , xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
PCA.flowchart <- ggplot(mtcars, aes(wt, mpg)) +
                 xlab(NULL) +
                 ylab(NULL) +
                 theme_void() +
                 geom_blank() + 
                 annotation_raster(img, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, interpolate = TRUE )
unlink(file.temp.name)
plot_grid(PCA.flowchart, PCA.curve, PCA.cluster, PCA.flowchart, labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2)
#---
# please see ino80.zn.R
#    and load the data in corresponding
#---

load("ino80.zn.Rdata")
fetal <- intersect(rownames(vsd.expr),toupper(fetal.genes))
fetal.heatmap       <- vsd.expr[fetal,]
ind <- apply(fetal.heatmap, 1, sd) == 0
subset <- fetal.heatmap[!ind,]
setwd("C:\\Users\\Yisong\\Desktop")
png(file.temp.name)
heatmap.result      <- heatmap.2(subset, col = my_palette, scale  = 'row', 
						         Rowv = TRUE,Colv = FALSE, density.info = 'none',
                                 key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
						         cexRow  = 1, cexCol = 1,srtCol = 30,
                                 distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						         hclustfun  = function(d) hclust(d, method = 'complete'),
						         dendrogram = 'no',margins = c(12,6),labRow = NA,
						        );
dev.off()
img             <- readPNG(file.temp.name)
PCA.heatmap     <- ggplot(mtcars, aes(wt, mpg)) +
                   theme_void() +
                   geom_blank() + 
                   annotation_raster(img , xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, interpolate = TRUE)
unlink(file = file.temp.name)

plot_grid(PCA.flowchart, PCA.curve, PCA.cluster, PCA.heatmap, labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2)

## end----


setwd("C:\\Users\\Yisong\\Desktop")

grid.arrange(PCA.flowchart, PCA.curve, PCA.cluster, PCA.heatmap, ncol=2, nrow =2)




library(ggplot2)
# Create a box plot
# Scatter plot
sp <- ggplot(mpg, aes(x = cty, y = hwy, colour = factor(cyl)))+ 
  geom_point(size=2.5)
sp
# Bar plot
bp <- ggplot(diamonds, aes(clarity, fill = cut)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle=70, vjust=0.5))

library(gridExtra)
grid.arrange(sp, bp, sp, g, ncol=2, nrow =2)

plot.A <- render_graph(graph, output = "SVG")

plot_grid(plot.B,plot.C, labels = c("A", "B"), ncol = 2, align = 'h')