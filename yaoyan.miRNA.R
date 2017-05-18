# @author Yisong zhen
# @since 2017-03-30
# @updated 2017-05-18
# @parent
#    yaoyan.miRNA.sh
# 
# ageing
# database
# http://senescence.info/
# I downloaded the human ageing related gene
# http://genomics.senescence.info/genes/microarray.php
# unzip and get the 

library(gplots)
library(xlsx)
library(Rsubread)
library(edgeR)
library(limma)
library(DESeq2)
library(genefilter)
library(RColorBrewer)
library(org.Mm.eg.db)
library(cluster)
library(factoextra)
library(clusterProfiler)
library(pathview)
library(sva)
library(systemPipeR)
library(rtracklayer)
library(stringr)
library(GenomicFeatures)
library(tidyverse)
library(magrittr)
library(GEOquery)
library(grid)
library(pheatmap)
library(VennDiagram)


# Agilen image raw data normlization method
#
# library(LVSmiRNA)

# read
setwd('/home/zhenyisong/biodata/wanglab/wangdata/yaoyan/rawdata/miRNA_bwa')
setwd('D:\\wangli_data\\Rdata')
load("yaoyan.miRNA.Rdata")
load("yaoyan.encode.Rdata")


miRNA.count.files  <- list.files(pattern = "Sample.*\\.txt$")
miRNA.count.matrix <- NULL
for( filename in miRNA.count.files) {
    column <- read.delim(filename, header = F, row.names = 1)
    miRNA.count.matrix <- cbind(miRNA.count.matrix , column$V2 )
}

rownames(miRNA.count.matrix) <- rownames(column)
colnames(miRNA.count.matrix) <- sub("\\.txt","",miRNA.count.files)

groups <- factor(c(rep(1,3),rep(2,3), rep(1,3),rep(2,3)), levels = 1:2, labels = c('Control','Patient'))
ages   <- factor(c(rep(1,6),rep(2,6)), levels = 1:2, labels = c('AF60','AF70'))
AF.miRNA.trans.data <- DGEList(counts = miRNA.count.matrix[,7:18], group = groups ) %>%
                       calcNormFactors(method = 'TMM') %>%
                       estimateCommonDisp() %>%
                       estimateTagwiseDisp(trend = 'movingave')

#
# depracated exactTest for simple analysis
#---
#AF.miRNA.edgeR               <- exactTest(AF.miRNA.trans.data)
#AF.miRNA.edgeR.pvalue        <- AF.miRNA.edgeR$table$PValue
#AF.miRNA.edgeR.p.adj         <- p.adjust(AF.miRNA.edgeR.pvalue, method = 'BH')

# introduce the ages parameter to cancel the batch effect?
#
design               <- model.matrix(~ 0 + groups + ages);
colnames(design)     <- c('Control','Patient','Age')
contrast.matrix      <- makeContrasts(Patient - Control, levels = design)
fit                  <- glmFit(AF.miRNA.trans.data, design)
model.result         <- glmLRT(fit, contrast = contrast.matrix[,1])
AF.miRNA.edgeR.pvalue        <- model.result$table$PValue
AF.miRNA.edgeR.p.adj         <- p.adjust(AF.miRNA.edgeR.pvalue, method = 'BH')


#---
#AF.miRNA.edgeR.sig.table     <- AF.miRNA.edgeR$table[which(AF.miRNA.edgeR.p.adj < 0.05),]
# this step output zero gene below the threshold
# I instead not correct the p.value using BH method
#---
AF.miRNA.edgeR.sig.table           <- model.result$table[which(AF.miRNA.edgeR.pvalue < 0.05),] %>%
                                      arrange(PValue)
rownames(AF.miRNA.edgeR.sig.table) <- rownames(model.result$table)[which(AF.miRNA.edgeR.pvalue < 0.05)]


# branched project
#---

gene.exprs.miRNA    <- DGEList(counts = miRNA.count.matrix[,c(1:9,13:15)]) %>% calcNormFactors()
dge.tmm.miRNA       <- t(t(gene.exprs.miRNA$counts) * gene.exprs.miRNA$samples$norm.factors) %>%
                       apply(2, as.integer)

rownames(dge.tmm.miRNA) <- rownames(column)
dge.tmm.miRNA.log       <- log(dge.tmm.miRNA + 1) 


sample.info              <- data.frame( treat  = c('AF40','AF40','AF40','AF50','AF50', 'AF50',
                                                   'AF60','AF60','AF60','AF70','AF70', 'AF70') )
dds.miRNA                <- DESeqDataSetFromMatrix( countData = dge.tmm.miRNA,
                                                    colData   = sample.info,
                                                    design    = ~ treat)
vsd.exprs.miRNA          <- varianceStabilizingTransformation(dds.miRNA, blind = FALSE) %>%
                            assay()

colnames(vsd.exprs.miRNA)   <- c('AF40-1','AF40-2','AF40-3','AF50-1','AF50-2', 'AF50-3',
                                 'AF60-1','AF60-2','AF60-3','AF70-1','AF70-2', 'AF70-3') 

aging.miRNA.sd              <- apply(vsd.exprs.miRNA, 1, sd)
hist(aging.miRNA.sd)
estimator.sd                <- shorth(aging.miRNA.sd)
aging.miRNA.filtered        <- vsd.exprs.miRNA[aging.miRNA.sd > 0.3,]

heatmap.result <- heatmap.2(  aging.miRNA.filtered  , col = greenred(75), scale  = 'row', 
						      Rowv = TRUE,Colv = FALSE, density.info = 'none',
                              key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
                              cexCol = 0.5,srtCol = 30,
                              distfun = function(d) as.dist(1 - cor(t(d), method = 'spearman')),
						      hclustfun  = function(d) hclust(d, method = 'complete'),
						      dendrogram = 'no', labRow = NA);

#aging.dist.miRNA            <- cor(t(aging.miRNA.filtered), method = "spearman") %>% as.dist()
aging.hclust.miRNA          <- hclust(d = aging.dist.miRNA, method = "complete")


aging.group          <- factor( c(rep(1:4, each = 3)), 
                                levels = 1:4, 
                                labels = c('AF40','AF50','AF60','AF70'))

miRNA.aging.data     <- DGEList( counts = miRNA.count.matrix[,c(1:9,13:15)], group = aging.group ) %>%
                                 calcNormFactors(method = 'TMM') %>%
                                 estimateCommonDisp() %>%
                                 estimateTagwiseDisp(trend = 'movingave')

design               <- model.matrix(~ 0 + aging.group);
colnames(design)     <- levels(aging.group)
contrast.matrix      <- makeContrasts(AF60 - AF40, AF70 - AF50, levels = design)
fit                  <- glmFit(miRNA.aging.data, design)
model.result.1       <- glmLRT(fit, contrast = contrast.matrix[,1])
model.result.2       <- glmLRT(fit, contrast = contrast.matrix[,2])
ageing.miRNA.edgeR.pvalue                <- model.result.1$table$PValue
ageing.miRNA.edgeR.sig.table.1           <- model.result.1$table[which(ageing.miRNA.edgeR.pvalue  < 0.05),] %>%
                                            arrange(PValue)
rownames(ageing.miRNA.edgeR.sig.table.1) <- rownames(model.result.1$table)[which(ageing.miRNA.edgeR.pvalue  < 0.05)]

ageing.miRNA.edgeR.pvalue                <- model.result.2$table$PValue
ageing.miRNA.edgeR.sig.table.2           <- model.result.2$table[which(ageing.miRNA.edgeR.pvalue  < 0.05),] %>%
                                            arrange(PValue)
rownames(ageing.miRNA.edgeR.sig.table.2) <- rownames(model.result.2$table)[which(ageing.miRNA.edgeR.pvalue  < 0.05)]

aging.names <- intersect(rownames(ageing.miRNA.edgeR.sig.table.1), rownames(ageing.miRNA.edgeR.sig.table.2))

colnames(dge.tmm.miRNA.log) <- c('AF40-1','AF40-2','AF40-3','AF50-1','AF50-2', 'AF50-3',
                                 'AF60-1','AF60-2','AF60-3','AF70-1','AF70-2', 'AF70-3') 

aging.specific.miRNA        <- dge.tmm.miRNA.log[aging.names,]

color.bar <- colorRampPalette(c("midnightblue", "grey", "mediumvioletred"))(100)
pheatmap( aging.specific.miRNA, cluster_rows = TRUE, cluster_cols = FALSE,
          clustering_distance_rows = 'correlation', color = color.bar,
          scale = 'row', fontsize_col = 10, cellwidth = 23, cellheight = 23)


# QC data aging
# GSE43556
#
# Boon RA, Iekushi K, Lechner S, Seeger T et al. 
# MicroRNA-34a regulates cardiac ageing and function. 
# Nature 2013 Mar 7;495(7439):107-10. PMID: 23426265
#---
setwd('/home/zhenyisong/biodata/cardiodata/GSE43556')
gse       <- getGEO(filename = "GSE43556_family.soft.gz")
gsmlist   <- GSMList(gse)
gpl       <- GPLList(gse)

probesets   <- Table(GPLList(gse)[[2]])$ID
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
                                      {tab <- Table(x)
                                       mymatch <- match(probesets,tab$ID_REF)
                                       return(tab$VALUE[mymatch])
                                     }))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
boxplot(data.matrix)
anno.df     <- dataTable(gpl$GPL7732)@table
gene.ind    <- match(probesets, anno.df$ID)
gene.symbol <- anno.df$'miRNA_ID'[gene.ind]
gene.symbol <- make.names(gene.symbol, unique = TRUE)
rownames(data.matrix) <- NULL
rownames(data.matrix) <- gene.symbol

aging.miRNA.pub   <- na.omit(data.matrix[,-c(1:8)]) %>% log
group             <- factor(rep(1:2, each = 4), levels = 1:2, labels = c('Old','Young'))
design            <- model.matrix(~ 0 + group)
colnames(design)  <- levels(group)
contrast.matrix   <- makeContrasts( Old - Young, levels = design)
fit               <- lmFit(aging.miRNA.pub, design)
fit2              <- contrasts.fit(fit,contrast.matrix)
fit2.miRNA        <- eBayes(fit2)

miRNA.result       <- topTable(  fit2.miRNA, 
                                number        = Inf, 
                                adjust.method = "BH", 
                                sort.by       = "p",
                                p             = 0.05);

miRNA.result['mmu.miR.34a',]
dim(miRNA.result)

hsa.miRNA.aging.matchnames   <- sub("mmu\\.","hsa-", rownames(miRNA.result))
hsa.miRNA.aging.matchnames   <- gsub("\\.","-",hsa.miRNA.aging.matchnames)
hsa.miRNA.aging.gsub         <- gsub("(hsa-miR-\\w+)-.*","\\1",hsa.miRNA.aging.matchnames,perl = TRUE)

# next is GSEA enrichment analysis

ageing.miRNA.edgeR.1.exprs                <- model.result.1$table$logFC
names(ageing.miRNA.edgeR.1.exprs)         <- seq( 1:length(ageing.miRNA.edgeR.1.exprs ))
ageing.miRNA.edgeR.1.sorted.exprs         <- sort(ageing.miRNA.edgeR.1.exprs, decreasing = TRUE)
aging.index.id                            <- names(ageing.miRNA.edgeR.1.exprs)[rownames(model.result.1$table) %in% hsa.miRNA.aging.matchnames  ]

aging2gene  <- data.frame( diseaseId = unlist(as.character(rep('aging.miRNA',length(aging.index.id )))), 
                                        geneId    = unlist(as.integer(aging.index.id)), check.names = TRUE)

aging.GSEA       <- GSEA(ageing.miRNA.edgeR.1.sorted.exprs, TERM2GENE = aging2gene, maxGSSize = 5000, pvalueCutoff = 1)
gseaplot(aging.GSEA, 'aging.miRNA')

#
# venns's digram
# AF60 ~ AF40
# AF70 ~ AF50
# Old ~ Young Mouse
#---
homo.AF60.AF40 <- gsub("(hsa-miR-\\w+)-.*","\\1", rownames(ageing.miRNA.edgeR.sig.table.1), perl = TRUE)
homo.AF70.AF50 <- gsub("(hsa-miR-\\w+)-.*","\\1", rownames(ageing.miRNA.edgeR.sig.table.2), perl = TRUE)
mouse.aging    <- hsa.miRNA.aging.gsub 

area.1           <- length(homo.AF60.AF40)
area.2           <- length(homo.AF70.AF50)
area.3           <- length(mouse.aging)
n.12             <- length( intersect(homo.AF60.AF40,homo.AF70.AF50) )
n.13             <- length( intersect(homo.AF60.AF40,mouse.aging) )
n.23             <- length( intersect(homo.AF70.AF50,mouse.aging) )
n.123            <- length( intersect(intersect(homo.AF70.AF50,mouse.aging), intersect(homo.AF60.AF40, mouse.aging)) )
dev.off()
draw.triple.venn( area.1, area.2, area.3, n.12, n.23, n.13, n.123,
                  alpha = rep(0.5, 3), fill = c('cyan','limegreen','navyblue'),
                  col = rep('gray97',3), lwd = rep(0,3), 
                  category = c('homo.AF60.AF40','homo.AF70.AF50','mouse.aging'),
                  cat.cex = 1.5)


# please check the Nature figure 1.D
# PMID: 23426265
# QC
#---

exprs.df        <- data.frame( avg   = apply(aging.miRNA.pub, 1, median), 
                               ratio = apply(aging.miRNA.pub[,1:4], 1, median) / apply(aging.miRNA.pub[,5:8], 1, median))
exprs.df['mmu.miR.34a',]
ggplot(data = exprs.df, aes( x= avg, y = ratio)) +
geom_point(alpha = 0.02) +
geom_point(aes(x =  exprs.df['mmu.miR.34a', 1] , y = exprs.df['mmu.miR.34a',2]), color = 'red')


#----
# QC miRNA pipeline
# RNA
# data was download from ENCODE
# please check yaoyan.miRNA.sh
# this preocessed data was save 
# as yaoyan.encode.Rdata
# for tissue specific selection
# 
#---

setwd('/home/zhenyisong/biodata/cardiodata/encodeHeart/humanHeart/miRNA/miRNA_bwa')

miRNA.count.files.encode  <- list.files(pattern = "ENC.*\\.txt$")
miRNA.count.matrix.encode <- NULL
column.encode             <- NULL
for( filename in miRNA.count.files.encode) {
    column.encode <- read.delim(filename, header = F, row.names = 1)
    miRNA.count.matrix.encode  <- cbind(miRNA.count.matrix.encode , column.encode$V2 )
}

rownames(miRNA.count.matrix.encode) <- rownames(column.encode)
#colnames(miRNA.count.matrix) <- sub("\\.txt","",miRNA.count.files)
miRNA.count.files.encode
colnames(miRNA.count.matrix.encode) <- c('parietal_lobe','temporal_lobe','spinal_cord','tongue','occipital_lobe',
                                         'spinal_cord','diencephalon','heart','skin_of_body','skeletal_muscle','temporal_lobe','uterus',
                                         'occipital_lobe','metanephros','parietal_lobe','lung','skeletal_muscle','skin_of_body',
                                         'heart','cerebellum', 'frontal_cortex','uterus','cerebellum','frontal_cortex',
                                         'lung','thyroid_gland','diencephalon','thyroid_gland','metanephros',
                                         'tongue','urinary_bladder','liver')

gene.exprs.miRNA.encode    <- DGEList(counts = miRNA.count.matrix.encode[,c(1:5,7:10,12,16,21,26,31,32)]) %>% calcNormFactors()
dge.tmm.miRNA.encode       <- t(t(gene.exprs.miRNA.encode$counts) * gene.exprs.miRNA.encode$samples$norm.factors) %>%
                              apply(2, as.integer)
miRNA.exprs.log.encode           <- log(dge.tmm.miRNA.encode + 1)
rownames(miRNA.exprs.log.encode) <- rownames(column.encode)
miRNA.exprs.sd.encode            <- apply(miRNA.exprs.log.encode, 1, sd)
hist(miRNA.exprs.sd.encode)
estimator.sd.encode                <- shorth(miRNA.exprs.sd.encode)
miRNA.exprs.filtered.encode        <- miRNA.exprs.log.encode[miRNA.exprs.sd.encode > 0.7,]

heatmap.result <- heatmap.2( miRNA.exprs.filtered.encode  , col = greenred(75), scale  = 'row', 
						         Rowv = TRUE,Colv = FALSE, density.info = 'none',
                                 key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
                                 cexCol = 0.8,srtCol = 30,
                                 distfun = function(d) as.dist(1 - cor(t(d), method = 'spearman')),
						         hclustfun  = function(d) hclust(d, method = 'complete'),
						         dendrogram = 'no',labRow = NA);


#---
# index of heart is 7
#---
tissue.max <- apply(miRNA.exprs.filtered.encode, 1, which.max)
tissue.max[tissue.max != 7] <- 0
tissue.max[tissue.max == 7] <- 1

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

tao.index    <- apply(miRNA.exprs.filtered.encode, 1, tao.func)
tissue.score <- tao.index * tissue.max


miRNA.heart.names <- rownames(miRNA.exprs.filtered.encode)[tissue.score > 0.8]
miRNA.heart.exprs <- miRNA.exprs.filtered.encode[tissue.score > 0.8,]

summary(tissue.score)

heatmap.result <- heatmap.2( miRNA.heart.exprs , col = greenred(75), scale  = 'row', 
						         Rowv = TRUE,Colv = FALSE, density.info = 'none',
                                 key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
                                 cexCol = 0.8,srtCol = 30,
                                 distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						         hclustfun  = function(d) hclust(d, method = 'complete'),
						         dendrogram = 'no',labRow = miRNA.heart.names );
#
# this script is from 
# http://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps
#---
draw_colnames_45 <- function (coln, gaps, ...) {
    coord <- pheatmap:::find_coordinates( length(coln), gaps)
    x     <- coord$coord - 0.5 * coord$size
    res   <- textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
    return(res)
}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x = "draw_colnames", value = "draw_colnames_45",
                  ns = asNamespace("pheatmap"))

color.bar <- colorRampPalette(c("green4", "yellow", "red"))(100)
pheatmap( miRNA.heart.exprs, cluster_rows = TRUE, cluster_cols = FALSE,
          clustering_distance_rows = 'correlation', color = color.bar, 
          scale = 'none', fontsize_col = 15)

setwd('/home/zhenyisong/biodata/wanglab/wangdata/Rdata')
save.image('yaoyan.encode.Rdata')
setwd('D:\\wangli_data\\Rdata')
save.image('yaoyan.miRNA.Rdata')
q('no')
