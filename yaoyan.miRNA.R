# @author Yisong zhen
# @since 2017-03-30
# @updated
# @parent
#    yaoyan.miRNA.sh

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


# read
setwd('/home/zhenyisong/biodata/wanglab/wangdata/yaoyan/rawdata/miRNA_bwa')
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
AF.miRNA.trans.data <- DGEList(counts = miRNA.count.matrix[,7:18], group = groups )
AF.miRNA.trans.data <- calcNormFactors(AF.miRNA.trans.data, method = 'TMM')
AF.miRNA.trans.data <- estimateCommonDisp(AF.miRNA.trans.data)
AF.miRNA.trans.data <- estimateTagwiseDisp(AF.miRNA.trans.data,trend = 'movingave')
# I do not know why it falied
#AF.miRNA.trans.data %>% calcNormFactors(method = 'TMM') %>% 
#                   estimateCommonDisp() %>% estimateTagwiseDisp(trend = "movingave")
# the result cannot be accepte by the next step

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
AF.miRNA.edgeR.sig.table     <- model.result$table[which(AF.miRNA.edgeR.pvalue < 0.05),]