library(affy)
library(annotate)
library(oligo)
library(mogene11sttranscriptcluster.db)
library(pd.mogene.1.0.st.v1)
library(hgu133a2.db)
library(limma)
library(xlsx)
library(GEOquery)

# Series GSE52317
# [MoGene-1_0-st] Affymetrix Mouse Gene 1.0 ST Array 
# Platform Design Info for Affymetrix MoGene-1_0-st-v1
# pd.mogene.1.0.st.v1

Ino80.mouse.family <- c('Ino80b','Ino80c','Ino80d','Ino80e','Ino80')
Ino80.human.family <- c('INO80B','INO80C','INO80D','INO80E','INO80',
                        'TFPT','NFRKB','RUVBL1','RUVBL2','ACTL6A',
                        'ACTR5','ACTR8','MCRS1','UCHC5','YY1')
Ino80.human.illumina.family <- c('ZNHIT4','C18ORF37','FLJ20309','CCDC95','INOC1',
                                  'TFPT','NFRKB','RUVBL1','RUVBL2','ACTL6A',
                                  'ACTR5','ACTR8','MCRS1','UCHC5','YY1')
#---
# mouse dataset I
#---
setwd('/home/zhenyisong/wanglab/zongna/GSE52317')
celPath  <- list.files( getwd(), pattern = "\\.CEL|\\.CEL\\.gz", 
                      full.names = TRUE, ignore.case = TRUE)
raw.data    <- read.celfiles(celPath, pkgname = "pd.mogene.1.0.st.v1")
rma.data    <- rma(raw.data,target = 'core')
exprs.data  <- exprs(rma.data)
annodb      <- "mogene11sttranscriptcluster.db"
probe.id    <- featureNames(rma.data)
gene.symbol <- as.character(lookUp(probe.id, annodb, "SYMBOL"))
# Name   <- as.character(lookUp(ID, annodb, "GENENAME"))
# Entrez <- as.character(lookUp(ID, annodb, "ENTREZID"))
# the above code is from 
# https://www.biostars.org/p/69597/
#---
rownames(exprs.data)  <- gene.symbol
f                     <- factor( c(1,1,2,2,3,3,3,4,4,5,5,5), levels = c(1:5),
                                 labels = c("WT","GATA4.WT","GATA4.KO","GATA6.WT","GATA6.KO"))
design                <- model.matrix(~ 0 + f)
colnames(design)      <- c("WT","GATA4.WT","GATA4.KO","GATA6.WT","GATA6.KO")
fit                   <- lmFit(exprs.data, design)
contrast.matrix       <- makeContrasts(GATA4.KO - GATA4.WT, GATA6.KO - GATA6.WT,levels = design)
fit2                  <- contrasts.fit(fit, contrast.matrix)
fit2                  <- eBayes(fit2)
result                <- topTable(fit2, coef = 1, adjust = "BH", number = Inf)

GATA4                 <- result[result$ID %in% Ino80.mouse.family,]
result                <- topTable(fit2, coef = 2, adjust = "BH", number = Inf)
GATA6                 <- result[result$ID %in% Ino80.mouse.family,]

setwd('/home/zhenyisong/wanglab/zongna/')
write.xlsx(GATA4, 'GSE52317.xlsx', sheetName  = "Sheet1")
write.xlsx(GATA4, 'GSE52317.xlsx', sheetName  = "Sheet2", append = TRUE)

#
# mouse dataset II
#---
setwd('/home/zhenyisong/wanglab/zongna/GSE30428')
celPath  <- list.files( getwd(), pattern = "\\.CEL|\\.CEL\\.gz", 
                      full.names = TRUE, ignore.case = TRUE)
raw.data    <- read.celfiles(celPath, pkgname = "pd.mogene.1.0.st.v1")
rma.data    <- rma(raw.data,target = 'core')
exprs.data  <- exprs(rma.data)
annodb      <- "mogene11sttranscriptcluster.db"
probe.id    <- featureNames(rma.data)
gene.symbol <- as.character(lookUp(probe.id, annodb, "SYMBOL"))

rownames(exprs.data)  <- gene.symbol
f                     <- factor( c(1,1,2,2,3,3,4,4,5,5,6,6), levels = c(1:6),
                                 labels = c("transaortic","sham_left","pulmonary_1w","pulmonary_3w",
                                            "pulmonary_6w",'sham_right' ))
design                <- model.matrix(~ 0 + f)
colnames(design)      <- c( "transaortic","sham_left",
                            "pulmonary_1w","pulmonary_3w",
                            "pulmonary_6w",'sham_right' )
fit                   <- lmFit(exprs.data, design)
contrast.matrix       <- makeContrasts(transaortic - sham_left,levels = design)
fit2                  <- contrasts.fit(fit, contrast.matrix)
fit2                  <- eBayes(fit2)
result                <- topTable(fit2, coef = 1, adjust = "BH", number = Inf)
hcm.left              <- result[result$ID %in% Ino80.mouse.family,]
setwd('/home/zhenyisong/wanglab/zongna/')
write.xlsx(hcm.left, 'GSE30428.xlsx', sheetName  = "Sheet1")

#
# human dataset 
#---
setwd('/home/zhenyisong/wanglab/zongna/GSE17800')
raw.data    <- ReadAffy();
rma.data    <- affy::rma(raw.data);
exprs.data  <- exprs(rma.data)
probes      <- rownames(exprs.data)
gene.symbol <- unlist(mget(probes, hgu133a2SYMBOL, ifnotfound = NA))
rownames(exprs.data) <- gene.symbol
f           <- factor( c( rep(1,40),rep(2,8) ), levels = c(1:2),
                       labels = c("patient","normal") )
design      <- model.matrix(~ 0 + f)
colnames(design)   <- c( "patient","normal")
fit                   <- lmFit(exprs.data, design)
contrast.matrix       <- makeContrasts(patient - normal,levels = design)
fit2                  <- contrasts.fit(fit, contrast.matrix)
fit2                  <- eBayes(fit2)
result                <- topTable(fit2, coef = 1, adjust = "BH", number = Inf)
dcm.left              <- result[result$ID %in% Ino80.human.family,]
setwd('/home/zhenyisong/wanglab/zongna/')
write.xlsx(dcm.left, 'GSE17800.xlsx', sheetName  = "Sheet1")

#
# human dataset
#---
setwd('/home/zhenyisong/wanglab/zongna/GSE60291')
raw.data    <- ReadAffy();
rma.data    <- affy::rma(raw.data);
exprs.data  <- exprs(rma.data)
probes      <- rownames(exprs.data)
gene.symbol <- unlist(mget(probes, hgu133a2SYMBOL, ifnotfound = NA))
rownames(exprs.data) <- gene.symbol
f           <- factor( c(1,2,1,2,1,2 ), levels = c(1:2),
                       labels = c("Control","ET.1"))
design      <- model.matrix(~ 0 + f)
colnames(design)   <- c( "Control","ET.1")
fit                   <- lmFit(exprs.data, design)
contrast.matrix       <- makeContrasts(Control - ET.1,levels = design)
fit2                  <- contrasts.fit(fit, contrast.matrix)
fit2                  <- eBayes(fit2)
result                <- topTable(fit2, coef = 1, adjust = "BH", number = Inf)
et1                   <- result[result$ID %in% Ino80.human.family,]
setwd('/home/zhenyisong/wanglab/zongna/')
write.xlsx(et1, 'GSE60291.xlsx', sheetName  = "Sheet1")

#
# human dataset
# Please see Using the GEOquery Package
# By Sean Davis
#---
setwd('/home/zhenyisong/wanglab/zongna/GSE32453')
gse <- getGEO(filename = "GSE32453_family.soft.gz")
gsmlist <- GSMList(gse)
gpl <- GPLList(gse)

probesets   <- Table(GPLList(gse)[[1]])$ID
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
                                      {tab <- Table(x)
                                       mymatch <- match(probesets,tab$ID_REF)
                                       return(tab$VALUE[mymatch])
                                     }))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
data.matrix <- log2(data.matrix)
rownames(data.matrix) <- probesets
colnames(data.matrix) <- names(gsmlist)
pdata <- data.frame(samples = names(gsmlist))
rownames(pdata) <- names(gsmlist)
pheno <- as(pdata,"AnnotatedDataFrame")
eset2 <- new('ExpressionSet',exprs = data.matrix,phenoData = pheno)
f           <- factor( c( rep(1,8),rep(2,5) ), levels = c(1:2),
                       labels = c("HCM","normal") )
design      <- model.matrix(~ 0 + f)
colnames(design)   <- c( "HCM","normal")
fit                   <- lmFit(eset2, design)
contrast.matrix       <- makeContrasts(HCM - normal,levels = design)
fit2                  <- contrasts.fit(fit, contrast.matrix)
fit2                  <- eBayes(fit2)
result                <- topTable(fit2, coef = 1, adjust = "BH", number = Inf)

dataTable(gpl$GPL6104)@table
# names(dataTable(gpl$GPL6104)@table)
anno.df     <- dataTable(gpl$GPL6104)@table
gene.ind    <- match(rownames(result), anno.df$ID)
gene.symbol <- anno.df$Symbol[gene.ind]
result$gene.symbol <- gene.symbol
et1                <- result[gene.symbol %in% Ino80.human.illumina.family,]
setwd('/home/zhenyisong/wanglab/zongna/')
write.xlsx(et1, 'GSE32453.xlsx', sheetName  = "Sheet1")

# GSE9800
# Human dataset
#---

setwd('/home/zhenyisong/wanglab/zongna/GSE9800')
gse <- getGEO(filename = "GSE9800_family.soft.gz")
gsmlist <- GSMList(gse)
gpl <- GPLList(gse)

probesets   <- Table(GPLList(gse)[[1]])$ID
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
                                      {tab <- Table(x)
                                       mymatch <- match(probesets,tab$ID_REF)
                                       return(tab$VALUE[mymatch])
                                     }))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
data.matrix <- log2(data.matrix)
rownames(data.matrix) <- probesets
colnames(data.matrix) <- names(gsmlist)
pdata <- data.frame(samples = names(gsmlist))
rownames(pdata) <- names(gsmlist)
pheno <- as(pdata,"AnnotatedDataFrame")
eset2 <- new('ExpressionSet',exprs = data.matrix,phenoData = pheno)
f           <- factor( c( rep(1,5),2, 2, c(3:6), 2, 7, 2, 2, rep(1, 3), 5, 8,rep(2,6), 1, 1, 2, 1), levels = c(1:8),
                       labels = c("normal","DCM",'myocarditis','sarcoidosis','ischemic','peripartal','alcoholic','HCM') )
f.2         <- factor( c( rep(1,5),2, 2, rep(3,4), 2, 3, 2, 2, rep(1, 3), 3, 3,rep(2,6), 1, 1, 2, 1), levels = c(1:3),
                       labels = c("normal","DCM",'other') )

design      <- model.matrix(~ 0 + f.2)
#colnames(design)   <- c("normal","DCM",'myocarditis','sarcoidosis','ischemic','peripartal','alcoholic','HCM')
colnames(design)   <- c("normal","DCM",'Other')
fit                   <- lmFit(eset2, design)
#contrast.matrix       <- makeContrasts(HCM - normal, DCM - normal,levels = design)
contrast.matrix       <- makeContrasts(Other - normal, DCM - normal,levels = design)
fit2                  <- contrasts.fit(fit, contrast.matrix)
fit2                  <- eBayes(fit2)
result1               <- topTable(fit2, coef = 1, adjust = "BH", number = Inf)
result2               <- topTable(fit2, coef = 2, adjust = "BH", number = Inf)
# dataTable(gpl$GPL887)@table
# names(dataTable(gpl$GPL887)@table)
anno.df     <- dataTable(gpl$GPL887)@table
gene.ind    <- match(rownames(result1), anno.df$ID)
gene.symbol <- anno.df$GENE_SYMBOL[gene.ind]
result1$gene.symbol <- gene.symbol
result2$gene.symbol <- gene.symbol
et1                <- result1[gene.symbol %in% Ino80.human.illumina.family,]
et2                <- result2[gene.symbol %in% Ino80.human.illumina.family,]
setwd('/home/zhenyisong/wanglab/zongna/')
write.xlsx(et1, 'GSE9800.xlsx', sheetName  = "Other")
write.xlsx(et2, 'GSE9800.xlsx', sheetName  = "DCM",append = TRUE)