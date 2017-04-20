# @data source
#     deposited at Fuwai SuperComputer
#     /home/zhenyisong/data/wanglilab/projects/yaoyan/rawdata
# @author Yisong Zhen
# @since 2017-03-08
# cd /home/zhenyisong/data/wanglilab/wangcode
# nohup R CMD BATCH yaoyan.R &
#---
library(gplots)
library(xlsx)
library(Rsubread)
library(edgeR)
library(limma)
library(DESeq2)
library(genefilter)
library(RColorBrewer)
library(gridExtra)
library(affy)
library(annotate)
library(org.Hs.eg.db)
library(VennDiagram)
library(Rsamtools)
library(clusterProfiler)
library(pathview)
library(tidyverse)
library(biomaRt)
library(cluster)
library(factoextra)
library(GEOquery)


setwd('/mnt/date/biodata/wanglab/wangdata/yaoyan/rawdata/rsubread')
setwd('D:\\wangli_data\\Rdata')
load("yaoyan.Rdata")

genome_ref.path      <- "/mnt/date/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
setwd('/home/zhenyisong/biodata/wanglab/wangdata/yaoyan/rawdata')
unlink('rsubread', force = TRUE, recursive = TRUE)
dir.create('rsubread')
output.path          <- '/home/zhenyisong/biodata/wanglab/wangdata/yaoyan/rawdata/rsubread'
setwd('/home/zhenyisong/biodata/wanglab/wangdata/yaoyan/rawdata/rsubread')
yaoyan.rawdata       <- list.files( path = "../",pattern = "\\.fq.\\gz$", all.files = FALSE,
                                    full.names = TRUE, recursive = TRUE, 
                                    ignore.case = FALSE, include.dirs = TRUE)
base.string          <-  'hg19_index'
buildindex( basename = base.string, reference = genome_ref.path )
mRNA.yaoyan.rawdata  <- yaoyan.rawdata[1:60]
reads.paths.1        <- mRNA.yaoyan.rawdata[grep("\\.1\\.clean\\.fq\\.gz",mRNA.yaoyan.rawdata)]
reads.paths.2        <- mRNA.yaoyan.rawdata[grep("\\.2\\.clean\\.fq\\.gz",mRNA.yaoyan.rawdata)]
read.file.names      <- basename(reads.paths.1)
read.file.names      <- sub("\\.\\d\\D+$", "\\1", read.file.names, perl = TRUE)
outputs.files        <- paste0(output.path,'/', read.file.names,'.bam')
align( index          = base.string, 
       readfile1      = reads.paths.1, 
       readfile2      = reads.paths.2, 
       input_format   = "gzFASTQ", 
       type           = 'rna', 
       output_file    = outputs.files, 
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 15, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

# merge 2 bam file into one
# file needed not to be sorted, if sorted, it will re-order it.
# I tested the old way to sort and merge, but this protocol is immature.
#-----
file.needed.sort <- outputs.files[7:30]
file.odd.bam     <- file.needed.sort[1:length(file.needed.sort) %% 2 == 1]
file.even.bam    <- file.needed.sort[1:length(file.needed.sort) %% 2 == 0]
merged.bams      <- c()
for( i in c(1:length(file.odd.bam)) ) {
   odd.name  <- basename(file.odd.bam[i])
   even.name <- basename(file.even.bam[i])
   file.name <- sub("_lane.*$", "\\1", odd.name, perl = TRUE)
   # old rascal, deprecated!
   #file.odd.temp  <- paste0(file.name,"_temp_1")
   #file.even.temp <- paste0(file.name,"_temp_2")
   #bam1 <- sortBam(odd.name, file.odd.temp);
   #bam2 <- sortBam(even.name, file.even.temp);
   file.name <- paste0(file.name,'.bam')
   #bam3 <- mergeBam(c(bam1,bam2), file.name, overwrite = TRUE)
   bam3 <- mergeBam(c(odd.name,even.name), file.name, overwrite = TRUE)
   merged.bams <- c(merged.bams,bam3)
}

unlink("*_temp_*\\.bam")

#save.image(file = 'yaoyan.Rdata')
#quit('no')

final.bams <- c(basename(outputs.files)[1:6], merged.bams)


yaoyan.mRNA.genes      <- featureCounts( final.bams, useMetaFeatures = TRUE,
                                         countMultiMappingReads = FALSE,
                                         strandSpecific         = 0, 
                                         isPairedEnd            = TRUE,
                                         requireBothEndsMapped  = TRUE,
                                         autosort               = TRUE,
                                         nthreads               = 15,
                                         annot.inbuilt = "hg19", allowMultiOverlap = TRUE)

#save.image(file = 'yaoyan.Rdata')
#quit('no')
gene         <- yaoyan.mRNA.genes
gene.counts  <- gene$counts[,7:18]
gene.ids     <- gene$annotation$GeneID

columns      <- c("ENTREZID","SYMBOL", "GENENAME");
GeneInfo     <- AnnotationDbi::select( org.Hs.eg.db, keys= as.character(gene.ids), 
                       keytype = "ENTREZID", columns = columns);
m            <- match(gene$annotation$GeneID, GeneInfo$ENTREZID);
Ann          <- cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                              GeneInfo[m, c("SYMBOL", "GENENAME")]);
             
Ann$Chr      <- unlist( lapply(strsplit(Ann$Chr, ";"), 
                        function(x) paste(unique(x), collapse = "|")))
Ann$Chr      <- gsub("chr", "", Ann$Chr)

gene.exprs   <- DGEList(counts = gene.counts, genes = Ann)
gene.exprs   <- calcNormFactors(gene.exprs)
dge.tmm      <- t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
dge.tmm.counts           <- apply(dge.tmm,2, as.integer)

sample.info              <- data.frame( treat  = c('Age60','Age60','Age60',
                                                   'Age60_AF','Age60_AF','Age60_AF',
                                                   'Age70','Age70','Age70',
                                                   'Age70_AF','Age70_AF','Age70_AF') )
dds                      <- DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                                    colData   = sample.info,
                                                    design    = ~ treat)
vsd                      <- varianceStabilizingTransformation(dds);
vsd.expr                 <- assay(vsd)

sds <- rowSds(vsd.expr)
sh  <- shorth(sds)
vsd.filtered.expr <- vsd.expr[sds > 0.3,]

heatmap(cor(vsd.filtered.expr, method = 'spearman'), cexCol = 0.8, cexRow = 0.8)

age.param        <- factor(c(rep(1,6),rep(2,6)),levels = 1:2 , labels = c('age_60','age_70'))
patient.param    <- factor(c(rep(1,3),rep(2,3),rep(1,3),rep(2,3)),levels = 1:2 , labels = c('Normal','Patient'))
design           <- model.matrix(~ 0 + patient.param + age.param);
colnames(design) <- c('Normal',"Patient",'Age')
contrast.matrix  <- makeContrasts(Patient - Normal, levels = design)
d.norm           <- voom(gene.exprs, design = design)
fit              <- lmFit(d.norm, design)
fit2             <- contrasts.fit(fit,contrast.matrix)
fit2             <- eBayes(fit2)
gene.result      <- topTable(  fit2, 
                               number        = Inf, 
                               adjust.method = "BH", 
                               sort.by       = "p"
                               );


# post_QC
# sex gene selection and 
# https://www.biostars.org/p/3650/
# Question: How To Get A List Of All The Genes On The Human Chromosome Y
# PMID: 25376837 
# Sex-biased gene expression and sexual conflict throughout development.
# ---
sample.names <- c( 'SR60-1', 'SR60-2', 'SR60-3', 'SR70-1','SR70-2', 'SR70-3',
                   'AF60-1', 'AF60-2', 'AF60-3', 'AF70-1', 'AF70-2', 'AF70-3')
sample.sex   <- factor(c(2,1,1,1,1,2,2,1,1,1,2,2), levels = 1:2, labels = c('Male','Female'))

mart    <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM( attributes = c("chromosome_name", "entrezgene", "hgnc_symbol"),
                  filters = "chromosome_name", values = "Y", mart = mart)
y.chromosome.gene <- unique(results$hgnc_symbol) %>% na.omit()
y.chromosome.gene <- y.chromosome.gene[ y.chromosome.gene != '']
vsd.expr.sex      <- vsd.expr[ Ann$SYMBOL %in% y.chromosome.gene, c(1:3,6:8,4,12,5,9:11)]
colnames(vsd.expr.sex)
colnames(vsd.expr.sex) <- sample.names
sds.sex <- rowSds(vsd.expr.sex)
summary(sds.sex)
vsd.expr.sex[sds.sex > 0.5,c(1,6,7,11)]
km.res <- kmeans(t(vsd.expr.sex[sds.sex > .7,]), 2, nstart = 100)


fviz_nbclust(t(vsd.expr.sex[sds.sex > .7,]),  kmeans, method = "silhouette")+
theme_classic()

fviz_cluster(km.res, data = t(vsd.expr.sex[sds.sex > .7,]),
             palette = c("#2E9FDF", "#FC4E07"),
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal()
)

# QC sample
# Confirm this sample is from human atrium
# cross-validated with rat heart samples
# public data was downloaded from 
# GEO database
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5266
# 
#---
setwd("D:\\wangli_data\\GSE5266")
gse       <- getGEO(filename = "GSE5266_family.soft.gz")
gsmlist   <- GSMList(gse)
gpl       <- GPLList(gse)

probesets   <- Table(GPLList(gse)[[1]])$ID
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
                                      {tab <- Table(x)
                                       mymatch <- match(probesets,tab$ID_REF)
                                       return(tab$VALUE[mymatch])
                                     }))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
boxplot(data.matrix)
anno.df     <- dataTable(gpl$GPL1355)@table
gene.ind    <- match(probesets, anno.df$ID)
gene.symbol <- anno.df$'Gene Symbol'[gene.ind]
gene.symbol <- make.names(gene.symbol, unique = TRUE)
rownames(data.matrix) <- NULL
rownames(data.matrix) <- gene.symbol
rat.sampel.sd         <- rowSds(data.matrix)
hist(rat.sampel.sd)
sh  <- shorth(rat.sampel.sd, na.rm = T)
rat.pca.matrix        <- data.matrix[rat.sampel.sd > 0.12, ]

#
# GSM119422 -- GSM119425    Rat atrium
# GSM119426 -- GSM119429    Rat ventricle
sum(is.na(rat.pca.matrix))
checke.na <- rowSums(rat.pca.matrix)
rat.pca.matrix <- rat.pca.matrix[! is.na(checke.na), ]
rat.chamber.PCA  <- prcomp(t(rat.pca.matrix))

plot(rat.chamber.PCA$x)
pr.var <- rat.chamber.PCA$sdev^2
pve    <- pr.var/sum(pr.var)

rotations <- order(rat.chamber.PCA$rotation[,1], decreasing = FALSE)
gene.len  <- 20
total.len <- length(rownames(rat.pca.matrix))
atrium.markers    <- unique(rownames(rat.pca.matrix)[rotations[1:gene.len]]) %>% toupper()
ventricle.markers <- unique(rownames(rat.pca.matrix)[rotations[(total.len - gene.len):total.len]] ) %>% toupper()

gene.exprs.log   <- rpkm(gene.exprs, log = T)
genes.atrium     <- gene.exprs.log[Ann$SYMBOL %in% atrium.markers ,]
genes.ventricle  <- gene.exprs.log[Ann$SYMBOL %in% ventricle.markers ,]

#chamber.df       <- data.frame(atrium = )

boxplot(genes.atrium[,1], genes.ventricle[,1])
ks.test(genes.atrium[,1], genes.ventricle[,1])





## test GSEA
# http://www.disgenet.org/ds/DisGeNET/results/all_gene_disease_associations.tsv.gz
# DisGeNET which contains 381056 associations, between 16666 genes and 13172 diseases 
#---

setwd("C:\\Users\\Yisong\\Documents\\Tencent Files\\2837113358\\FileRecv\\all_gene_disease_associations.tsv")
gda <- read.delim('all_gene_disease_associations.tsv', header = TRUE,  comment.char = "#")
disease2gene <- gda[, c("diseaseId", "geneId")]
disease2name <- gda[, c("diseaseId", "diseaseName")]
set.seed(123)
AF.pathway.genelist        <- gene.result$logFC
names(AF.pathway.genelist) <- gene.result$GeneID
AF.pathway.genelist        <- sort(AF.pathway.genelist, decreasing = TRUE)
AF.GSEA <- GSEA(AF.pathway.genelist, TERM2GENE = disease2gene, TERM2NAME = disease2name) 

gsea.result <- summary(AF.GSEA)
gsea.result$Des

# Inon Channel list
setwd("E:\\FuWai\\寿实验室\\ion_channel")
ion.list         <- read.csv('IonChannelsListDownloadForward.csv')
ion.geneset.id   <- ion.list$HGNC.ID
ion.pathway.list <-  data.frame( diseaseId = unlist(rep('Ion.pathway', length(ion.geneset.id ))), 
                                             geneId   = unlist(ion.geneset.id) )
gsea.result      <- GSEA(AF.pathway.genelist ,TERM2GENE = ion.pathway.list, maxGSSize = 1000 )

# 

# negative result
#---end


# explorative analysis
# using log RPKM 
#--

gene.exprs.log   <- rpkm(gene.exprs, log = T)
gene.exprs.median <- NULL
gene.exprs.median <- cbind( apply(gene.exprs.log[,c(1:3,7:9)],1,median),
                            apply(gene.exprs.log[,c(4:6,10:12)],1,median))
rownames(gene.exprs.median) <- make.names(Ann$SYMBOL, unique = TRUE)
gene.exprs.median           <- cbind(gene.exprs.median, rep(0,dim(gene.exprs.median)[1]))
colnames(gene.exprs.median) <- c('Health','Patient','Type')

setwd("E:\\FuWai\\寿实验室\\ion_channel")
ion.list         <- read.csv( 'IonChannelsListDownloadForward.csv', 
                              stringsAsFactors = FALSE)
ion.genename     <- ion.list$Human.gene.name

gene.exprs.median[row.names(gene.exprs.median) %in% ion.genename,3] <- 1
gene.exprs.median         <- as.data.frame(gene.exprs.median)

vsmc <- ggplot(data = gene.exprs.median, aes(x = Health, y = Patient, color = as.factor(Type), alpha = as.factor(Type)) ) +
       geom_point( ) +
       scale_alpha_manual(values = c(1/200, 1), guide = FALSE) + 
       scale_colour_manual(name = 'gene groups',values = c("black","red"), 
                           labels = c('non-related genes','inon gene family')) + 

       geom_abline(intercept = 0.58, slope = 1, size = 1, alpha = 1/5) +
       geom_abline(intercept = -0.58, slope = 1, size = 1, alpha = 1/5) +
       theme(legend.position = c(0.8, 0.2),legend.title.align = 0.5)

theme_4ppt <- theme()

gene.exprs[row.names(gene.exprs.median) %in% ion.genename,]
gene.exprs.log[row.names(gene.exprs.median) %in% ion.genename,]
ion.gene.fc <- gene.exprs.median[(row.names(gene.exprs.median) %in% ion.genename) &
                     apply(as.matrix(gene.exprs.median), 1 , function(x) abs(x[1] - x[2])) > 0.58,]
dim(ion.gene.fc)[1]
length(ion.genename)

over_used.inogenes <- rownames(ion.gene.fc)

ion.gene.list      <- apply(gene.exprs.median[(row.names(gene.exprs.median) %in% ion.genename),], 1 , function(x) x[1] - x[2])
ion.gene.list      <- sort(ion.gene.list, decreasing = TRUE)

# subgroup of Inwardly rectifying potassium channels
# this list is extract manually from 

subgroup.1.genes <- c( 'KCNJ1','KCNJ2','KCNJ12', 'KCNJ4', 'KCNJ14','KCNJ3',
                       'KCNJ6','KCNJ9','KCNJ5','KCNJ10','KCNJ15','KCNJ16',
                       'KCNJ8','KCNJ11','KCNJ13')

# Calcium-activated potassium channels
subgroup.2.genes <- c('KCNMA1', 'KCNN1', 'KCNN2','KCNN3','KCNN4', 'KCNT1','KCNT2','KCNU1')

# Two P domain potassium channels
subgroup.3.genes <- c( 'KCNK1', 'KCNK2', 'KCNK3', 'KCNK4', 'KCNK5', 'KCNK6',
                       'KCNK7','KCNK9','KCNK10','KCNK12','KCNK13','KCNK15',
                       'KCNK16', 'KCNK17', 'KCNK18')

# Voltage-gated potassium channels
subgroup.4.genes <- c( 'KCNA1', 'KCNA2','KCNA3', 'KCNA4', 'KCNA5','KCNA6','KCNA7',
                       'KCNA10', 'KCNB1', 'KCNB2', 'KCNC1','KCNC2', 'KCNC3', 'KCNC4',
                       'KCND1', 'KCND2', 'KCND3', 'KCNF1', 'KCNG1', 'KCNG2', 'KCNG3',
                       'KCNG4', 'KCNQ1', 'KCNQ2', 'KCNQ3', 'KCNQ4', 'KCNQ5','KCNV1',
                       'KCNV2', 'KCNS1', 'KCNS2', 'KCNS3', 'KCNH1', 'KCNH5', 'KCNH2',
                       'KCNH6', 'KCNH7', 'KCNH8', 'KCNH3', 'KCNH4')

setwd('D:\\wangli_data\\Rdata')
save.image(file = "yaoyan.Rdata")
quit("no")




