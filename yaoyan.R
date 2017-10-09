# @data source
#     deposited at Fuwai SuperComputer
#     /home/zhenyisong/data/wanglilab/projects/yaoyan/rawdata
# @author Yisong Zhen
# @since 2017-03-08
# cd /home/zhenyisong/data/wanglilab/wangcode
# nohup R CMD BATCH yaoyan.R &
#---

# supplement
# http://bioconductor.org/help/workflows/rnaseqGene/#time-course-experiments
#---

# GSEA 
# ClusterProfiler
# 
# http://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
# https://guangchuangyu.github.io/clusterProfiler/
# ---
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
library(hgu133plus2.db)
library(hgu133a.db)
library(hgu133b.db)


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



# QC motivation
# atrium ~ ventricle gene distribution
# 
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
gene.len  <- 400
total.len <- length(rownames(rat.pca.matrix))
atrium.markers    <- unique(rownames(rat.pca.matrix)[rotations[1:gene.len]]) %>% toupper()
ventricle.markers <- unique(rownames(rat.pca.matrix)[rotations[(total.len - gene.len):total.len]] ) %>% toupper()

# this is to check whether the distribution is 

rat.atrium.exprs.a <- data.matrix[gene.symbol %in% atrium.markers, 1:4] %>% apply(1, median)
rat.atrium.exprs.v <- data.matrix[gene.symbol %in% ventricle.markers, 1:4] %>% apply(1, median)

boxplot(rat.atrium.exprs.a, rat.atrium.exprs.v, names = c('atrium', 'ventricle'))
ks.test(rat.atrium.exprs.a, rat.atrium.exprs.v)

## > ks.test(rat.atrium.exprs.a, rat.atrium.exprs.v)
## 
## 	Two-sample Kolmogorov-Smirnov test
## 
## data:  rat.atrium.exprs.a and rat.atrium.exprs.v
## D = 1, p-value = 0.001099
## alternative hypothesis: two-sided
## 

boxplot(rat.atrium.exprs.a, rat.atrium.exprs.v, xlab = c('atrium', 'ventricle'))

rat.ventricle.exprs.a <- data.matrix[gene.symbol %in% atrium.markers, 5:8]
rat.ventricle.exprs.v <- data.matrix[gene.symbol %in% ventricle.markers, 5:8]

ks.test(rat.ventricle.exprs.a, rat.ventricle.exprs.v)

## > ks.test(rat.ventricle.exprs.a, rat.ventricle.exprs.v)
## 
## 	Two-sample Kolmogorov-Smirnov test
## 
## data:  rat.ventricle.exprs.a and rat.ventricle.exprs.v
## D = 0.91667, p-value = 0.005495
## alternative hypothesis: two-sided

boxplot(rat.ventricle.exprs.a, rat.ventricle.exprs.v, xlab = c('atrium', 'ventricle'))



gene.exprs.log   <- rpkm(gene.exprs, log = T)
genes.atrium     <- gene.exprs.log[Ann$SYMBOL %in% atrium.markers, c(1:3, 6:8)] %>% apply(1, median)
genes.ventricle  <- gene.exprs.log[Ann$SYMBOL %in% ventricle.markers, c(1:3, 6:8)]  %>% apply(1,median)

#chamber.df       <- data.frame(atrium = )

boxplot(as.vector(genes.atrium), as.vector(genes.ventricle), names = c('atrium','ventricle'))
ks.test(as.vector(genes.atrium), as.vector(genes.ventricle))

genes.atrium     <- gene.exprs.log[Ann$SYMBOL %in% atrium.markers, -c(1:3,6:8)] %>% apply(1,median)
genes.ventricle  <- gene.exprs.log[Ann$SYMBOL %in% ventricle.markers, -c(1:3,6:8)] %>% apply(1,median)
boxplot(as.vector(genes.atrium), as.vector(genes.ventricle), names = c('atrium', 'venctricle'))
ks.test(as.vector(genes.atrium), as.vector(genes.ventricle))





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

gseaplot(gsea.result, 'Ion.pathway' )

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
names(ion.gene.list) <- Ann$GeneID[Ann$SYMBOL %in% names(ion.gene.list)]

# subgroup of Inwardly rectifying potassium channels
# this list is extract manually from 
# data source
# IUPHAR/BPS database
# http://www.guidetopharmacology.org/
#----


# Potassium channels :Inwardly rectifying potassium channels
subgroup.1.genes <- c( 'KCNJ1','KCNJ2','KCNJ12', 'KCNJ4', 'KCNJ14','KCNJ3',
                       'KCNJ6','KCNJ9','KCNJ5','KCNJ10','KCNJ15','KCNJ16',
                       'KCNJ8','KCNJ11','KCNJ13')

# Potassium channels : Calcium-activated potassium channels
subgroup.2.genes <- c('KCNMA1', 'KCNN1', 'KCNN2','KCNN3','KCNN4', 'KCNT1','KCNT2','KCNU1')

# Potassium channels : Two P domain potassium channels
subgroup.3.genes <- c( 'KCNK1', 'KCNK2', 'KCNK3', 'KCNK4', 'KCNK5', 'KCNK6',
                       'KCNK7','KCNK9','KCNK10','KCNK12','KCNK13','KCNK15',
                       'KCNK16', 'KCNK17', 'KCNK18')

# Potassium channels : Voltage-gated potassium channels
subgroup.4.genes <- c( 'KCNA1', 'KCNA2','KCNA3', 'KCNA4', 'KCNA5','KCNA6','KCNA7',
                       'KCNA10', 'KCNB1', 'KCNB2', 'KCNC1','KCNC2', 'KCNC3', 'KCNC4',
                       'KCND1', 'KCND2', 'KCND3', 'KCNF1', 'KCNG1', 'KCNG2', 'KCNG3',
                       'KCNG4', 'KCNQ1', 'KCNQ2', 'KCNQ3', 'KCNQ4', 'KCNQ5','KCNV1',
                       'KCNV2', 'KCNS1', 'KCNS2', 'KCNS3', 'KCNH1', 'KCNH5', 'KCNH2',
                       'KCNH6', 'KCNH7', 'KCNH8', 'KCNH3', 'KCNH4')


# Voltage-gated ion channels : CatSper and Two-Pore channels

subgroup.5.genes <- c('CATSPER1', 'CATSPER2', 'CATSPER3','CATSPER4', 'TPCN1', 'TPCN2')

# Voltage-gated ion channels : Cyclic nucleotide-regulated channels

subgroup.6.genes <- c('CNGA1','CNGA2', 'CNGA3', 'CNGA4', 'CNGB1', 'CNGB3', 'HCN1','HCN2','HCN3','HCN4')

# Ryanodine receptor

subgroup.7.genes <- c('RYR1','RYR2','RYR3')

# Transient Receptor Potential channels

subgroup.8.genes <- c( 'TRPA1','TRPC1', 'TRPC2', 'TRPC3','TRPC4','TRPC5', 'TRPC6','TRPC7',
                       'TRPM1', 'TRPM2', 'TRPM3', 'TRPM4', 'TRPM5', 'TRPM6', 'TRPM7','TRPM8',
                       'MCOLN1', 'MCOLN2', 'MCOLN3', 'PKD2', 'PKD2L1', 'PKD2L2', 'TRPV1',
                       'TRPV2', 'TRPV3', 'TRPV4','TRPV5','TRPV6' )

# Voltage-gated calcium channels

subgroup.9.genes <- c( 'CACNA1S','CACNA1C','CACNA1D', 'CACNA1F', 'CACNA1A', 'CACNL1A5',
                       'CACNA1E','CACNA1G','CACNA1H', 'CACNA1I')

#  Voltage-gated proton channel
subgroup.10.genes <- c('HVCN1')

# Voltage-gated sodium channels

subgroup.11.genes <- c( 'SCN1A', 'SCN2A', 'SCN3A','SCN4A','SCN5A','SCN8A',
                        'SCN9A', 'SCN10A','SCN11A')


# Ligand-gated ion channels : 5-HT3 receptors
subgroup.12.genes <- c('HTR3A','HTR3B','HTR3C','HTR3D','HTR3E')

# Acid-sensing (proton-gated) ion channels (ASICs)

subgroup.13.genes <- c('ASIC1','ASIC2','ASIC3')

#  Epithelial sodium channels (ENaC)
subgroup.14.genes <- c('SCNN1A','SCNN1B','SCNN1D','SCNN1G')


#  GABAA receptors
subgroup.15.genes  <- c( 'GABRA1','GABRA2','GABRA3','GABRA4','GABRA5','GABRA6',
                         'GABRB1', 'GABRB2', 'GABRB3','GABRG1','GABRG2','GABRG3',
                         'GABRD','GABRE','GABRQ','GABRP','GABRR1','GABRR2','GABRR3')

# Glycine receptors
subgroup.16.genes <- c( 'GLRA1', 'GLRA2', 'GLRA3', 'GLRA4', 'GLRB')


# Ionotropic glutamate receptors

subgroup.17.genes  <- c( 'GRIA1', 'GRIA2','GRIA3', 'GRIA4', 'GRID1', 'GRID2','GRIK1',
                         'GRIK2', 'GRIK3', 'GRIK4', 'GRIK5', 'GRIN1', 'GRIN2A',
                         'GRIN2B', 'GRIN2C', 'GRIN2D', 'GRIN3A', 'GRIN3B' )


# IP3 receptor
subgroup.18.genes <- c('ITPR1','ITPR2','ITPR3')


# Nicotinic acetylcholine receptors
# 	Q03481 this is dropped
#   nicotinic acetylcholine receptor α8 subunit (avian)*
#   dropped
#---
subgroup.19.genes  <- c( 'CHRNA1','CHRNA2','CHRNA3','CHRNA4','CHRNA5','CHRNA6',
                         'CHRNA7', 'CHRNA9','CHRNA10', 'CHRNB1', 'CHRNB2', 'CHRNB3',
                         'CHRNB4', 'CHRNG', 'CHRND', 'CHRNE')


# P2X receptors
subgroup.20.genes <- c('P2RX1','P2RX2', 'P2RX3', 'P2RX4','P2RX5','P2RX6','P2RX7')

# ZAC
subgroup.21.genes <- c('ZACN')

# Aquaporins
subgroup.22.genes <- c('MIP','AQP1','AQP2','AQP3','AQP4','AQP5','AQP6','AQP7','AQP8','AQP9','AQP10')


# Chloride channels : ClC family
subgroup.23.genes <- c('CLCN1','CLCN2','CLCNKA','CLCNKB','CLCN3','CLCN4',
                        'CLCN5','CLCN6','CLCN7')

# CFTR

subgroup.24.genes <- c('CFTR')

# Calcium activated chloride channel

subgroup.25.genes <- c('ANO1')

# Maxi chloride channel
# droped : Maxi Cl-

# Volume regulated chloride channels
# dropped : VRAC

# Connexins and Pannexins

subgroup.26.genes <- c( 'GJE1','GJB7','GJB2', 'GJB6','GJC3','GJB4','GJB3','GJB5',
                        'GJD3', 'GJB1','GJD2','GJA4','GJA5','GJD4','GJA1','GJC1',
                        'GJA3','GJC2', 'GJA8','GJA9','GJA10','PANX1','PANX2','PANX3')
# Sodium leak channel, non-selective

subgroup.27.genes <- c('NALCN')





kcn.family.genes <- c(subgroup.1.genes, subgroup.2.genes, subgroup.3.genes, subgroup.4.genes)
kcn.family.genes <- subgroup.2.genes
kcn.family.ID    <- Ann$GeneID[Ann$SYMBOL %in% kcn.family.genes]
kcn.pathway.list <-  data.frame( diseaseId = unlist(rep('kcn.pathway', length(kcn.family.ID ))), 
                                             diseaseName   = unlist(kcn.family.ID))
gsea.result      <- GSEA(ion.gene.list ,TERM2GENE = kcn.pathway.list, maxGSSize = 1000, minGSSize = 8 )
# check if the kcn.family genes contained in the ion family
# kcn.family.genes %in% ion.genename
# reverse is to check the gene index in ion.genename
# ion.genename %in% kcn.family.genes


# annotate the public array data
#-----
setwd("D:\\wangli_data\\GSE79768")
raw.data    <- ReadAffy();
rma.data    <- affy::rma(raw.data);
exprs.data  <- exprs(rma.data)
probes      <- rownames(exprs.data)
gene.symbol <- unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound = NA))
rownames(exprs.data) <- gene.symbol

setwd("E:\\FuWai\\寿实验室\\ion_channel")
ion.list         <- read.csv( 'IonChannelsListDownloadForward.csv', 
                               stringsAsFactors = FALSE)
human.ion.genename     <- ion.list$Human.gene.name

AF.patient.exprs <- exprs.data[, seq(from = 2,  to = 14, by = 2)] %>% apply(1, median)
SR.patient.exprs <- exprs.data[, seq(from = 16, to = 26, by = 2)] %>% apply(1, median)
ion.group        <- rep(1, length(AF.patient.exprs))
ion.group[row.names(exprs.data) %in% human.ion.genename ] <- 2

gene.exprs.median         <- data.frame(Patient = AF.patient.exprs, Health = SR.patient.exprs, Type = ion.group )

vsmc <- ggplot(data = gene.exprs.median, aes(x = Health, y = Patient, color = as.factor(Type), alpha = as.factor(Type)) ) +
       geom_point( ) +
       scale_alpha_manual(values = c(1/200, 1), guide = FALSE) + 
       scale_colour_manual(name = 'gene groups',values = c("black","red"), 
                           labels = c('non-related genes','inon gene family')) + 

       geom_abline(intercept = 0.58, slope = 1, size = 1, alpha = 1/5) +
       geom_abline(intercept = -0.58, slope = 1, size = 1, alpha = 1/5) +
       theme(legend.position = c(0.8, 0.2),legend.title.align = 0.5)

genes.atrium     <- AF.patient.exprs[rownames(AF.patient.exprs ) %in% atrium.markers ,] %>%
                    apply(1,median)
genes.ventricle  <- AF.patient.exprs[rownames(AF.patient.exprs ) %in% ventricle.markers ,] %>%
                    apply(1,median)

#chamber.df       <- data.frame(atrium = )

boxplot(genes.atrium, genes.ventricle)
t.test(genes.atrium, genes.ventricle)




# public data cross-validation
# GSE2240
# PMID : 15877233   15817885
# title:
# ---

setwd("D:\\wangli_data\\GSE2240\\chipA")
raw.data      <- ReadAffy();
rma.data      <- affy::rma(raw.data);
exprs.data.A  <- exprs(rma.data)
probes        <- rownames(exprs.data.A)
gene.symbol   <- unlist(mget(probes, hgu133aSYMBOL, ifnotfound = NA))
rownames(exprs.data.A) <- gene.symbol

setwd("D:\\wangli_data\\GSE2240\\chipB")

raw.data      <- ReadAffy();
rma.data      <- affy::rma(raw.data);
exprs.data.B  <- exprs(rma.data)
probes        <- rownames(exprs.data.B)
gene.symbol   <- unlist(mget(probes, hgu133bSYMBOL, ifnotfound = NA))
rownames(exprs.data.B) <- gene.symbol


# 31:35 non-failing left ventricle
# 1:10 AF
# 11:30 SR
# 31:35 normal left ventricle

GSE2240.atrium.A    <- exprs.data.A[rownames(exprs.data.A) %in% atrium.markers,c(11:30)]    %>% apply(1,mean)
GSE2240.ventrile.A  <- exprs.data.A[rownames(exprs.data.A) %in% ventricle.markers,c(11:30)] %>% apply(1,mean)

boxplot(GSE2240.atrium.A, GSE2240.ventrile.A, names = c('atrium', 'ventricle'))
t.test(GSE2240.atrium.A, GSE2240.ventrile.A)

GSE2240.atrium.B    <- exprs.data.B[rownames(exprs.data.B) %in% atrium.markers, c(11:30)]    %>% apply(1,mean)
GSE2240.ventrile.B  <- exprs.data.B[rownames(exprs.data.B) %in% ventricle.markers,c(11:30)]  %>% apply(1,mean)

t.test(GSE2240.atrium.B ,GSE2240.ventrile.B)




# download the ENCODE project
# human total RNA samples RNA-seq (heart)
# download the files.txt, total 6 samples
# edit and erase the non-related data source
# and save as download.txt
# please checke the cheatsheet
# curl --help
# nohup xargs -n 1 curl -O -C - -L < download.txt &
# these raw data used as the reference data set
# to check if my QC procedure is correct!!!
# raw data from Encode is deposited here
#  /home/zhenyisong/biodata/cardiodata/encodeHeart/humanHeart
#  any reference please check yaoyanpostQC.mRNA.R
#  raw scripts: yaoyan.postQC.mRNA.R
# GSE87943 Human heart Tissue obtained from a healthy, 34 year old Caucasian male donor
# GSE78567 	Homo sapiens heart tissue female fetal (19 weeks)
# GSE88367 	Homo sapiens right atrium auricular region tissue female adult (53 years)
# GSE88247 	Homo sapiens heart left ventricle tissue female adult (53 years)
# GSE88271  Homo sapiens heart left ventricle tissue female adult (51 year)
#-----

setwd('D:\\wangli_data\\Rdata')
load('yaoyan.qc.Rdata')


human.sex.qc.matrix <- yaoyan.qc.gene.exprs.log[Ann$SYMBOL %in% y.chromosome.gene, ]

colnames(human.sex.qc.matrix)
colnames(human.sex.qc.matrix) <- c( 'male_adult_34Y','female_fetal_19W','female_fetal_19W.2',
                                    'female_adult_53Y','female_adult_53Y.x','female_adult_51Y')
sds.sex <- rowSds(human.sex.qc.matrix)
summary(sds.sex)
human.sex.qc.matrix[sds.sex > 0.5,]
km.res <- kmeans(t(human.sex.qc.matrix[sds.sex > .7,]), 2, nstart = 100)


fviz_nbclust( t(human.sex.qc.matrix[sds.sex > .7,]),  kmeans, method = "silhouette") +
              theme_classic()

fviz_cluster(km.res, data = t(human.sex.qc.matrix[sds.sex > .7,]),
             palette = c("#2E9FDF", "#FC4E07"),
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal()
)


# human AV comparision
# 	Homo sapiens right atrium auricular region tissue female adult (53 years)
# Sample GSM2343865
# 	
# GSE88367 4
human.atrium.qc.matrix <- yaoyan.qc.gene.exprs.log[Ann$SYMBOL %in% atrium.markers, 4]
human.ventricle.qc.matrix <- yaoyan.qc.gene.exprs.log[Ann$SYMBOL %in% ventricle.markers, 4]

boxplot(human.atrium.qc.matrix, human.ventricle.qc.matrix, names = c('atrium','ventricle'))

ks.test(human.atrium.qc.matrix, human.ventricle.qc.matrix)
#---


#
# aging validation
# human aging dataset was downloaded from
# http://www.senescence.info/
# aging information
#
#---
setwd('E:\\FuWai\\wangli.lab\\yaoyan')
human.aging.pub   <- read.csv('genage_human.csv', header = TRUE, stringsAsFactors = FALSE)
human.aging.names <- human.aging.pub$symbol

gene         <- yaoyan.mRNA.genes
gene.counts  <- gene$counts[,c(1:9, 13:15)]
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

aging.rpkm.mRNA           <- rpkm(gene.exprs, log = TRUE)
rownames(aging.rpkm.mRNA) <- make.names(Ann$SYMBOL, unique = TRUE)

aging.group  <- factor(rep(1:4, each = 3), levels = 1:4, labels = c('AF40','AF50','AF60','AF70'));
design       <- model.matrix(~ 0 + aging.group );
colnames(design) <- levels(aging.group)
contrast.matrix  <- makeContrasts( AF60 - AF40, AF70 - AF40, levels = design)
d.norm           <- voom(gene.exprs, design = design)
fit              <- lmFit(d.norm, design)
fit2             <- contrasts.fit(fit,contrast.matrix)
fit2             <- eBayes(fit2)

aging.result.1     <- topTable( fit2, coef    = 1,
                                number        = Inf, 
                                adjust.method = "BH", 
                                sort.by       = "p");

aging.result.1.names     <- filter(aging.result.1, P.Value < 0.05, abs(logFC) > 0.58) %>%
                            dplyr::select(SYMBOL) %>% na.omit()

aging.result.2     <- topTable( fit2, coef    = 2,
                                number        = Inf, 
                                adjust.method = "BH", 
                                sort.by       = "p");
aging.result.2.names     <- filter(aging.result.2, P.Value < 0.05, abs(logFC) > 0.58) %>%
                            dplyr::select(SYMBOL) %>% na.omit()

#GSEA
aging.whole.set.1         <- aging.result.1$logFC
names(aging.whole.set.1)  <- aging.result.1$GeneID
aging.whole.set.1         <- sort(aging.whole.set.1, decreasing = TRUE)

human.aging.geneID        <- aging.result.1$GeneID[aging.result.1$SYMBOL %in% human.aging.names]

aging.pathway.list <-  data.frame( Term = unlist(rep('aging.pathway', length(human.aging.geneID))), 
                                             Name   = unlist(human.aging.geneID))
aging.1.result      <- GSEA(aging.whole.set.1 ,TERM2GENE = aging.pathway.list, maxGSSize = 1000, minGSSize = 8,pvalueCutoff = 1 )

human.gsea.figure1B    <- gseaplot(aging.1.result, "aging.pathway")


#
# AF40 ~ AF70
#---
aging.whole.set.2         <- aging.result.2$logFC
names(aging.whole.set.2)  <- aging.result.2$GeneID
aging.whole.set.2         <- sort(aging.whole.set.2, decreasing = TRUE)

human.aging.geneID        <- aging.result.1$GeneID[aging.result.2$SYMBOL %in% human.aging.names]

aging.pathway.list <-  data.frame( Term = unlist(rep('aging.pathway', length(human.aging.geneID))), 
                                             Name   = unlist(human.aging.geneID))
aging.2.result      <- GSEA(aging.whole.set.2 ,TERM2GENE = aging.pathway.list, maxGSSize = 1000, minGSSize = 8, pvalueCutoff = 1 )

gseaplot(aging.2.result, "aging.pathway")

intersect(union(aging.result.1.names[[1]],aging.result.2.names[[1]]),human.aging.names)
intersect(aging.result.1.names[[1]],aging.result.2.names[[1]])

# read the sample inforamtion from
# Excel file deposited in 
#---
setwd('E:\\FuWai\\wangli.lab\\yaoyan\\AmericanBootcamp')
human.patient.meta <- read.xlsx( 'tables.xlsx', sheetName = 'T1',
                                  startRow = 3, endRow = 21 )
setwd('D:\\wangli_data\\Rdata')
save.image(file = "yaoyan.Rdata")
quit("no")




