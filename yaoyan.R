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
library(ggplot2)
library(gridExtra)
library(affy)
library(annotate)
library(org.Hs.eg.db)
library(VennDiagram)
library(Rsamtools)
library(clusterProfiler)
library(pathview)


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
GeneInfo     <- select( org.Hs.eg.db, keys= as.character(gene.ids), 
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
setwd("E:\\FuWai\\สูสตั้สา\\ion_channel")
ion.list         <- read.csv('IonChannelsListDownloadForward.csv')
ion.geneset.id   <- ion.list$HGNC.ID
ion.pathway.list <-  data.frame( diseaseId = unlist(rep('Ion.pathway', length(ion.geneset.id ))), 
                                             geneId   = unlist(ion.geneset.id) )
gsea.result      <- GSEA(AF.pathway.genelist ,TERM2GENE = ion.pathway.list, maxGSSize = 1000 )

# 

# negative result
#---end


save.image(file = "yaoyan.Rdata")
quit("no")




