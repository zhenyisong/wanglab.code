# the raw data were generated from next-seq 500
# SE or PE?

library(gplots)
library(xlsx)
library(Rsubread)
library(edgeR)
library(limma)
library(DESeq2)
library(gplots)
library(genefilter)
library(RColorBrewer)
library(org.Rn.eg.db)
library(cluster)
library(factoextra)
library(clusterProfiler)
library(pathview)
library(sva)


setwd('/home/zhenyisong/biodata/wanglab/zongna/Ino80Data/rsubread')
setwd("D:\\wangli_data\\Rdata")
load('zonnagIno80.Rdata')
save.image(file = 'zonnagIno80.Rdata')
quit('no')


rat.genome_ref.path   <- '/mnt/date/igenomes/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa'
rat.gtf               <- '/mnt/date/igenomes/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf'

#---test code
# for sva analysis
# to be implemented


#---
# this code was branched from xinli.multiple.R
# see Fuwi Supercomputer
# zongna
#---

setwd('/home/zhenyisong/biodata/wanglab/zongna/Ino80Data/rawdata')
reads.files.names    <- list.files(pattern = '*.fastq$')
rat.reads.paths      <- paste0(getwd(),'/',reads.files.names)
setwd('/home/zhenyisong/biodata/wanglab/zongna/Ino80Data')
unlink('rsubread', force = TRUE, recursive = TRUE)
dir.create('rsubread')
output.path          <- '/home/zhenyisong/biodata/wanglab/zongna/Ino80Data/rsubread'
setwd('/home/zhenyisong/biodata/wanglab/zongna/Ino80Data/rsubread')
rat.outputs.files    <- paste0(output.path,'/', reads.files.names,'.bam')
rat.base.string      <- 'rn6_index'

buildindex( basename = rat.base.string, reference = rat.genome_ref.path )
align( index         = rat.base.string, 
       readfile1     = rat.reads.paths, 
       input_format  = "FASTQ", 
       type          = 'rna',
       output_file   = rat.outputs.files, 
       output_format = "BAM", 
       nthreads      = 15, 
       indels        = 1,
       maxMismatches = 3,
       phredOffset   = 33,
       unique        = T )

gene.rat           <- featureCounts( rat.outputs.files, useMetaFeatures = TRUE, 
                                      annot.ext  = rat.gtf, isGTFAnnotationFile = TRUE,
                                      nthreads   = 15, allowMultiOverlap = TRUE)

gene         <- gene.rat
gene.counts  <- gene$counts
gene.ids     <- gene$annotation$GeneID

Ann          <- gene$annotation[, c("GeneID", "Chr","Length")]

Ann$Chr      <- unlist( lapply(strsplit(Ann$Chr, ";"), 
                    function(x) paste(unique(x), collapse = "|")))
Ann$Chr      <- gsub("chr", "", Ann$Chr)

gene.exprs   <- DGEList(counts = gene.counts, genes = Ann)
gene.exprs   <- calcNormFactors(gene.exprs)

rpkm.genes           <- rpkm(gene.exprs, log = TRUE)
colnames(rpkm.genes) <- c('Control1','Control2','Ino801','Ino802')
rpkm.lfc.genes       <- cbind( rpkm.genes[,3] - rpkm.genes[,1], rpkm.genes[,4]- rpkm.genes[,2],
                               (rpkm.genes[,3] + rpkm.genes[,4])/2 - (rpkm.genes[,1] + rpkm.genes[,2])/2 )
colnames(rpkm.lfc.genes) <- c('Ino801 - Control1','Ino802 -Control2','mean(Ino801+Ino802) - mean(C1+C2)')

setwd("C:\\Users\\Yisong\\Desktop")
write.xlsx(rpkm.lfc.genes, file = 'zongna.lfc.xlsx', sheetName = 'log_fold_change')
write.xlsx(rpkm.genes, file = 'zongna.lfc.xlsx',sheetName = 'log_transf_data', append = TRUE)

#----
# QC
# PCA
#---

dge.tmm    <- t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
dge.tmm.counts           <- apply(dge.tmm,2, as.integer)
sample.info              <- data.frame( treat  = c('Control','Control',
                                                   'Treat','Treat') )
dds                      <- DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                                    colData   = sample.info,
                                                    design    = ~ treat)
vsd                      <- varianceStabilizingTransformation(dds);
vsd.expr                 <- assay(vsd)
#rownames(vsd.expr)       <- gene.exprs$genes$SYMBOL
rownames(vsd.expr)       <- NULL
colnames(vsd.expr)       <- c('Control-1','Control-2',
                              'Ino80-1','Ino80-2')  
sds <- rowSds(vsd.expr)
sh  <- shorth(sds)
vsd.filtered.expr <- vsd.expr[sds > 0.3,]
my.palette        <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
heatmap( cor(vsd.filtered.expr, method = 'spearman'), 
         cexCol = 0.8, cexRow = 0.8, col = my.palette)
pr <- prcomp(t(vsd.expr))
plot(pr$x[,1:2], col = 'white', main = 'Sample PCA plot', type = "p")
text(pr$x[,1], pr$x[,2], labels = colnames(vsd.expr), cex = 0.7)
pr.var <- pr$sdev^2
pve    <- pr.var/sum(pr.var)
plot(pve, xlab = 'PCA', ylab = 'Variance explained', type = 'b')


#
# DGE analysis
#---
cell.line  <- factor(c(1,2,1,2),levels = 1:2, labels = c('cell1','cell2'))
groups     <- factor(c(1,1,2,2),levels = 1:2, labels = c('Control','Ino80'))
design     <- model.matrix(~ 0 + groups + cell.line);
colnames(design) <- c('Control',"Ino80",'Cellline')
contrast.matrix  <- makeContrasts(Ino80 - Control, levels = design)
d.norm           <- voom(gene.exprs, design = design)
fit              <- lmFit(d.norm, design)
fit2             <- contrasts.fit(fit,contrast.matrix)
fit2             <- eBayes(fit2)
gene.result      <- topTable(  fit2, 
                               number        = Inf, 
                               adjust.method = "BH", 
                               sort.by       = "p"
                               );

setwd("C:\\Users\\Yisong\\Desktop")
write.xlsx(gene.result, file = 'zongnaIno80.benchdeleted.xlsx')
# Her data is not consistancy, and use p.value as the cutoff will
# lead to topTable Error!!!

