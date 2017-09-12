# @author  Yisong Zhen
# @since   2017-09-12
# @update 
# I re-live this script according to the request by zongna;
# 
# the raw data were generated from next-seq 500
# SE or PE?

# nohup R CMD BATCH /home/zhenyisong/biodata/wanglab/wangcode/zongna.Ino80.R &


pkgs         <- c( 'tidyverse','Rsubread', 'edgeR', 'limma','DiagrammeR',
                   'gplots', 'genefilter', 'DESeq2', 'RColorBrewer',
                   'org.Rn.eg.db', 'cluster', 'factoextra','clusterProfiler',
                   'grid', 'gridExtra',
                   'pathview','sva')
load.lib     <- lapply(pkgs, require, character.only = TRUE)

rat.genome.path    <- file.path('/mnt/date/igenomes/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa')
rat.gtf            <- file.path('/mnt/date/igenomes/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf')
Rdata.output.dir   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/Rdata')
rsubread.index.lib <- file.path('/mnt/date/igenomes/rsubread')

"
setwd(Rdata.output.dir)
load('zonnagIno80.Rdata')
"

"
setwd('D:\\wangli_data\\Rdata')
load('zonnagIno80.Rdata')
"

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
rat.base.string      <- 'rn6'


setwd(rsubread.index.lib)
align( index         = rat.base.string, 
       readfile1     = rat.reads.paths, 
       input_format  = 'FASTQ', 
       type          = 'rna',
       output_file   = rat.outputs.files, 
       output_format = 'BAM', 
       nthreads      = 6, 
       indels        = 1,
       maxMismatches = 3,
       phredOffset   = 33,
       unique        = T )

zongna.Ino80.rat    <- featureCounts( rat.outputs.files, useMetaFeatures = TRUE, 
                                      annot.ext  = rat.gtf, isGTFAnnotationFile = TRUE,
                                      nthreads   = 6, allowMultiOverlap = TRUE)


setwd(Rdata.output.dir)
save.image(file = 'zonnagIno80.Rdata')
quit('no')

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

sample.info              <- data.frame( treat  = c('Control','Control',
                                                   'Treat','Treat') )
zongna.vsd               <- zongna.Ino80.rat %$% counts %>%
                            DESeqDataSetFromMatrix( countData = .,
                                                    colData   = sample.info,
                                                    design    = ~ treat) %>%
                            DESeq() %>% getVarianceStabilizedData()
                          
rownames(zongna.vsd)       <- NULL
colnames(zongna.vsd)       <- c('Control-1','Control-2',
                                'Ino80-1','Ino80-2')  
sds <- rowSds(zongna.vsd)
sh  <- shorth(sds)
vsd.filtered.exprs <- zongna.vsd[sds > 0.06,]
my.palette         <- rev(colorRampPalette(brewer.pal(10, 'RdBu'))(10))

zongna.table <- grid.newpage() %>% {cor(vsd.filtered.exprs)} %>%
                tableGrob(rows = colnames(vsd.filtered.exprs)) %>%
                grid.draw()

pheatmap( cor(vsd.filtered.exprs), cluster_rows = T, cluster_cols = T,
          color = my.palette, scale = 'none', fontsize_col = 8, 
          fontsize_row = 8, cellwidth = 30, cellheight = 30,
          labels_row = colnames(vsd.filtered.exprs), 
          labels_col = colnames(vsd.filtered.exprs),
          display_numbers = FALSE, number_color = 'orange',
          fontsize_number = 10)


"
heatmap( cor(vsd.filtered.exprs, method = 'spearman'), 
         cexCol = 0.8, cexRow = 0.8, col = my.palette)
"
pr <- prcomp(t(zongna.vsd))
plot(pr$x[,1:2], col = 'white', main = 'Sample PCA plot', type = 'p')
text(pr$x[,1], pr$x[,2], labels = colnames(zongna.vsd), cex = 0.7)
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

