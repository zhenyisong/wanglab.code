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
                   'grid', 'gridExtra', 'openxlsx', 'rat2302.db', 'annotate',
                   'affy',
                   'pathview','sva')
load.lib     <- lapply(pkgs, require, character.only = TRUE)

rat.genome.path    <- file.path('/mnt/date/igenomes/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa')
rat.gtf            <- file.path('/mnt/date/igenomes/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf')
Rdata.output.dir   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/Rdata')
rsubread.index.lib <- file.path('/mnt/date/igenomes/rsubread')

#---
# how to get the rat lncRNA annotation file from noncode database
#http://www.noncode.org/datadownload/NONCODEv5_rat_rn6_lncRNA.gtf.gz
#md5sum NONCODEv5_rat_rn6_lncRNA.gtf.gz
#30789eb01ab7f659788659401cabf02d  NONCODEv5_rat_rn6_lncRNA.gtf.gz
#md5sum rn6_UCSC_NONCODEv5_total.GTF
#60c419ee723ab60bf91c47869255edee  rn6_UCSC_NONCODEv5_total.GTF
#---
rat.lncRNAmRNA.gtf  <- file.path('/mnt/date/genomelib/annotation/rn6_UCSC_NONCODEv5_total.GTF')


setwd(Rdata.output.dir)
load('zonnagIno80.Rdata')


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
"
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


# GSEA method for wangli suggestions;
# data source 

#
# GSE82339   SRP076227 PE, Stranded;
#---

output.path             <- '/home/zhenyisong/biodata/wanglab/zongna/Ino80Data/rsubread'
setwd('/home/zhenyisong/biodata/wanglab/zongna/pubData/SRP076227')
GSE82339.files.names      <- list.files(pattern = '*.fastq$')
GSE82339.reads.paths      <- paste0(getwd(),'/',GSE82339.files.names)
GSE82339.1.files          <- GSE82339.reads.paths[grep('_1.fastq',GSE82339.reads.paths)]
GSE82339.2.files          <- GSE82339.reads.paths[grep('_2.fastq',GSE82339.reads.paths)]

GSE82339.output.files     <- basename(GSE82339.1.files) %>% 
                             sub(pattern = '_1.fastq', replacement = '') %>%
                             paste0(output.path,'/', . ,'.bam')


base.string      <- 'rn6'

setwd(rsubread.index.lib)
align( index          = base.string, 
       readfile1      = GSE82339.1.files , 
       readfile2      = GSE82339.2.files, 
       input_format   = 'FASTQ', 
       type           = 'rna',
       output_file    = GSE82339.output.files, 
       output_format  = 'BAM',
       PE_orientation = 'fr', 
       nthreads       = 20, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

#---
# check the strandness
# I used the following protocol to generate
# the rn6_UCSC.bed
# this file is used to judge the stradness
#
# letter source
# https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome/x8BxjfISECY
# One of our engineers suggests that you could try lifting the Ensembl 
# Genes track (table ensGene) from rn5 to rn6 to get some approximate 
# gene names. One way to go about this would be to follow these steps:

# 1. Visit the UCSC Table Browser at http://genome.ucsc.edu/cgi-bin/hgTables.
# 2. Select the following options:
#Clade: Mammal
#Genome: Rat
#Assembly: Mar. 2012 (RGSC 5.0/rn5)
#Group: Genes and Gene Predictions
#Track: Ensembl Genes
#Table: ensGene
#Region: genome
#Output Format: BED - browser extensible data
#Output File: ensGene_from_rn5.bed

#3. Click "get output".
#4. On the next page, leave the default setting to create 
# "one BED record per whole gene" and click "get BED".
#5. Your computer should now download a file titled 
# "ensGene_from_rn5.bed". The next step is to load this 
# into the LiftOver tool for conversion to the rn6 assembly.
#6. Browse to the LiftOver tool at 
# http://genome.ucsc.edu/cgi-bin/hgLiftOver (also available in the top menu under "Tools").
#7. Select the rat rn5 assembly as the "Original Assembly" and rat rn6 as the "New Assembly".
#8. Click the "choose file" button and select the "ensGene_from_rn5.bed" 
# file that you previously downloaded. Click "submit file".
#9. In the Results section, there should be a link to "View Conversions". 
# Click that link to download a BED file that contains the coordinates 
# of the ensGene track as lifted to the rn6 genome.
#10. Open the Custom Track tool at http://genome.ucsc.edu/cgi-bin/hgCustom.
#11. If you already have other custom tracks loaded, you will need 
# to click on the "add custom tracks" button. Otherwise, proceed to step 12.
#12. Select the rat rn6 genome assembly at the top of the "Add Custom Tracks" page, 
# then click the top "Choose File" button (next to "Paste URLs or data")
#  and select the file of lifted coordinates that you downloaded in step 9.
#13. Click "submit".

# You should now have a custom track loaded for the rn6 genome assembly 
# that contains an automatic conversion of the Ensembl genes track from rn5. 
# Please note that as this is an automated conversion, it will not have gone 
# through the same examination and testing that our standard track releases do.
#  It may still be enough, however, to provide you with the gene names that you are looking for.
#infer_experiment.py -r /mnt/date/genomelib/annotation/rn6_UCSC.bed \
#-i /home/zhenyisong/biodata/wanglab/zongna/Ino80Data/rsubread/SRR3643464.bam
#---
GSE82339.genes  <- featureCounts( GSE82339.output.files, 
                                      useMetaFeatures        = TRUE,
                                      countMultiMappingReads = FALSE,
                                      strandSpecific         = 0, 
                                      isPairedEnd            = TRUE,
                                      requireBothEndsMapped  = TRUE,
                                      autosort               = TRUE,
                                      nthreads               = 20,
                                      annot.ext              = rat.gtf,
                                      isGTFAnnotationFile    = TRUE, 
                                      GTF.featureType        = 'exon',
                                      allowMultiOverlap      = TRUE)


# GSE98575   SRP106502 SE,
# rat
#--

setwd('/home/zhenyisong/biodata/wanglab/zongna/pubData/SRP106502')
GSE98575.files.names    <- list.files(pattern = '*.fastq$') %>%
                           paste0(getwd(),'/',.)
GSE98575.outputs.files  <- basename(GSE98575.files.names) %>% 
                           sub(pattern = '.fastq', replacement = '') %>%
                           paste0(output.path,'/', . ,'.bam')
setwd(rsubread.index.lib)
align( index         = base.string, 
       readfile1     = GSE98575.files.names, 
       input_format  = 'FASTQ', 
       type          = 'rna',
       output_file   = GSE98575.outputs.files, 
       output_format = 'BAM', 
       nthreads      = 10, 
       indels        = 1,
       maxMismatches = 3,
       phredOffset   = 33,
       unique        = T )

GSE98575.genes    <- featureCounts( GSE98575.outputs.files, useMetaFeatures = TRUE, 
                                    annot.ext  = rat.gtf, isGTFAnnotationFile = TRUE,
                                    nthreads   = 6, allowMultiOverlap = TRUE)

# 
# GSE48111
#---
setwd('/home/zhenyisong/biodata/wanglab/zongna/pubData/GSE48111')
GSE48111.exprs           <- ReadAffy() %>% affy::rma() %>%
                            exprs()
GSE48111.gene.symbol     <- rownames(GSE48111.exprs) %>% 
                            {unlist(mget(., rat2302SYMBOL, ifnotfound = NA))} %>%
                            make.names(unique = TRUE)

setwd(Rdata.output.dir)
save.image(file = 'zonnagIno80.Rdata')
quit('no')

# completed processing all raw data
"




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
colnames(design)     <- c('Control','Ino80','Cellline')
contrast.matrix      <- makeContrasts(Ino80 - Control, levels = design)
zongnaIno80.result   <- zongna.Ino80.rat %$% counts %>% 
                        DGEList(genes = zongna.Ino80.rat$annotation) %>% 
                        calcNormFactors() %>% 
                        voom(design = design) %>%
                        lmFit(design) %>% 
                        contrasts.fit(contrast.matrix) %>%
                        eBayes() %>%
                        topTable( coef          = 1,
                                  number        = Inf, 
                                  adjust.method = 'BH', 
                                  sort.by       = 'none')

rat.GeneID <- mapIds( org.Rn.eg.db, keys = rownames(zongna.Ino80.rat$counts), 
                      column = 'ENTREZID', keytype = 'SYMBOL',multiVals = 'first' )
# GSEA analysis of 
raw.excel.path <- file.path('E:\\FuWai\\wangli.lab\\liangping\\cardiachypertrophyGSEA.xlsx')
#raw.excel.path <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/lianglab/cardiachypertrophyGSEA.xlsx')
cardio.ann.df  <- read.xlsx( raw.excel.path, sheet = 'mapping', 
                             startRow = 1, colNames = TRUE) %>%
                  filter(!(is.na(HumanSymbol))) %>%
                  filter(!(is.na(RatEnrezID)))

zongnaIno80.logFC           <- zongnaIno80.result$logFC
names(zongnaIno80.logFC)    <- rat.GeneID
zongnaIno80.logFC           <- sort(zongnaIno80.logFC, decreasing = TRUE)
cartrophy.geneset           <- data.frame( diseaseId = rep('cardiotropy', nrow(cardio.ann.df)), 
                                           geneId    = cardio.ann.df$RatEnrezID, check.names = TRUE)
zongnaIno80.gsea.result     <- zongnaIno80.logFC  %>% 
                               GSEA( TERM2GENE = cartrophy.geneset, 
                                     maxGSSize = 5000, pvalueCutoff = 1)
zongnaIno80.gsea.plot       <- gseaplot(zongnaIno80.gsea.result, 'cardiotropy')

setwd('C:\\Users\\Yisong\\Desktop')
write.xlsx(zongnaIno80.result , file = 'zongnaIno80.20170913.xlsx')

#---
# GSEA method used by wangli
# just for fun
# the three set to screen the ineresting results
#
#---

rat.GeneID <- mapIds( org.Rn.eg.db, keys = rownames(zongna.Ino80.rat$counts), 
                      column = 'ENTREZID', keytype = 'SYMBOL',multiVals = 'first' )

groups     <- factor(c(1,1,1,2,2,2),levels = 1:2, labels = c('V0','V48'))
design     <- model.matrix(~ 0 + groups)
colnames(design)     <- levels(groups)
contrast.matrix      <- makeContrasts(V48 - V0, levels = design)

GSE98575.result   <- GSE98575.genes %$% counts[,c(1:3,13:15)] %>% 
                     DGEList(genes = GSE98575.genes$annotation) %>% 
                     calcNormFactors() %>% 
                     voom(design = design) %>%
                     lmFit(design) %>% 
                     contrasts.fit(contrast.matrix) %>%
                     eBayes() %>%
                     topTable( coef          = 1,
                               number        = Inf, 
                               adjust.method = 'BH', 
                               sort.by       = 'none') %>%
                      mutate(EntrezID = rat.GeneID) %>%
                      arrange(desc(logFC))
GSE98575.logFC         <- GSE98575.result %>% dplyr::select(logFC) %>%
                          unlist()
names(GSE98575.logFC)  <- GSE98575.result %>% dplyr::select(EntrezID) %>%
                          unlist()


GSE82339.result  <- GSE82339.genes %$% counts[,c(1,3)] %>% 
                    DGEList(genes = GSE82339.genes$annotation) %>% 
                    calcNormFactors() %>% rpkm(log = TRUE)
GSE82339.logFC         <- GSE82339.result[,2] - GSE82339.result[,1]
names(GSE82339.logFC)  <- rat.GeneID
GSE82339.logFC         <- sort(GSE82339.logFC, decreasing = TRUE)

groups     <- factor(c(1,1,1,2,2,2), levels = 1:2, labels = c('dmso_90','dmsoPE_2880'))
design     <- model.matrix(~ 0 + groups)
colnames(design)      <- levels(groups)
contrast.matrix       <- makeContrasts(dmsoPE_2880 - dmso_90, levels = design)
rat.GeneID.microarray <- mapIds( org.Rn.eg.db, keys = GSE48111.gene.symbol, 
                                 column = 'ENTREZID', keytype = 'SYMBOL',multiVals = 'first' )
GSE48111.result   <- GSE48111.exprs[,c(1:3,15:17)] %>% 
                     lmFit(design) %>% 
                     contrasts.fit(contrast.matrix) %>%
                     eBayes() %>%
                     topTable( coef          = 1,
                               number        = Inf, 
                               adjust.method = 'BH', 
                               sort.by       = 'none') %>%
                      mutate(EntrezID = rat.GeneID.microarray) %>%
                      arrange(desc(logFC))
GSE48111.logFC         <- GSE48111.result %>% dplyr::select(logFC) %>%
                          unlist()
names(GSE48111.logFC)  <- GSE48111.result %>% dplyr::select(EntrezID) %>%
                          unlist()



cell.line  <- factor(c(1,2,1,2),levels = 1:2, labels = c('cell1','cell2'))
groups     <- factor(c(1,1,2,2),levels = 1:2, labels = c('Control','Ino80'))
design     <- model.matrix(~ 0 + groups + cell.line);
colnames(design)     <- c('Control','Ino80','Cellline')
contrast.matrix      <- makeContrasts(Ino80 - Control, levels = design)
zongna.GSEA.set      <- zongna.Ino80.rat %$% counts %>% 
                        DGEList(genes = zongna.Ino80.rat$annotation) %>% 
                        calcNormFactors() %>% 
                        voom(design = design) %>%
                        lmFit(design) %>% 
                        contrasts.fit(contrast.matrix) %>%
                        eBayes() %>%
                        topTable( coef          = 1,
                                  number        = Inf, 
                                  adjust.method = 'BH', 
                                  sort.by       = 'none') %>%
                         mutate(EntrezID = rat.GeneID) %>%
                         filter(P.Value < 0.05) %>%
                         filter(abs(logFC) > 0.58) %>%
                         dplyr::select(EntrezID)

Ino80.set                <- data.frame( diseaseId = rep('Ino80', length(zongna.GSEA.set)), 
                                        geneId    = zongna.GSEA.set, check.names = TRUE)

GSE98575.gsea.result     <- GSE98575.logFC  %>% 
                            GSEA( TERM2GENE = Ino80.set, 
                                  maxGSSize = 5000, pvalueCutoff = 1)
GSE98575.gsea.plot       <- gseaplot(GSE98575.gsea.result, 'Ino80')  

GSE48111.gsea.result     <- GSE48111.logFC  %>% 
                            GSEA( TERM2GENE = Ino80.set, 
                                  maxGSSize = 5000, pvalueCutoff = 1)
GSE48111.gsea.plot       <- gseaplot(GSE48111.gsea.result, 'Ino80') 


GSE82339.gsea.result     <- GSE82339.logFC  %>% 
                            GSEA( TERM2GENE = Ino80.set, 
                                  maxGSSize = 5000, pvalueCutoff = 1)
GSE82339.gsea.plot       <- gseaplot(GSE82339.gsea.result, 'Ino80')           