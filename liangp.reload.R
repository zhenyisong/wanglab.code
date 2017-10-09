# @author Yisong Zhen
# @since  2017-06-22
# @update 2017-08-15
#
# @rawdata PE model, Human
# contact with Dr. Liang, Ping, his paired-end project
# the preprocessing was completed by yaofang
# nohup R CMD BATCH /home/zhenyisong/biodata/wanglab/wangcode/liangping.restart.R &

# data were demultiplexed by 
# using the parameter

# the rawdata is not completed
# find ./ -type f -name "*LP*"
# @since 2017-06-20
#---

# task assigned by Liang.P
# email message in hotmail account
#
# 一松，
# 
# 请按照以下分组进行数据分析，
# 
# 一、CIB1课题
# 1、比较WT H9 Basal（3、7）与H9-CIB1-KO Basal（1、5）
# 2、比较WT H9 ＋ PMA（4、8）与H9-CIB1-KO PMA（2、6）
# 
# 二、TRPC1课题
# 1、比较WT H9 Basal（3、7）与H9-TRPC1-KO Basal（9、11）
# 2、比较WT H9 ＋ PMA（4、8）与H9-TRPC1-KO PMA（10、12）
# 
# 三、心肌肥大模型验证
# 比较WT H9 Basal（3、7）与WT H9 ＋ PMA（4、8）
# 
# 另外，1和2、3和4、5和6、7和8、9和10、11和12均为配对的样本
# （准备的相同的心肌细胞，1组为Basal，另外1组为PMA处理）。
# 

# efficient R programming
# package management
pkgs <- c( 'tidyverse','Rsubread','org.Hs.eg.db','edgeR',
           'limma', 'DESeq2', 'genefilter','grid',
           'openxlsx','pheatmap','gridExtra','ggrepel',
           'QuasR','annotate','clusterProfiler',
           'cowplot','tools',
           'GGally','RColorBrewer','RDAVIDWebService',
           'cluster','factoextra',"ggpubr",
           'BSgenome.Hsapiens.UCSC.hg38',
           'BSgenome.Hsapiens.UCSC.hg38.Rbowtie')
#install.package(pkgs)
# this will install all necessary packages
# required in this project
#---
#install.lib <- lapply(pkgs, library, character.only = TRUE)
install.lib <- lapply(pkgs, require, character.only = TRUE)
# data transfer 
# in window 10
#---

x1.runing.path <- file.path('D:\\wangli_data\\Rdata')
setwd(x1.runing.path)
load('liangp.Rdata')

#deprecated
#liangp.output.path     <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/lianglab/nextdata')
#setwd(liangp.output.path)
#load('liangp.graphs.Rdata')
#---

#

#--

# read the raw data into the working enviroment
# and classsify the data into read group, 1, 2
#---

human.genome.path      <- file.path('/mnt/date/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa')
liangp.output.path     <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/lianglab/nextdata')
liangp.QC.path         <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/lianglab/nextdata/multiQC')

human.rawdata.path.1   <- file.path('/mnt/date/Sequencing/FastQ/20170605_LF_RNAseq_LPBD')
liangp.files.1         <- list.files( path = human.rawdata.path.1, pattern = "h\\_.*\\.gz$", all.files = FALSE,
                                      full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = F) %>%
                                      `[`(-c(13:24))

human.rawdata.path.2   <- file.path('/mnt/date/Sequencing/FastQ/170616_LF_RNAseq')
liangp.files.2         <- list.files( path = human.rawdata.path.2, pattern = "LP-2.*$", all.files = FALSE,
                                      full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = F)

human.rawdata.path.3   <- file.path('/mnt/date/Sequencing/FastQ/20170621_LF_RNAseq')
liangp.files.3         <- list.files( path = human.rawdata.path.3, pattern = "LP-1.*$", all.files = FALSE,
                                      full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = F)

liangp.files           <- c(liangp.files.1, liangp.files.2, liangp.files.3)
# add this line at 2017-07-21
#--
rawdata.figureprints   <- md5sum(liangp.files)
read.1.files           <- liangp.files[grep("R1",liangp.files)]
read.2.files           <- liangp.files[grep("R2",liangp.files)]

human.output.filenames <- basename(read.1.files) %>% sub(pattern = '_R1_001.fastq.gz', replacement = '') %>%
                          paste0(liangp.output.path,'/', . ,'.bam')
sample.names           <- basename(read.1.files) %>% sub(pattern = '_R1_001.fastq.gz', replacement = '')

# preQC by nextflow
# cd /home/zhenyisong/biodata/wanglab/wangdata/lianglab/nextdata/nextflowQC
# ln -s /mnt/date/Sequencing/FastQ/20170605_LF_RNAseq_LPBD/h_* 
# rm h_vsmc*
# ln -s /mnt/date/Sequencing/FastQ/170616_LF_RNAseq/LP-2* ./
# ln -s /mnt/date/Sequencing/FastQ/20170621_LF_RNAseq/LP-1* ./
# ln -s /mnt/s *.gz |wc -l
# 24
# nextflow.sh -p 'liangp' -g hg38 -r "*{R1,R2}*.fastq.gz"

# new full whole QC procedure
# plase check the QC results
# /home/zhenyisong/biodata/wanglab/wangdata/lianglab/nextdata/nextflowQC
#---

# preQC by QuasR
#
sampleFile      <- tempfile(pattern = "zhen3temp", tmpdir = tempdir())
sample.file     <- data.frame( FileName1  = unlist(read.1.files),
                               FileName2  = unlist(read.2.files), 
                               SampleName = unlist(sample.names))
write_tsv(sample.file, path = sampleFile)
genome          <- 'BSgenome.Hsapiens.UCSC.hg38'
cluster         <- makeCluster(18)

dir.create(file.path(liangp.QC.path), showWarnings = FALSE)
setwd(file.path(liangp.QC.path))

# this running code if checkOlny set to FALSE
# will lead to at the first time to
# prepare the in-house bowtie index files
# for exmaple, to create the ‘BSgenome.Hsapiens.UCSC.hg38.Rbowtie'
# library(BSgenome.Hsapiens.UCSC.hg38.Rbowtie)
# after this time ,checkOnly = FALSE
#---

liangp.qPorject <- qAlign( sampleFile,
                           genome,
                           auxiliaryFile = NULL,
                           aligner = "Rbowtie",
                           maxHits = 1,
                           paired  = 'fr',
                           splicedAlignment = FALSE,
                           snpFile = NULL,
                           bisulfite = "no",
                           alignmentParameter = NULL,
                           projectName = "qProject",
                           alignmentsDir = liangp.QC.path,
                           lib.loc  = NULL,
                           cacheDir = NULL,
                           clObj = cluster,
                           checkOnly = F)
                           
qQCReport( liangp.qPorject, pdfFilename = 'liangp.QC.pdf', 
           useSampleNames = TRUE, clObj = cluster)

save.image('QClinagp.Rdata')
quit('no')
#-- read data end


# build genome index
#
setwd(liangp.output.path)
base.string          <-  'hg38'
buildindex( basename = base.string, reference = human.genome.path)
align( index          = base.string, 
       readfile1      = read.1.files, 
       readfile2      = read.2.files, 
       input_format   = "gzFASTQ", 
       type           = 'rna',
       output_file    = human.output.filenames, 
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 15, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

# take  gene count
#
liangp.genes      <- featureCounts( human.output.filenames, useMetaFeatures = TRUE,
                                    countMultiMappingReads = FALSE,
                                    strandSpecific         = 0, 
                                    isPairedEnd            = TRUE,
                                    requireBothEndsMapped  = TRUE,
                                    autosort               = TRUE,
                                    nthreads               = 15,
                                    annot.inbuilt = "hg38", allowMultiOverlap = TRUE)




gene         <- liangp.genes
gene.counts  <- gene$counts
gene.ids     <- gene$annotation$GeneID


# this need to be updated to pipeline
# use the inner_join call?
# which would be much flexible
#---

columns      <- c("ENTREZID","SYMBOL", "GENENAME");
GeneInfo     <- AnnotationDbi::select( org.Hs.eg.db, keys= as.character(gene.ids), 
                        keytype = "ENTREZID", columns = columns);
m            <- match(gene$annotation$GeneID, GeneInfo$ENTREZID);
Ann          <- cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                       GeneInfo[m, c("SYMBOL", "GENENAME")]);

Ann$Chr      <- unlist( lapply(strsplit(Ann$Chr, ";"), 
                        function(x) paste(unique(x), collapse = "|")))
Ann$Chr      <- gsub("chr", "", Ann$Chr)

# post-QC check for grouping
# P3 project, cardiac hypertrophy
# QC check
#---
sample.p3.liangp  <- data.frame( treat = c('WT.H9.Basal','WT.H9.Basal','WT.H9.PMA','WT.H9.PMA') )

liangp.p3.vst     <- liangp.genes %$% counts[,7:10] %>%
                     DESeqDataSetFromMatrix(colData = sample.p3.liangp, design = ~ treat) %>%
                     DESeq() %>%
                     varianceStabilizingTransformation() %>% assay()
colnames(liangp.p3.vst) <- c('WT.H9.Basal-3','WT.H9.Basal-7','WT.H9.PMA-4','WT.H9.PMA-8')
sds <- rowSds(liangp.p3.vst)
sh  <- shorth(sds)
(sh)
# [1] 0.2169877
liangp.p3.pca          <- liangp.p3.vst %>% subset(sds > 0.3) %>%
                          t() %>% prcomp

pca.p3.groups          <- c('WT.H9.Basal-3','WT.H9.Basal-7','WT.H9.PMA-4','WT.H9.PMA-8')


# deprecated, 
# basic ploting function
"
dev.off()
heatmap( cor(liangp.p3.vst.filtered), 
         labRow = pca.p3.groups, labCol = pca.p3.groups)
plot(liangp.p3.pca$x[,c(1,2)], ylim = c(-80,45),xlim = c(-90,180))
text( liangp.p3.pca$x[,1], liangp.p3.pca$x[,2], 
      labels = pca.p3.groups , cex = 0.5, pos = 3, adj = c(0,1))
"

liangp.p3.pve    <- liangp.p3.pca$sdev^2/sum(liangp.p3.pca$sdev^2)

pve.p3.df        <- data.frame( variance = liangp.p3.pve, 
                                pca      = 1:length(liangp.p3.pve))


# ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
# please use the above key words to search hits
#---
pca.p3.labels        <- paste0('PCA-',c(1:length(liangp.p3.pve)))
names(pca.p3.labels) <- c(1:length(liangp.p3.pve))
p3.screeplot <- ggplot(pve.p3.df) +
        xlab('Principle Component') +
        ylab('Proportion of Variance Explained') +
        scale_x_discrete( breaks = c(1:length(liangp.p3.pve)), 
                          labels = pca.p3.labels, 
                          limits = as.character(c(1:length(liangp.p3.pve)))) +
        geom_point(aes(x = pca, y = variance), size = 3) +
        geom_line(aes(x = pca, y = variance), size = 0.8) +
        scale_linetype_discrete() +
        theme(legend.position = "none") +
        theme_grey()

p3.screeplot

color.bar <- colorRampPalette(c("midnightblue", "grey", "mediumvioletred"), alpha = TRUE)(10) 
draw_colnames_45 <- function (coln, gaps, ...) {
    coord <- pheatmap:::find_coordinates( length(coln), gaps)
    x     <- coord$coord - 0.5 * coord$size
    res   <- textGrob( coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), 
                       vjust   = 0.5, hjust = 1, rot = 45, gp = gpar(...))
    return(res)
}  

assignInNamespace(x = "draw_colnames", value = "draw_colnames_45",
                  ns = asNamespace("pheatmap"))          
p3.sample.heatmap <- liangp.p3.vst %>% subset(sds > 0.3) %>%
                     cor() %>%
                     pheatmap( cluster_rows = T, cluster_cols = T,
                               color = c(color.bar, alpha = 0.1), scale = 'none', fontsize_col = 9, 
                               fontsize_row = 9, cellwidth = 40, cellheight = 40,
                               labels_row = pca.p3.groups, labels_col = pca.p3.groups,
                               display_numbers = TRUE, number_color = 'orange',
                               fontsize_number = 10)
p3.sample.heatmap

color.pal          <- rainbow(5)
p3.sample.corrmap  <- liangp.p3.vst %>% subset(sds > 0.3) %>% 
                      ggcorr( method  =  c("pairwise", "spearman"), 
                              label = TRUE, limits = c(0.7, 1))

p3.sample.corrmap
# this is deprecated
# only one point in cetroid
# please check this reponse
# https://groups.google.com/forum/#!topic/r-help-archive/SFGlrypOGS0
#---

"
p3.cluster.estimation  <- liangp.p3.vst %>%
                          t() %>%
                          fviz_nbclust(pam, method = 'silhouette') +
                          theme_classic()
"


p3.pam.cluster <- liangp.p3.vst %>% subset(sds > 0.3) %>%
                  t() %>% pam(diss = FALSE, k = 3) %>% 
                  fviz_cluster( pallete = c('darkmagenta','gold','firebrick'),
                                repel   = TRUE,
                                xlab    = 'PCAd.1, 68%',
                                ylab    = 'PCAd.2, 21%',
                                legend  = NULL,
                                ggtheme = theme_classic())
ggpar(p3.pam.cluster, legend = 'none', legend.title = NULL)
                                                 


pca.maps <- function(pca.index, pca.res, pca.group) {
     pca.check.df <- data.frame( pcaX = pca.res$x[,1], 
                                 pcaY = pca.res$x[,pca.index])
     graph <- ggplot(data = pca.check.df, aes(x = pcaX, y = pcaY)) + 
              geom_point(col = 'blue', size = 2, show.legend = F) +
              geom_text_repel(label = as.character(pca.group)) +
              xlab('PCAd.1') +
              ylab( paste0('PCAd.',pca.index,sep = '') )
              theme_classic()
     invisible(return(graph))  
}

p3.pca.num <- c('W3','W7','P4','P8')
all.PCA.grobs <- map( c(2:4), pca.maps, 
                      pca.res   = liangp.p3.pca,
                      pca.group = p3.pca.num)
pca.inOne <- grid.arrange(grobs = all.PCA.grobs, nrow = 2, ncol = 2)



#pca.maps(pca.df = liangp.pca,pca.group = pca.num, 3)
# here map can be substitued with lapply
# see purrr
# and other R package and function usage
#---

# read the curated cardiac hypertrophy genes
#---

# data curation
# date source
# https://www.ncbi.nlm.nih.gov/pubmed/18629135?dopt=Citation
# PMID: 18629135 
# RawExcel:
# 

raw.excel.path <- file.path("E:\\FuWai\\wangli.lab\\liangping\\cardiachypertrophyGSEA.xlsx")
#raw.excel.path <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/lianglab/cardiachypertrophyGSEA.xlsx')
cardio.ann.df  <- read.xlsx( raw.excel.path, sheet = 'mapping', 
                             startRow = 1, colNames = TRUE) %>%
                  filter(!(is.na(HumanSymbol)))


# limma analysis for DGE in Project 3
# use pipeline to fast and clearify study
#---
p3.reload.group   <- factor( c( 'WT.H9.Basal','WT.H9.Basal','WT.H9.PMA','WT.H9.PMA'), 
                             levels = c( 'WT.H9.Basal','WT.H9.PMA'));
#cell.line;
p3.reload.design             <- model.matrix(~ 0 + p3.reload.group)
colnames(p3.reload.design)   <- levels(p3.reload.group)
reload.contrast.matrix       <- makeContrasts(WT.H9.Basal - WT.H9.PMA, levels = p3.reload.design)

liangp.p3.result             <- liangp.genes %$% counts[,7:10] %>% DGEList(genes = Ann) %>% 
                                calcNormFactors() %>% 
                                voom(design = p3.reload.design) %>%
                                lmFit(p3.reload.design) %>% 
                                contrasts.fit(reload.contrast.matrix) %>%
                                eBayes() %>%
                                topTable( coef          = 1,
                                          number        = Inf, 
                                          adjust.method = "BH", 
                                          sort.by       = "p")


# output of Project 3 in excel
setwd("C:\\Users\\Yisong\\Desktop")
write.xlsx(liangp.p3.result, file = "liangp.p3.xlsx", colNames = TRUE, borders = "columns")

liangp.p3.logFC             <- liangp.p3.result$logFC
names(liangp.p3.logFC)      <- liangp.p3.result$GeneID
liangp.p3.logFC             <- sort(liangp.p3.logFC, decreasing = TRUE)
cartrophy.geneset           <- data.frame( diseaseId = rep('cardiotropy', nrow(cardio.ann.df)), 
                                           geneId    = cardio.ann.df$HumanEntrezID, check.names = TRUE)
p3.gsea.result              <- liangp.p3.logFC %>% 
                               GSEA( TERM2GENE = cartrophy.geneset, 
                                     maxGSSize = 5000, pvalueCutoff = 1)
p3.gsea.plot                <- gseaplot(p3.gsea.result, 'cardiotropy')

# single out the sample 4 index = 3
#                       7, index = 4

p3.rpkm.log <- liangp.genes %$% counts[,7:10] %>% DGEList(genes = Ann) %>% 
                                calcNormFactors() %>% rpkm(log = TRUE)
liangp.sample4.gsea        <- (p3.rpkm.log[,4] - apply(p3.rpkm.log[,1:2], 1, mean)) %>%
                              sort(decreasing = T) %>% 
                              GSEA( TERM2GENE = cartrophy.geneset, 
                                    maxGSSize = 5000, pvalueCutoff = 1)


(liangp.sample4.gsea)
# the above indicate that the
# cardiac hypertrophy pathyway was stimulated
# then I tried to confirm using the whole gene epxression level;
#---

p3.rpkm.log <- liangp.genes %$% counts[,7:10] %>% DGEList(genes = Ann) %>% 
                                calcNormFactors() %>% rpkm(log = TRUE) %>%
                                as.data.frame() %>%
                                filter(rownames(.) %in% cardio.ann.df$HumanEntrezID)
cardiotropy.tidy     <- data.frame( gene.exprs  = c(apply(p3.rpkm.log[,1:2], 1, mean),
                                                    apply(p3.rpkm.log[,3:4], 1, mean)),
                                    status.d    = factor( rep(c(1,2), 
                                                              each = nrow(p3.rpkm.log)),
                                                          levels = 1:2,
                                                          labels = c('normal','disease')))

p3.cardio.boxplot   <- ggplot(data = cardiotropy.tidy, aes(x = status.d, y = gene.exprs)) + 
                       geom_violin(trim = FALSE) +
                       geom_jitter(position = position_jitter(.2) ) +
                       stat_summary(fun.data = 'mean_sdl', fun.args = list(mult = 1),
                                     geom = 'pointrange', color = 'red')
                                                           
#norm.data    <-   apply(p3.rpkm.log[,1:2], 1, mean)
#disease.data <-   apply(p3.rpkm.log[,2:4], 1, mean)
#wilcox.test( norm.data, disease.data, paired = T)                                                  
#boxplot(gene.exprs ~ status.d , cardiotropy.tidy)
nonparam.result      <- wilcox.test( gene.exprs ~ status.d , 
                                    data = cardiotropy.tidy, paired = T)
(nonparam.result)

# manual check the cardiac hypertrophy genes
# curated base upone the recent review
# PMID: 24663494
# PLoS One. 2014 Mar 24;9(3):e92903. 
# A systematic review of fetal genes as biomarkers of 
# cardiac hypertrophy in rodent models of diabetes.
# another
# title:
# A personalized, multi-omics approach identifies genes 
# involved in cardiac hypertrophy and heart failure
# http://www.biorxiv.org/content/early/2017/03/24/120329.full.pdf+html
# 
# PMID: 23258295 
# title:
# Molecular basis of physiological heart growth: fundamental concepts and new players.
# 
# ---

hg38.cardiotrophy.geneID        <- c('4878','4879','3280','4625','58')
names(hg38.cardiotrophy.geneID) <- c('NPPA', 'NPPB','HES1','MYH7','ACTA1')

liangp.p3.counts  <- liangp.genes %$% counts[,7:10] %>%
                     DESeqDataSetFromMatrix(colData = sample.p3.liangp, design = ~ treat) %>%
                     DESeq() %>% counts(normalized = TRUE) %>% apply(c(1,2), as.integer)

(liangp.p3.counts[hg38.cardiotrophy.geneID,])

cardiotrophy.markers           <- liangp.p3.counts[hg38.cardiotrophy.geneID,] 
dimnames(cardiotrophy.markers) <- list( names(hg38.cardiotrophy.geneID),
                                        c( 'WT.H9.Basal-3','WT.H9.Basal-7',
                                           'WT.H9.PMA-4','WT.H9.PMA-8') )

cardiotrophy.marker.table      <- tableGrob( cardiotrophy.markers)
grid.newpage()
grid.draw(cardiotrophy.marker.table)

# deleting 4/8 respectively
#---

# set deleing sample id
# 3 is the sample-4
# 4 is corresponding to sample 8
#---
sample.id   <- 3 # this is 8
p3.rpkm.log <- liangp.genes %$% counts[,7:10] %>% DGEList(genes = Ann) %>% 
                                calcNormFactors() %>% rpkm(log = TRUE) %>%
                                as.data.frame() %>%
                                filter(rownames(.) %in% cardio.ann.df$HumanEntrezID)
cardiotropy.tidy.del <- data.frame( gene.exprs  = c(apply(p3.rpkm.log[,1:2], 1, mean),
                                                    p3.rpkm.log[,sample.id]),
                                    status.d    = factor( rep(c(1,2), 
                                                              each = nrow(p3.rpkm.log)),
                                                          levels = 1:2,
                                                          labels = c('normal','disease')))
nonparam.result.del  <- wilcox.test( gene.exprs ~ status.d , 
                                     data = cardiotropy.tidy.del, paired = T)
(nonparam.result.del)


# 

liangp.p3.log     <- liangp.genes %$% counts[,7:10] %>%
                     DESeqDataSetFromMatrix(colData = sample.p3.liangp, design = ~ treat) %>%
                     DESeq() %>% counts(normalized = TRUE) %>% apply(c(1,2), as.integer) %>%
                     add(1) %>% log



cardiotrophy.markers.log <- liangp.p3.log[hg38.cardiotrophy.geneID,] 
cardiotrophy.markers.df  <- data.frame( avg   = apply(cardiotrophy.markers.log[,-3], 1, median), 
                                        ratio = cardiotrophy.markers.log[,4] / 
                                                apply(cardiotrophy.markers.log[,-3], 1, median))
PMA.model.df     <- data.frame( avg   = apply(liangp.p3.log[,-3], 1, median), 
                                ratio = liangp.p3.log[,4] / apply(liangp.p3.log[,-3], 1, median)) %>%
                    filter(complete.cases(.) & !is.infinite(rowSums(.))) %>%
                    filter((ratio <= 4 & avg >= .4))

# remove illagal letter, which can call the function in package
# library(IDPmisc)
# NaRV.omit(df)
# is much better?
#---


p3.sample.dotplot <- ggplot(data = PMA.model.df, aes( x= avg, y = ratio)) +
                     geom_point(alpha = 0.1) +
                     geom_point( data  = cardiotrophy.markers.df, aes(x =  avg , y = ratio), 
                                 color = 'red', size = 2, shape = 16) +
                     geom_text_repel( data  = cardiotrophy.markers.df, aes(x =  avg , y = ratio),
                                      label = names(hg38.cardiotrophy.geneID), 
                                      color = 'green3', size = 3.5, 
                                      fontface = "bold") +
                     xlab('average expression counts (log)') + 
                     ylab('ratio PMA/Basal (log)') +
                     theme_classic()
(p3.sample.dotplot)
# output the result in Excel file
# using new package
#---
#setwd("C:\\Users\\Yisong\\Desktop")
#contrast.list <- c('WT.H9.Basal - CIB1.KO.Basal',  'WT.H9.PMA - CIB1.KO.PMA',
#                   'WT.H9.Basal - TRPC1.KO.Basal', 'WT.H9.PMA - TRPC1.KO.PMA',
#                   'WT.H9.PMA   - TRPC1.KO.PMA',    'WT.H9.Basal - WT.H9.PMA')
#
#liangp.wb <- createWorkbook("Yisong")
#
#for( i in 1:6) {
#   addWorksheet(liangp.wb, contrast.list[i])
#   liangp.result <-  topTable( liangp.results, coef = i, adjust = "BH", 
#                               number = Inf, sort.by = "P")
#   writeData(liangp.wb, i, liangp.result)
#}
#saveWorkbook(liangp.wb, file = "liangp.xlsx", overwrite = TRUE)
#
# end of Project II
# ---


# Project II
# TRPC1
# 二、TRPC1课题
# 1、比较WT H9 Basal（3、7）与H9-TRPC1-KO Basal（9、11）
# 2、比较WT H9 ＋ PMA（4、8）与H9-TRPC1-KO PMA（10、12）
#---

sample.p2.1.liangp  <- data.frame( treat = c('WT.H9.Basal','WT.H9.Basal','H9.TRPC1.KO','H9.TRPC1.KO') )

liangp.p2.1.vst     <- liangp.genes %$% counts[,c(7,8,4,3)] %>%
                       DESeqDataSetFromMatrix(colData = sample.p2.1.liangp, design = ~ treat) %>%
                       DESeq() %>%
                       varianceStabilizingTransformation() %>% assay()
colnames(liangp.p2.1.vst) <- c('WT.H9.Basal-3','WT.H9.Basal-7','H9.TRPC1.KO-9','H9.TRPC1.KO-11')
sds <- rowSds(liangp.p2.1.vst)
sh  <- shorth(sds)
(sh)


liangp.p2.1.pca          <- liangp.p2.1.vst %>% subset(sds > 0.1) %>%
                            t() %>% prcomp

pca.p2.1.groups          <- c('WT.H9.Basal-3','WT.H9.Basal-7','H9.TRPC1.KO-9','H9.TRPC1.KO-11')

plot( liangp.p2.1.pca$x[,c(1,2)], 
      xlim = c(-90,100), ylim = c(-60, 75), col = 'red', pch = 16 )
text( liangp.p2.1.pca$x[,1], liangp.p2.1.pca$x[,2], 
      labels = pca.p2.1.groups , cex = 0.75, pos = 3, adj = c(0,1))


#--- END project p2.1

sample.p2.2.liangp  <- data.frame( treat = c('WT.H9.PMA','WT.H9.PMA','H9.TRPC1.PMA','H9.TRPC1.PMA') )

liangp.p2.2.vst     <- liangp.genes %$% counts[,c(9,10,5,6)] %>%
                       DESeqDataSetFromMatrix(colData = sample.p2.2.liangp, design = ~ treat) %>%
                       DESeq() %>%
                       varianceStabilizingTransformation() %>% assay()
colnames(liangp.p2.2.vst)     <- c('WT.H9.PMA-4','WT.H9.PMA-8','H9.TRPC1.PMA-10','H9.TRPC1.PMA-12')

sds <- rowSds(liangp.p2.2.vst)
sh  <- shorth(sds)
(sh)

liangp.p2.2.pca          <- liangp.p2.2.vst %>% subset(sds > 0.3) %>%
                            t() %>% prcomp
liangp.p2.2.pve         <- liangp.p2.2.pca$sdev^2/sum(liangp.p2.2.pca$sdev^2)
(liangp.p2.2.pve)

pca.p2.2.groups          <- c('WT.H9.PMA-4','WT.H9.PMA-8','H9.TRPC1.PMA-10','H9.TRPC1.PMA-12')

plot( liangp.p2.2.pca$x[,c(1,2)], 
      xlim = c(-80,200), ylim = c(-100, 80), col = 'green', pch = 16)
text( liangp.p2.2.pca$x[,1], liangp.p2.2.pca$x[,2], 
      labels = pca.p2.2.groups , cex = 0.75, pos = 3, adj = c(0,1))


# start project 1
# 一、CIB1课题
# 1、比较WT H9 Basal（3、7）与H9-CIB1-KO Basal（1、5）
# 2、比较WT H9 ＋ PMA（4、8）与H9-CIB1-KO PMA（2、6）
#---

sample.p1.1.liangp  <- data.frame( treat = c('WT.H9.Basal','WT.H9.Basal','H9.CIB1.KO','H9.CIB1.KO') )

liangp.p1.1.vst     <- liangp.genes %$% counts[,c(7,8,12,1)] %>%
                       DESeqDataSetFromMatrix(colData = sample.p1.1.liangp, design = ~ treat) %>%
                       DESeq() %>%
                       varianceStabilizingTransformation() %>% assay()
colnames(liangp.p1.1.vst) <- c('WT.H9.Basal-3','WT.H9.Basal-7','H9.CIB1.KO-1','H9.CIB1.KO-5')
sds <- rowSds(liangp.p1.1.vst)
sh  <- shorth(sds)
(sh)


liangp.p1.1.pca          <- liangp.p1.1.vst %>% subset(sds > 0.15) %>%
                            t() %>% prcomp

liangp.p1.1.pve         <- liangp.p1.1.pca$sdev^2/sum(liangp.p1.1.pca$sdev^2)
(liangp.p1.1.pve)

pca.p1.1.groups          <- c('WT.H9.Basal-3','WT.H9.Basal-7','H9.CIB1.KO-1','H9.CIB1.KO-5')

plot( liangp.p1.1.pca$x[,c(1,2)], 
      xlim = c(-70,135), ylim = c(-70, 50), col = 'blue', pch = 16 )
text( liangp.p1.1.pca$x[,1], liangp.p1.1.pca$x[,2], 
      labels = pca.p1.1.groups , cex = 0.75, pos = 3, adj = c(0,1))

#--- END p-1.1

sample.p1.2.liangp  <- data.frame( treat = c('WT.H9.PMA','WT.H9.PMA','H9.CIB1.PMA','H9.CIB1.PMA') )

liangp.p1.2.vst     <- liangp.genes %$% counts[,c(9,10,11,2)] %>%
                       DESeqDataSetFromMatrix(colData = sample.p1.2.liangp, design = ~ treat) %>%
                       DESeq() %>%
                       varianceStabilizingTransformation() %>% assay()
colnames(liangp.p1.2.vst) <- c('WT.H9.PMA-4','WT.H9.PMA-8','H9.CIB1.PMA-2','H9.CIB1.PMA-6')
sds <- rowSds(liangp.p1.2.vst)
sh  <- shorth(sds)
(sh)


liangp.p1.2.pca          <- liangp.p1.2.vst %>% subset(sds > 0.2) %>%
                            t() %>% prcomp
liangp.p1.2.pve          <- liangp.p1.2.pca$sdev^2/sum(liangp.p1.2.pca$sdev^2)
(liangp.p1.2.pve)

pca.p1.2.groups          <- c('WT.H9.PMA-4','WT.H9.PMA-8','H9.CIB1.PMA-2','H9.CIB1.PMA-8')

plot( liangp.p1.2.pca$x[,c(1,2)], 
      xlim = c(-80,190), ylim = c(-60, 110), col = 'orange', pch = 16 )
text( liangp.p1.2.pca$x[,1], liangp.p1.2.pca$x[,2], 
      labels = pca.p1.2.groups , cex = 0.75, pos = 3, adj = c(0,1))


# now perform the DGE analysis
# see liangp email:
# 课题1：TRPC1基因敲除对于心肌细胞功能及PMA诱导的心肌肥厚中的作用研究
# 1、WT H9 Basal：#3，#7
# 2、WT H9 ＋ PMA：#4，#8
# 3、H9-TRPC1-KO Basal：#9，#11
# 4、H9-TRPC1-KO PMA：#10，#12
#---

# in responding to liangp email Aug, 07
# I changed the contrast matrix
# 2017-08-07
#--
wp1.reload.group               <- factor( c( 'WT.H9.Basal','WT.H9.Basal','WT.H9.PMA','WT.H9.PMA',
                                            'H9.TRPC1.KO','H9.TRPC1.KO','H9.TRPC1.PMA','H9.TRPC1.PMA'), 
                                          levels = c( 'WT.H9.Basal','WT.H9.PMA','H9.TRPC1.KO','H9.TRPC1.PMA'));
wp1.reload.design              <- model.matrix(~ 0 + wp1.reload.group)

colnames(wp1.reload.design)    <- levels(wp1.reload.group)
wp1.contrast.matrix            <- makeContrasts( H9.TRPC1.KO - H9.TRPC1.PMA, WT.H9.Basal - WT.H9.PMA,
                                                 levels = wp1.reload.design)
liangp.wp1.result              <- liangp.genes %$% counts[,c(7:10,3:6)] %>% 
                                  DGEList(genes = Ann) %>% 
                                  calcNormFactors() %>% 
                                  voom(design = wp1.reload.design) %>%
                                  lmFit(wp1.reload.design) %>% 
                                  contrasts.fit(wp1.contrast.matrix ) %>%
                                  eBayes()
liangp.wp1.rpkm                <- liangp.genes %$% counts[,c(7:10,3:6)] %>% 
                                  DGEList(genes = Ann) %>% 
                                  calcNormFactors() %>% rpkm()
colnames(liangp.wp1.rpkm)      <- c('WT.H9.Basal-3','WT.H9.Basal-7','WT.H9.PMA-4','WT.H9.PMA-8',
                                    'H9.TRPC1.KO-11','H9.TRPC1.KO-9','H9.TRPC1.PMA-11','H9.TRPC1.PMA-12')
liangp.wp1.result.1            <- topTable( liangp.wp1.result , 
                                            coef          = 1,
                                            number        = Inf, 
                                            adjust.method = 'BH', 
                                            sort.by       = 'none') %>%
                                  cbind(liangp.wp1.rpkm) %>%
                                  arrange(P.Value)


liangp.wp1.result.2            <- topTable( liangp.wp1.result , 
                                            coef          = 2,
                                            number        = Inf, 
                                            adjust.method = 'BH', 
                                            sort.by       = 'none') %>%
                                  cbind(liangp.wp1.rpkm) %>%
                                  arrange(P.Value)

#---
# 课题2：CIB1基因敲除对于心肌细胞功能及PMA诱导的心肌肥厚中的作用研究
# 1、WT H9 Basal：#3，#7
# 2、WT H9 ＋ PMA：#4，#8
# 3、H9-CIB1-KO Basal：#1，#5
# 4、H9-CIB1-KO PMA：#2，#6  
#---

wp2.reload.group               <- factor( c( 'WT.H9.Basal','WT.H9.Basal','WT.H9.PMA','WT.H9.PMA',
                                            'H9.CIB1.KO','H9.CIB1.KO','H9.CIB1.PMA','H9.CIB1.PMA'), 
                                          levels = c('WT.H9.Basal','WT.H9.PMA','H9.CIB1.KO','H9.CIB1.PMA'));
wp2.reload.design              <- model.matrix(~ 0 + wp2.reload.group)

colnames(wp2.reload.design)    <- levels(wp2.reload.group)
wp2.contrast.matrix            <- makeContrasts( H9.CIB1.KO - H9.CIB1.PMA, WT.H9.Basal - WT.H9.PMA,
                                                 levels = wp2.reload.design)
liangp.wp2.result              <- liangp.genes %$% counts[,c(7:10,5,12,2,11)] %>% 
                                  DGEList(genes = Ann) %>% 
                                  calcNormFactors() %>% 
                                  voom(design = wp2.reload.design) %>%
                                  lmFit(wp2.reload.design) %>% 
                                  contrasts.fit(wp2.contrast.matrix ) %>%
                                  eBayes()
liangp.wp2.rpkm                <- liangp.genes %$% counts[,c(7:10,5,12,2,11)] %>% 
                                  DGEList(genes = Ann) %>% 
                                  calcNormFactors() %>% rpkm()
colnames(liangp.wp2.rpkm )     <- c( 'WT.H9.Basal-3','WT.H9.Basal-7','WT.H9.PMA-4','WT.H9.PMA-8',
                                     'H9.CIB1.KO-11','H9.CIB1.KO-1','H9.CIB1.PMA-6','H9.CIB1.PMA-2')
liangp.wp2.result.1            <- topTable( liangp.wp2.result , 
                                            coef          = 1,
                                            number        = Inf, 
                                            adjust.method = 'BH', 
                                            sort.by       = 'none') %>%
                                  cbind(liangp.wp2.rpkm) %>%
                                  arrange(P.Value)

liangp.wp2.result.2            <- topTable( liangp.wp2.result , 
                                            coef          = 2,
                                            number        = Inf, 
                                            adjust.method = 'BH', 
                                            sort.by       = 'none') %>%
                                 cbind(liangp.wp2.rpkm) %>%
                                 arrange(P.Value)

# export the result to two sheets
#---

setwd('C:\\Users\\Yisong\\Desktop')
liangp.p1.wb <- createWorkbook()
liangp.p2.wb <- createWorkbook()
addWorksheet(liangp.p1.wb, 'H9.TRPC1.KO - WT.H9.Basal')
addWorksheet(liangp.p1.wb, 'H9.TRPC1.PMA - WT.H9.PMA')
addWorksheet(liangp.p2.wb, 'H9.CIB1.KO - WT.H9.Basal')
addWorksheet(liangp.p2.wb, 'WT.H9.PMA - H9.CIB1.PMA')

writeData(liangp.p1.wb, sheet = 1, liangp.wp1.result.1 )
writeData(liangp.p1.wb, sheet = 2, liangp.wp1.result.2 )
writeData(liangp.p2.wb, sheet = 1, liangp.wp2.result.1 )
writeData(liangp.p2.wb, sheet = 2, liangp.wp2.result.2 )
saveWorkbook(liangp.p1.wb, 'liangp.p1.2017-08-07.xlsx', overwrite = TRUE)
saveWorkbook(liangp.p2.wb, 'liangp.p2.2017-08-07.xlsx', overwrite = TRUE)

# 1、Basal水平的比较：WT与KO比较
# 2、PMA水平的比较：WT与KO比较
# 之后，我们选择在2中的阳性、
# 但在1中不是阳性的基因/信号通路作为候选基因/信号通路
# see email from liangp
# Basal level p > 0.05
# PMA level   p < 0.05
# !!!! 
# contrast.matrix
# WP1, H9.TRPC1.KO - WT.H9.Basal, H9.TRPC1.PMA - WT.H9.PMA,
# WP2, H9.CIB1.KO - WT.H9.Basal, H9.CIB1.PMA - WT.H9.PMA,
#
#---

p1.df <- inner_join( liangp.wp1.result.1 %>% 
                     filter(P.Value > 0.05) %>%
                     filter(SYMBOL != ''),
                     liangp.wp1.result.2 %>% 
                     filter(P.Value < 0.05) %>%
                     filter(SYMBOL != ''),
                     by = 'SYMBOL')

p2.df <- inner_join( liangp.wp2.result.1 %>% 
                     filter(P.Value > 0.05) %>%
                     filter(SYMBOL != ''),
                     liangp.wp2.result.2 %>% 
                     filter(P.Value < 0.05) %>%
                     filter(SYMBOL != ''),
                     by = 'SYMBOL')

setwd('C:\\Users\\Yisong\\Desktop')
liangp.wb <- createWorkbook()
addWorksheet(liangp.wb, 'TRPC1-Basal')
addWorksheet(liangp.wb, 'CIB1-Basal')

writeData(liangp.wb, sheet = 1, p1.df )
writeData(liangp.wb, sheet = 2, p2.df )
saveWorkbook(liangp.wb, 'liangp.2017-08-09.xlsx', overwrite = TRUE)



# in response to liangp and wangli at
# 2017-08-14
# I sent an email to confirm this 
# contrast.matrix
# WP1,
# 1) H9.TRPC1.KO - H9.TRPC1.PMA, 2) WT.H9.Basal - WT.H9.PMA
# 2) positive, 1) negative
# WP2,
# H9.CIB1.KO - H9.CIB1.PMA, WT.H9.Basal - WT.H9.PMA,
#---



p1.1.df <- liangp.wp1.result.1 %>% 
           filter(P.Value < 0.05) %>%
           filter(abs(logFC) > 0.58) %>%
           filter(SYMBOL != '')

p1.2.df <- liangp.wp1.result.1 %>% 
           filter(P.Value < 0.05) %>%
           filter(abs(logFC) > 0.58) %>%
           filter(SYMBOL != '')


p2.1.df <- liangp.wp2.result.1 %>% 
           filter(P.Value < 0.05) %>%
           filter(abs(logFC) > 0.58) %>%
           filter(SYMBOL != '')

p2.2.df <- liangp.wp2.result.1 %>% 
           filter(P.Value < 0.05) %>%
           filter(abs(logFC) > 0.58) %>%
           filter(SYMBOL != '')

p1.geneID   <- inner_join( liangp.wp1.result.1 %>% 
                         filter(P.Value > 0.05) %>%
                         filter(abs(logFC) < 0.58) %>%
                         filter(SYMBOL != ''),
                         liangp.wp1.result.2 %>% 
                         filter(P.Value < 0.05) %>%
                         filter(abs(logFC) > 0.58) %>%
                         filter(SYMBOL != ''),
                         by = 'SYMBOL') %>% dplyr::select(GeneID.x) %>%
                         unlist

p1.geneID.df <- inner_join( liangp.wp1.result.1 %>% 
                         filter(P.Value > 0.05) %>%
                         filter(abs(logFC) < 0.58) %>%
                         filter(SYMBOL != ''),
                         liangp.wp1.result.2 %>% 
                         filter(P.Value < 0.05) %>%
                         filter(abs(logFC) > 0.58) %>%
                         filter(SYMBOL != ''),
                         by = 'SYMBOL')

p2.geneID <- inner_join( liangp.wp2.result.1 %>% 
                         filter(P.Value > 0.05) %>%
                         filter(abs(logFC) < 0.58) %>%
                         filter(SYMBOL != ''),
                         liangp.wp2.result.2 %>% 
                         filter(P.Value < 0.05) %>%
                         filter(abs(logFC) > 0.58) %>%
                         filter(SYMBOL != ''),
                         by = 'SYMBOL') %>% dplyr::select(GeneID.x) %>%
                         unlist

p2.geneID.df <- inner_join( liangp.wp2.result.1 %>% 
                         filter(P.Value > 0.05) %>%
                         filter(abs(logFC) < 0.58) %>%
                         filter(SYMBOL != ''),
                         liangp.wp2.result.2 %>% 
                         filter(P.Value < 0.05) %>%
                         filter(abs(logFC) > 0.58) %>%
                         filter(SYMBOL != ''),
                         by = 'SYMBOL') 

liangp.p1.kegg.tidy      <- enrichKEGG( p1.geneID, 
                                        organism = 'human', 
                                        pvalueCutoff  = 0.05, 
                                        pAdjustMethod = 'none',
                                        qvalueCutoff  = 1) %>%
                           summary() %$%
                           {data.frame( kegg.pvalue  = -log(pvalue),
                                        kegg.pathway = Description )}
liangp.p1.kegg.ggplot   <- liangp.p1.kegg.tidy %>% 
                           ggplot( aes( x = reorder(kegg.pathway, kegg.pvalue), 
                                                       y = kegg.pvalue)) + 
                           geom_bar( stat = 'identity', width = 0.4, 
                                     position = position_dodge(width = 0.1), size = 20) +
                           theme(axis.text.x = element_text(angle = 60,hjust = 1, size = 8),
                                 axis.text.y = element_text(hjust = 1, size = 10)) +
                           ylab('-log(pvalue)') + 
                           xlab('liangp.p1.KEGG.ggplot') + 
                           coord_flip()
liangp.p1.GO.tidy       <- enrichGO( gene  = p1.geneID,
                                     OrgDb = org.Hs.eg.db,
                                     ont   = "CC",
                                     pAdjustMethod = 'none',
                                     qvalueCutoff  = 1) %>%
                            summary() %$%
                            {data.frame( Golabels  = paste(ID, Description, sep = ' '),
                                         LogPvalue = - log(pvalue) ) } %>%
                            ggplot( aes( x = reorder(Golabels, LogPvalue), 
                                                       y = LogPvalue)) + 
                            geom_bar( stat = 'identity', width = 0.4, 
                                     position = position_dodge(width = 0.1), size = 20) +
                            theme(axis.text.x = element_text(angle = 60,hjust = 1, size = 8),
                                 axis.text.y = element_text(hjust = 1, size = 10)) +
                            ylab('-log(pvalue)') + 
                            xlab('liangp.p1.GO.ggplot2') + 
                            coord_flip()

# p2 ggplot
liangp.p2.kegg.tidy      <- enrichKEGG( p2.geneID, 
                                        organism = 'human', 
                                        pvalueCutoff  = 0.05, 
                                        pAdjustMethod = 'none',
                                        qvalueCutoff  = 1) %>%
                           summary() %$%
                           {data.frame( kegg.pvalue  = -log(pvalue),
                                        kegg.pathway = Description )}

liangp.p2.kegg.ggplot   <- liangp.p2.kegg.tidy %>% 
                           ggplot( aes( x = reorder(kegg.pathway, kegg.pvalue), 
                                                       y = kegg.pvalue)) + 
                           geom_bar( stat = 'identity', width = 0.4, 
                                     position = position_dodge(width = 0.1), size = 20) +
                           theme(axis.text.x = element_text(angle = 60,hjust = 1, size = 8),
                                 axis.text.y = element_text(hjust = 1, size = 10)) +
                           ylab('-log(pvalue)') + 
                           xlab('liangp.p2.KEGG.ggplot') + 
                           coord_flip()
liangp.p2.GO.tidy       <- enrichGO( gene  = p2.geneID,
                                     OrgDb = org.Hs.eg.db,
                                     ont   = "CC",
                                     pAdjustMethod = 'none',
                                     qvalueCutoff  = 1) %>%
                            summary() %$%
                            {data.frame( Golabels  = paste(ID, Description, sep = ' '),
                                         LogPvalue = - log(pvalue) ) } %>%
                            ggplot( aes( x = reorder(Golabels, LogPvalue), 
                                                       y = LogPvalue)) + 
                            geom_bar( stat = 'identity', width = 0.4, 
                                     position = position_dodge(width = 0.1), size = 20) +
                            theme(axis.text.x = element_text(angle = 60,hjust = 1, size = 8),
                                 axis.text.y = element_text(hjust = 1, size = 10)) +
                            ylab('-log(pvalue)') + 
                            xlab('liangp.p2.GO.ggplot2') + 
                            coord_flip()


liangp.p1.kegg.df        <- enrichKEGG( p1.geneID, 
                                        organism = 'human', 
                                        pvalueCutoff  = 0.05, 
                                        pAdjustMethod = 'none',
                                        qvalueCutoff  = 1) %>%
                            summary()
liangp.p1.GO.df          <- enrichGO( gene  = p1.geneID,
                                     OrgDb = org.Hs.eg.db,
                                     ont   = "CC",
                                     pAdjustMethod = 'none',
                                     qvalueCutoff  = 1) %>%
                            summary()
liangp.p2.kegg.df        <- enrichKEGG( p2.geneID, 
                                        organism = 'human', 
                                        pvalueCutoff  = 0.05, 
                                        pAdjustMethod = 'none',
                                        qvalueCutoff  = 1) %>%
                            summary()
liangp.p2.GO.df          <- enrichGO( gene  = p2.geneID,
                                     OrgDb = org.Hs.eg.db,
                                     ont   = "CC",
                                     pAdjustMethod = 'none',
                                     qvalueCutoff  = 1) %>%
                            summary()

setwd('C:\\Users\\Yisong\\Desktop')
liangp.wb <- createWorkbook()
addWorksheet(liangp.wb, 'H9.TRPC1.KO - H9.TRPC1.PMA')
addWorksheet(liangp.wb, 'WT.H9.Basal - WT.H9.PMA-I')

addWorksheet(liangp.wb, 'H9.CIB1.KO - H9.CIB1.PMA')
addWorksheet(liangp.wb, 'WT.H9.Basal - WT.H9.PMA-II')




addWorksheet(liangp.wb, 'p1.kegg')
addWorksheet(liangp.wb, 'p1.GO')
addWorksheet(liangp.wb, 'p2.kegg')
addWorksheet(liangp.wb, 'p2.GO')

addWorksheet(liangp.wb, 'p1.positive.set')
addWorksheet(liangp.wb, 'p2.positive.set')

writeData(liangp.wb, sheet = 1, p1.1.df )
writeData(liangp.wb, sheet = 2, p1.2.df )
writeData(liangp.wb, sheet = 3, p2.1.df )
writeData(liangp.wb, sheet = 4, p2.2.df )
writeData(liangp.wb, sheet = 5, liangp.p1.kegg.df )
writeData(liangp.wb, sheet = 6, liangp.p1.GO.df )
writeData(liangp.wb, sheet = 7, liangp.p2.kegg.df )
writeData(liangp.wb, sheet = 8, liangp.p2.GO.df )
writeData(liangp.wb, sheet = 9, p1.geneID.df )
writeData(liangp.wb, sheet = 10, p2.geneID.df )
saveWorkbook(liangp.wb, 'liangp.2017-08-16.xlsx', overwrite = TRUE)


# full samples normliaztion and vsn
#---
whole.samples.exprs <- liangp.genes %$% counts %>% 
                       DGEList(genes = Ann) %>% 
                       calcNormFactors() %>% rpkm()
samples.mean.rlog  <- cbind( WT.H9.Basal  = apply(whole.samples.exprs[,7:8], 1, mean),
                             WT.H9.PMA    = apply(whole.samples.exprs[,9:10], 1, mean),
                             H9.TRPC1.PMA = apply(whole.samples.exprs[,5:6], 1, mean),
                             H9.CIB1.PMA  = apply(whole.samples.exprs[,c(2,11)], 1, mean),
                             H9.CIB1.KO   = apply(whole.samples.exprs[,c(1,12)], 1, mean),
                             H9.TRPC1.KO  = apply(whole.samples.exprs[,3:4], 1, mean) ) %>%
                      apply(c(1,2), as.integer) %>% rlog()
(cardio.ann.df)

cardiotropy.data <- match(cardio.ann.df$HumanSymbol, Ann$SYMBOL)  %>% 
                    {samples.mean.rlog[.,]} 
rownames(cardiotropy.data) <- cardio.ann.df$HumanSymbol
cardiotropy.data           <- na.omit(cardiotropy.data)
sd.filter.cardiotrophy     <- apply(cardiotropy.data, 1, sd)

sd.filter.whole            <- apply(samples.mean.rlog, 1, sd)

liangp.final.4.groups      <- colnames(cardiotropy.data[,c(1,2,3,6)])
liangp.final.4.pca         <- cardiotropy.data[,c(1,2,3,6)] %>%
                              t() %>% prcomp 

plot( liangp.final.4.pca$x[,c(1,2)], 
      xlim = c(-5,7), ylim = c(-3,3),
      pch = 16, col = 'blue')
text( liangp.final.4.pca$x[,1], liangp.final.4.pca$x[,2], 
      labels = liangp.final.4.groups, 
      cex    = 0.5, pos = 3, adj = c(0,1))




color.bar <- colorRampPalette(c("blue4", "white", "springgreen4"))(10) 
trpc1.pheatmap <- cardiotropy.data[,-c(4:5)]  %>% cor(., method = 'spearman') %>%
                  pheatmap( cluster_rows = T, cluster_cols = T,
                            color = c(color.bar), scale = 'none', fontsize_col = 8, 
                            fontsize_row = 8, cellwidth = 30, cellheight = 30,
                            labels_row = colnames(cardiotropy.data)[-c(4:5)], 
                            labels_col = colnames(cardiotropy.data)[-c(4:5)],
                            display_numbers = TRUE, number_color = 'orange',
                            fontsize_number = 10)
cib1.pheatmap <- cardiotropy.data[,c(1,2,4,5)]  %>% cor(., method = 'spearman') %>%
                  pheatmap( cluster_rows = T, cluster_cols = T,
                            color = c(color.bar), scale = 'none', fontsize_col = 8, 
                            fontsize_row = 8, cellwidth = 30, cellheight = 30,
                            labels_row = colnames(cardiotropy.data)[c(1,2,4,5)], 
                            labels_col = colnames(cardiotropy.data)[c(1,2,4,5)],
                            display_numbers = TRUE, number_color = 'orange',
                            fontsize_number = 10)

wp.kegg.func    <- . %>% 
                  filter(P.Value < 0.05 & abs(logFC) > 0.58) %>%
                  dplyr::select(GeneID) %>% unlist %>%
                  enrichKEGG( organism = "hsa")
dotplot.func    <- .%>% dotplot( x = "count", showCategory = 20, colorBy = "qvalue")
wp.kegg.plots   <- map( list( liangp.wp1.result.1, liangp.wp1.result.2, 
                              liangp.wp2.result.1, liangp.wp2.result.2), wp.kegg.func) %>%
                   map(dotplot.func) %>% {plot_grid( plotlist = .[1:3], labels = c('A','B','C'),
                                                    ncol = 2, nrow = 2)}

# GO analysis include
# BP
# CC
# MF
#---
wp.GO.func      <- . %>% 
                   filter(P.Value < 0.05 & abs(logFC) > 0.58) %>%
                   dplyr::select(GeneID) %>% unlist %>%
                   enrichGO( OrgDb = 'org.Hs.eg.db', ont = "CC")
dotplot.func   <- . %>% dotplot( showCategory = 10)
wp.GO.plots    <- map( list( liangp.wp1.result.1, liangp.wp1.result.2, 
                             liangp.wp2.result.1, liangp.wp2.result.2), wp.GO.func) %>%
                  map(dotplot.func) %>% arrangeGrob(grobs = ., nrow = 2, ncol = 2) %>%
                  grid.arrange()
                 
ggdraw() + draw_plot(wp.GO.plots[[1]], x = 0, y = 0.5, width = .5, height = .5) +
           draw_plot(wp.GO.plots[[2]], x = .5, y = .5, width = .5, height = .5) +
           draw_plot(wp.GO.plots[[3]], x = 0, y = 0, width = .5, height = .5)  +
           draw_plot(wp.GO.plots[[4]], x = .5, y = 0, width = .5, height = .5) 

#---
# heatmap reload
# 

select.names.func <- . %>% filter(P.Value < 0.05 & logFC > 0.58) %>% dplyr::select(SYMBOL)
all.gene.names    <- list( liangp.wp1.result.1, liangp.wp1.result.2, 
                           liangp.wp2.result.1, liangp.wp2.result.1) %>%
                     map(select.names.func) %>% unlist %>% na.omit
dge.samples.rlog  <- match(all.gene.names, Ann$SYMBOL) %>% samples.mean.rlog[.,]
dge.samples.sd    <- apply(dge.samples.rlog, 1, sd)
dge.samples.rlog  <- dge.samples.rlog[dge.samples.sd > 0.3,]

heatmap.2( dge.samples.rlog , col = greenred(75),scale  = 'row', 
						     Rowv = TRUE,Colv = FALSE, density.info = 'none',key = TRUE, trace = 'none', 
						     cexCol = 1.5,distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						     hclustfun = function(d) hclust(d, method = 'complete'),
						     dendrogram = 'row',margins = c(12,9),labRow = NA, srtCol = 30)


pheatmap( dge.samples.rlog, cluster_rows = T, cluster_cols = F, 
          scale = 'row', clustering_distance_rows = 'correlation' )







# {plot_grid( plotlist = ., 
#                            labels   = c('A','B','C','D'),
#                            ncol     = 2, nrow = 2)}

 
# alternative way, successfully
# arrangeGrob(grobs = wp.GO.plots, nrow = 4, ncol = 1 ) %>% grid.arrange()
#---

#---
# depreacteed, 
# may need register first, otherwise throw out the errors
#
#david.url <- "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/"
#david     <- DAVIDWebService$new(email = "zhenyisong@fuwaihospital.org", url = david.url)
#setAnnotationCategories(david, 'KEGG_PATHWAY')
#liangp.wp1.result.1 %>% 
#filter(P.Value < 0.05 & abs(logFC) > 0.58) %>%
#dplyr::select(GeneID) %>% unlist %>%
#addList(david, .,
#idType="ENTREZ_GENE_ID",
#listName = "clusterProfiler", listType="Gene")
#

#---
# sample-gene correlation analysis
#--
save.image('liangp.Rdata')
quit('no')