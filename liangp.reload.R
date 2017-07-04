# @author Yisong Zhen
# @since  2017-06-22
# @update 2017-06-29
#
# @rawdata
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

columns      <- c("ENTREZID","SYMBOL", "GENENAME");
GeneInfo     <- select( org.Hs.eg.db, keys= as.character(gene.ids), 
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

color.bar <- colorRampPalette(c("midnightblue", "grey", "mediumvioletred"))(10) 
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
                     pheatmap( ., cluster_rows = T, cluster_cols = T,
                               color = color.bar, scale = 'none', fontsize_col = 8, 
                               fontsize_row = 8, cellwidth = 20, cellheight = 20,
                               labels_row = pca.p3.groups, labels_col = pca.p3.groups)
p3.sample.heatmap

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
save.image('liangp.Rdata')
quit('no')