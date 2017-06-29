# @author Yisong Zhen
# @since  2017-06-22
# @update 2017-06-28
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

library(tidyverse)
library(Rsubread)
library(org.Hs.eg.db)
library(annotate)
library(edgeR)
library(limma)
library(DESeq2)
library(magrittr)
library(genefilter)
library(openxlsx)
library(pheatmap)
library(gridExtra)
library(grid)
library(ggrepel)
library(purrr)
library(QuasR)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg38.Rbowtie)
# data transfer 
# in window 10
#---

#x1.runing.path <- file.path('D:\\wangli_data\\Rdata')
#setwd(x1.runing.path)
#load('liangp.Rdata')
#

#--

# read the raw data into the working enviroment
# and classsify the data into read group, 1, 2
#

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
#

# preQC by QuasR
#
sampleFile      <- '/tmp/temp.zhen3'
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
# 
sample.info.liangp  <- data.frame( treat = c( 'CIB1.KO.Basal','CIB1.KO.PMA','TRPC1.KO.Basal','TRPC1.KO.Basal',
                                               'TRPC1.KO.PMA','TRPC1.KO.PMA','WT.H9.Basal','WT.H9.Basal',
                                               'WT.H9.PMA','WT.H9.PMA','CIB1.KO.PMA','CIB1.KO.Basal') )
liangp.exprs.vst    <- gene.counts %>%
                       DESeqDataSetFromMatrix(colData = sample.info.liangp, design = ~ treat) %>%
                       DESeq() %>%
                       varianceStabilizingTransformation() %>% assay()

sds <- rowSds(liangp.exprs.vst)
sh  <- shorth(sds)
liangp.exprs.vst.fitlered <- liangp.exprs.vst[sds > 0.6,]
liangp.pca <- prcomp(t(liangp.exprs.vst.fitlered))

pca.groups  <- c( 'CIB1.KO.Basal-5', 'CIB1.KO.PMA-6',  'TRPC1.KO.Basal-11',
                  'TRPC1.KO.Basal-9','TRPC1.KO.PMA-10','TRPC1.KO.PMA-12',
                  'WT.H9.Basal-3',   'WT.H9.Basal-7',  'WT.H9.PMA-4',
                  'WT.H9.PMA-8',     'CIB1.KO.PMA-2',  'CIB1.KO.Basal-1')

pca.num    <- c(5,6,11,9,10,12,3,7,4,8,2,1)

# deprecated, 
# basic ploting function
"
dev.off()
heatmap( cor(liangp.exprs.vst.fitlered), 
         labRow = pca.groups, labCol = pca.groups)

plot(liangp.pca$x[,c(1,2)])
text( liangp.pca$x[,c(1,2)], liangp.pca$x[,c(1,2)], 
      labels = pca.num , cex = 0.5, pos = 3, adj = c(0,1))
"

liangp.pve    <- liangp.pca$sdev^2/sum(liangp.pca$sdev^2)

pve.liangp.df <- data.frame(variance = liangp.pve, pca = 1:length(liangp.pve))
ggplot(pve.liangp.df) +
        xlab('Principle Component') +
        ylab('Proportion of Variance Explained') +
        scale_x_continuous( breaks = c(1:length(liangp.pve)), 
                            labels = as.character(c(1:length(liangp.pve)), 
                            limits = as.character(c(1:length(liangp.pve))))) +
        geom_point(aes(x = pca, y = variance), size = 3) +
        geom_line(aes(x = pca, y = variance), size = 0.8) +
        scale_linetype_discrete() +
        theme(legend.position="none") +
        theme_classic()

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
pheatmap( cor(liangp.exprs.vst.fitlered), cluster_rows = T, cluster_cols = T,
          color = color.bar, scale = 'none', fontsize_col = 8, 
          fontsize_row = 8, cellwidth = 15, cellheight = 15,
          labels_row = pca.groups, labels_col = pca.groups)

pca.maps <- function(pca.index, pca.df, pca.group) {
     pca.check.df <- data.frame( pcaX = pca.df$x[,1], 
                                 pcaY = pca.df$x[,pca.index])
     graph <- ggplot(data = pca.check.df, aes(x = pcaX, y = pcaY)) + 
              geom_point(col = 'blue', size = 2, show.legend = F) +
              geom_text_repel(label = as.character(pca.group)) +
              theme_classic()
     invisible(return(graph))  
}

#pca.maps(pca.df = liangp.pca,pca.group = pca.num, 3)
# here map can be substitued with lapply
# see purrr
# and other R package and function usage
#---
all.in.one <- map( c(2:5), pca.maps, 
                      pca.df    = liangp.pca,
                      pca.group = pca.num)
grid.arrange(grobs = all.in.one, nrow = 2, ncol = 2)



# limma analysis for DGE
# use pipeline to fast and clearify study
#---
restart.group   <- factor( c( 'CIB1.KO.Basal','CIB1.KO.PMA','TRPC1.KO.Basal','TRPC1.KO.Basal',
                              'TRPC1.KO.PMA','TRPC1.KO.PMA','WT.H9.Basal','WT.H9.Basal',
                              'WT.H9.PMA','WT.H9.PMA','CIB1.KO.PMA','CIB1.KO.Basal'), 
                              levels = c( 'CIB1.KO.Basal','CIB1.KO.PMA','WT.H9.Basal',
                                          'WT.H9.PMA','TRPC1.KO.Basal','TRPC1.KO.PMA'));
#cell.line;
restart.design           <- model.matrix(~ 0 + restart.group)
colnames(restart.design) <- levels(restart.group)
restart.contrast.matrix  <- makeContrasts( WT.H9.Basal - CIB1.KO.Basal,  WT.H9.PMA - CIB1.KO.PMA,
                                           WT.H9.Basal - TRPC1.KO.Basal, WT.H9.PMA - TRPC1.KO.PMA,
                                           WT.H9.PMA   - TRPC1.KO.PMA,   WT.H9.Basal - WT.H9.PMA, 
                                           levels = restart.design)

liangp.results   <- gene.counts %>% DGEList(genes = Ann) %>% calcNormFactors() %>% 
                    voom(design = restart.design) %>%
                    lmFit(restart.design) %>% contrasts.fit(restart.contrast.matrix) %>%
                    eBayes() 

# output the result in Excel file
# using new package
#---
setwd("C:\\Users\\Yisong\\Desktop")
contrast.list <- c('WT.H9.Basal - CIB1.KO.Basal',  'WT.H9.PMA - CIB1.KO.PMA',
                   'WT.H9.Basal - TRPC1.KO.Basal', 'WT.H9.PMA - TRPC1.KO.PMA',
                   'WT.H9.PMA   - TRPC1.KO.PMA',    'WT.H9.Basal - WT.H9.PMA')

liangp.wb <- createWorkbook("Yisong")

for( i in 1:6) {
   addWorksheet(liangp.wb, contrast.list[i])
   liangp.result <-  topTable( liangp.results, coef = i, adjust = "BH", 
                               number = Inf, sort.by = "P")
   writeData(liangp.wb, i, liangp.result)
}
saveWorkbook(liangp.wb, file = "liangp.xlsx", overwrite = TRUE)

save.image('liangp.Rdata')
quit('no')