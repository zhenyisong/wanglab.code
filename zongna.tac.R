# author Yisong Zhen
# since 2017-07-07
# update
#---

pkgs <- c( 'tidyverse','Rsubread','org.Hs.eg.db',
           'edgeR', 'org.Mm.eg.db',
           'limma', 'DESeq2', 'genefilter','grid',
           'openxlsx','pheatmap','gridExtra','ggrepel',
           'QuasR','annotate','clusterProfiler',
           'cluster','factoextra',"ggpubr",
           'RColorBrewer', 'devtools',
           'BSgenome.Mmusculus.UCSC.mm10',
           'BSgenome.Rnorvegicus.UCSC.rn6',
           'BSgenome.Mmusculus.UCSC.mm10.Rbowtie',
           'BSgenome.Rnorvegicus.UCSC.rn6.Rbowtie')
# install.package(pkgs)
# this will install all necessary packages
# required in this project
#---
#install.lib <- lapply(pkgs, library, character.only = TRUE)
install.lib <- lapply(pkgs, require, character.only = TRUE)


# transfer to local data
#---
x1.runing.path <- file.path('D:\\wangli_data\\Rdata')
setwd(x1.runing.path)
load('zongna.tac.Rdata')

# please check the xinli
# all lib indexes were created in specific Rsubread version
# please make sure before calling featureCounts
#---
rsubread.index.dir  <- file.path('/mnt/date/igenomes/rsubread')

zongna.rawdata      <- file.path('/mnt/date/Sequencing/FastQ/20170705_LF_RNAseq_new')
zongna.QC.dir       <- file.path('/home/zhenyisong/biodata/wanglab/zongna/mouseTac/quasQC')
zongna.output.dir   <- file.path('/home/zhenyisong/biodata/wanglab/zongna/mouseTac/rsubread')
zongna.files        <- list.files( path = zongna.rawdata, pattern = "m\\_z\\_.*gz$", 
                                   all.files = FALSE, full.names = TRUE, 
                                   recursive = FALSE, ignore.case = FALSE, include.dirs = F) 
read.1.files     <- zongna.files[grep("R1",zongna.files)]
read.2.files     <- zongna.files[grep("R2",zongna.files)]
zongna.filenames <- basename(read.1.files) %>% 
                    sub(pattern = '_R1_001.fastq.gz', replacement = '') %>%
                    paste0(zongna.output.dir,'/', . ,'.bam')
sample.names     <- basename(read.1.files) %>% 
                    sub(pattern = '_R1_001.fastq.gz', replacement = '')
sampleFile       <- tempfile(pattern = "zhen3temp", tmpdir = tempdir())
sample.file      <- data.frame( FileName1  = read.1.files,
                                FileName2  = read.2.files, 
                                SampleName = sample.names)
write_tsv(sample.file, path = sampleFile)
genome          <- 'BSgenome.Mmusculus.UCSC.mm10'
#genome          <- 'BSgenome.Rnorvegicus.UCSC.rn6'
cluster         <- makeCluster(18)
dir.create( file.path(zongna.QC.dir), showWarnings = FALSE, recursive = TRUE)
dir.create( file.path(zongna.output.dir), showWarnings = FALSE, recursive = TRUE)
setwd(file.path(zongna.QC.dir))

zongna.qProject <- qAlign( sampleFile,
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
                           alignmentsDir = zongna.QC.dir,
                           lib.loc  = NULL,
                           cacheDir = NULL,
                           clObj = cluster,
                           checkOnly = F)

qQCReport( zongna.qProject, pdfFilename = 'zongnaTac.QC.pdf', 
           useSampleNames = TRUE, clObj = cluster)

setwd(rsubread.index.dir)
base.string          <-  'mm10'
align( index          = base.string, 
       readfile1      = read.1.files, 
       readfile2      = read.2.files, 
       input_format   = "gzFASTQ", 
       type           = 'rna',
       output_file    = zongna.filenames, 
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 15, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )


zongna.tac.genes  <- featureCounts( zongna.filenames, useMetaFeatures = TRUE,
                                    countMultiMappingReads = FALSE,
                                    strandSpecific         = 0, 
                                    isPairedEnd            = TRUE,
                                    requireBothEndsMapped  = TRUE,
                                    autosort               = TRUE,
                                    nthreads               = 15,
                                    annot.inbuilt = "mm10", allowMultiOverlap = TRUE)


# this need to be updated to pipeline
# use the inner_join call?
# which would be much flexible
#---
gene         <- zongna.tac.genes
gene.counts  <- gene$counts
gene.ids     <- gene$annotation$GeneI

columns      <- c("ENTREZID","SYMBOL", "GENENAME")
GeneInfo     <- select( org.Mm.eg.db, keys= as.character(gene.ids), 
                        keytype = "ENTREZID", columns = columns)
m            <- match(gene$annotation$GeneID, GeneInfo$ENTREZID)
Ann          <- cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                       GeneInfo[m, c("SYMBOL", "GENENAME")])

Ann$Chr      <- unlist( lapply(strsplit(Ann$Chr, ";"), 
                        function(x) paste(unique(x), collapse = "|")))
Ann$Chr      <- gsub("chr", "", Ann$Chr)


sample.zongna.tac  <- data.frame( treat = c('sham','sham','sham','tac','tac','tac') )

zongna.tac.vst     <- zongna.tac.genes %$% counts %>%
                      DESeqDataSetFromMatrix(colData = sample.zongna.tac, design = ~ treat) %>%
                      DESeq() %>%
                      varianceStabilizingTransformation() %>% assay()
colnames(zongna.tac.vst) <- c('sham.1','sham.2','sham.3','tac.1','tac.2','tac.3')
sds <- rowSds(zongna.tac.vst)
sh  <- shorth(sds)
(sh)
# [1] 0.07010637
#---

"
zongna.svd          <- zongna.tac.vst %>% subset(sds > 0.1) %>%
                       t() %>% scale() %>% svd

plot(zongna.svd$v[,c(1,2)], ylim = c(-.4, .5),xlim = c(.4,.42))
text( zongna.svd$v[,1], zongna.svd$v[,2], 
      labels = c('sham.1','sham.2','sham.3','tac.1','tac.2','tac.3'), 
      cex = 0.5, pos = 3, adj = c(0,1))
"

# old PCA method
zongna.tac.pca         <- zongna.tac.vst %>% subset(sds > 0.1) %>%
                          t() %>% prcomp
plot(zongna.tac.pca$x[,c(1,2)],xlim = c(-40,45), ylim = c(-35,30))
text( zongna.tac.pca$x[,1], zongna.tac.pca$x[,2], 
      labels = c('sham.1','sham.2','sham.3','tac.1','tac.2','tac.3'), 
      cex = 0.8, pos = 3, adj = c(1,2))


# DGE analysis
#---
zongna.tac.groups            <- factor(rep(c('sham','tac'), each = 3), labels = c('sham','tac'))
zongna.design                <- model.matrix(~ 0 + zongna.tac.groups)
colnames(zongna.design)      <- levels(zongna.tac.groups)
zongna.contrast.matrix       <- makeContrasts(tac - sham, levels = zongna.design)
zongna.tac.result            <- zongna.tac.genes %$% counts %>% DGEList(genes = Ann) %>% 
                                calcNormFactors() %>% 
                                voom(design = zongna.design) %>%
                                lmFit(zongna.design) %>% 
                                contrasts.fit(zongna.contrast.matrix) %>%
                                eBayes() %>%
                                topTable( coef          = 1,
                                          number        = Inf, 
                                          adjust.method = "BH", 
                                          sort.by       = "p")
setwd("C:\\Users\\Yisong\\Desktop")
write.xlsx(zongna.tac.result, file = "zongna.tac.xlsx", colNames = TRUE, borders = "columns")


# post QC check
#---


session_info()
setwd(zongna.output.dir)
save.image('zongna.tac.Rdata')
quit('no')