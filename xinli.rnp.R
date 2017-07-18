# @author Yisong Zhen
# @since 2017-07-10
# in response to wang/yuan request
# this old paper
#
# PMID: 28473716
# In this method section
# RNA immunoprecipitation, sequencing and analysis.
# RIP not RNP
#---

# positve markers(controls)
# mouse
# EntrezID Symbol Full.Name
# 21345 Tagln  Sm22a
# 12797 Cnn1   calponin 
# 11459 Acta1  actin, alpha 1, skeletal muscle chr8
# 17880 Myh11  myosin, heavy polypeptide 11, smooth muscle
# ====
# negative control
# 21956 Tnnt2   cTnT
# 21954 Tnni3   troponin I, cardiac 3 chr7
# 230899 Nppa   natriuretic peptide type A /Anp
# 14433 Gapdh   glyceraldehyde-3-phosphate dehydrogenase
# 22142 Tuba1a  tubulin, alpha 1A
# 11461 Actb    actin, beta


pkgs <- c( 'tidyverse','Rsubread','org.Mm.eg.db','edgeR',
           'limma', 'DESeq2', 'genefilter','grid',
           'openxlsx','pheatmap','gridExtra','ggrepel',
           'QuasR','annotate','clusterProfiler',
           'GGally','RColorBrewer','ggbio',
           'cluster','factoextra',"ggpubr",
           'Rsamtools', 'devtools',
           'TxDb.Mmusculus.UCSC.mm10.knownGene',
           'Mus.musculus', 'biovizBase',
           'BSgenome.Hsapiens.UCSC.hg38',
           'BSgenome.Hsapiens.UCSC.hg38.Rbowtie',
           'BSgenome.Mmusculus.UCSC.mm10',
           'BSgenome.Mmusculus.UCSC.mm10.Rbowtie',
           'BSgenome.Rnorvegicus.UCSC.rn6',
           'BSgenome.Rnorvegicus.UCSC.rn6.Rbowtie')
#install.package(pkgs)
# this will install all necessary packages
# required in this project
#---
#install.lib <- lapply(pkgs, library, character.only = TRUE)
install.lib <- lapply(pkgs, require, character.only = TRUE)

# import data from my X1 desktop
#---
x1.runing.path <- file.path('D:\\wangli_data\\Rdata')
setwd(x1.runing.path)
load('xinliRNP.Rdata')

# copycat from
#---
rsubread.index.lib <- file.path('/mnt/date/igenomes/rsubread')
# all reference genomes deposit path
#---
hg38.genome.file   <- file.path('/mnt/date/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa')
mm10.genome.file   <- file.path('/mnt/date/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa')
rn6.genome.file    <- file.path('/mnt/date/igenomes/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa')

xinliRNP.output.dir     <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/xinli/rsubread')
xinliRNP.QC.dir         <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/xinli/quasR')

xinliRNP.rawdata         <- file.path('/mnt/date/Sequencing/FastQ/20170705_LF_RNAseq_new')
xinliRNP.files           <- list.files( path = xinliRNP.rawdata, pattern = "m_vsmc_.*$", 
                                        all.files = FALSE, full.names  = TRUE, 
                                        recursive = FALSE, ignore.case = FALSE, include.dirs = F)
read.1.files             <- xinliRNP.files[grep("R1",xinliRNP.files)]
read.2.files             <- xinliRNP.files[grep("R2",xinliRNP.files)]

xinliRNP.output.filenames<- basename(read.1.files) %>% 
                            sub(pattern = '_R1_001.fastq.gz', replacement = '') %>%
                            paste0(xinliRNP.output.dir,'/', . ,'.bam')
xinliRNP.sample.names    <- basename(read.1.files) %>% 
                            sub(pattern = '_R1_001.fastq.gz', replacement = '')

sampleFile      <- tempfile(pattern = "zhen3temp", tmpdir = tempdir())
sample.file     <- data.frame( FileName1  = read.1.files,
                               FileName2  = read.2.files, 
                               SampleName = xinliRNP.sample.names)
write_tsv(sample.file, path = sampleFile)
genome          <- 'BSgenome.Mmusculus.UCSC.mm10'
cluster         <- makeCluster(18)

dir.create(xinliRNP.QC.dir, showWarnings = FALSE, recursive = TRUE)
dir.create(xinliRNP.output.dir, showWarnings = FALSE, recursive = TRUE)
setwd(xinliRNP.QC.dir)

xinliRNP.qPorject <- qAlign( sampleFile,
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
                             alignmentsDir = xinliRNP.QC.dir,
                             lib.loc  = NULL,
                             cacheDir = NULL,
                             clObj = cluster,
                             checkOnly = F)

qQCReport( xinliRNP.qPorject, pdfFilename = 'xinliRNP.QC.pdf', 
           useSampleNames = TRUE, clObj = cluster)

setwd(rsubread.index.lib)
base.string          <-  'mm10'
align( index          = base.string, 
       readfile1      = read.1.files, 
       readfile2      = read.2.files, 
       input_format   = "gzFASTQ", 
       type           = 'rna',
       output_file    = xinliRNP.output.filenames, 
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 20, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

xinliRNP.genes      <- featureCounts( xinliRNP.output.filenames, useMetaFeatures = TRUE,
                                       countMultiMappingReads = FALSE,
                                       strandSpecific         = 0, 
                                       isPairedEnd            = TRUE,
                                       requireBothEndsMapped  = TRUE,
                                       autosort               = TRUE,
                                       nthreads               = 20,
                                       annot.inbuilt = "mm10", allowMultiOverlap = TRUE)

#-- caculate the and output the result
#-- you mush restart here every time you
#-- replicate your results
#---
gene         <- xinliRNP.genes
gene.counts  <- gene$counts
gene.ids     <- gene$annotation$GeneID


# this need to be updated to pipeline
# use the inner_join call?
# which would be much flexible
#---

columns      <- c("ENTREZID","SYMBOL", "GENENAME");
GeneInfo     <- select( org.Mm.eg.db, keys= as.character(gene.ids), 
                        keytype = "ENTREZID", columns = columns);
m            <- match(gene$annotation$GeneID, GeneInfo$ENTREZID);
Ann          <- cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                       GeneInfo[m, c("SYMBOL", "GENENAME")]);

Ann$Chr      <- unlist( lapply(strsplit(Ann$Chr, ";"), 
                        function(x) paste(unique(x), collapse = "|")))
Ann$Chr      <- gsub("chr", "", Ann$Chr)

# this the end to get the Ann, the gene annotation for this 
# batch of RNA-seq data analysis
#---


# in responding to wangli analysis protocol
# caculating RPKM
# wangli intend to analysis the between-difference
# thocs5 ~ thocs2
# I normalized the raw RPKM data
# I assume that all RNA binding read counts is matched
# total number is the same
#---
xinliRIP.rpkm <- xinliRNP.genes %$% counts %>% DGEList(genes = Ann) %>%
                 rpkm(normalized.lib.sizes = TRUE, log = FALSE)
colnames(xinliRIP.rpkm) <- c('input','thocs2','thocs5')

xinliRIP.rpkm.df <- xinliRIP.rpkm %>% as.data.frame %>% 
                    mutate(GeneSymbol = Ann$SYMBOL, EntrezID = Ann$GeneID) %>%
                    dplyr::select(GeneSymbol, EntrezID, input, thocs2, thocs5)


setwd("C:\\Users\\Yisong\\Desktop")
write.xlsx(xinliRIP.rpkm.df , file = "rip.rpkm.xlsx", colNames = TRUE, borders = "columns")

# QC check for her positive and negative control
#---
control.groups.df <- data.frame( gene.name = c('Tagln','Cnn1','Acta1', 'Myh11', 
                                               'Tnnt2', 'Tnni3', 'Nppa', 'Gapdh','Tuba1a','Actb'),
                                 entrez.id = c(21345,12797,11459,17880,
                                               21956,21954,230899,14433,22142,11461),
                                 class     = c(rep('pos',4),rep('neg',6)),
                                 chr.loc   = c('chr9','chr9','chr8','chr16',
                                               'chr1','chr7','chr4','chr6','chr15','chr5') )
control.match.id     <- match( as.character(control.groups.df$gene.name), 
                               xinliRIP.rpkm.df$GeneSymbol )
control.groups.ratio <- control.groups.df %>% mutate(ratio2 = xinliRIP.rpkm.df[control.match.id,'thocs2']/
                                                              xinliRIP.rpkm.df[control.match.id,'input'],
                                                     ratio5 = xinliRIP.rpkm.df[control.match.id,'thocs5']/
                                                              xinliRIP.rpkm.df[control.match.id,'input'])

control.ratio.table  <- tableGrob(control.groups.ratio, rows = NULL)
grid.newpage()
grid.draw(control.ratio.table)

# the gene name for secretory gene type                       
# and genename for contractile gene list 
# was curtousy from Yupeng's curation
# see his and wangli feedback email
#---

# read the curationdata
# 
setwd("E:\\FuWai\\wangli.lab\\vsmc_analysis")
vsmc.contractile.df     <- read_tsv('VSMC_contract_geneTSS1k.bed',  col_names = FALSE)
vsmc.secretion.df       <- read_tsv('VSMC_synthetic_geneTSS1k.bed', col_names = FALSE)



vsmc.dge.table          <- read_tsv('P1_P10_DEG_2v2.xls',  col_names = TRUE) %>%
                           filter(P.Value < 0.05 & logFC > 0.58) %>% 
                           arrange(desc(logFC))
percentage              <- 0.2
dge.num                 <- {nrow(vsmc.dge.table) * percentage} %>% as.integer
vsmc.contractile.id <- vsmc.dge.table[1:dge.num,1] %>% 
                       unlist %>% unique %>%
                       match(Ann$SYMBOL) %>% na.omit
vsmc.secretion.id   <- vsmc.dge.table[(nrow(vsmc.dge.table) - dge.num):nrow(vsmc.dge.table),1] %>%
                       unlist %>% unique %>%
                       match(Ann$SYMBOL) %>% na.omit


# data set from 
vsmc.contractile.id <- unique(unlist(vsmc.contractile.df[,4])) %>%
                       match(Ann$SYMBOL) %>% na.omit
vsmc.secretion.id   <- unique(unlist(vsmc.secretion.df[,4])) %>%
                       match(Ann$SYMBOL) %>% na.omit

#
vsmc.genes.class    <- c( rep('con',length(vsmc.contractile.id)),
                          rep('sec',length(vsmc.secretion.id)) ) %>%
                       as.factor
# please notice that this is not Entrez GeneID
# this the row index
#---
vsmc.ann.id         <- c(vsmc.contractile.id, vsmc.secretion.id)

xinliRIP.input.rpkm        <- xinliRIP.rpkm[,1]
names(xinliRIP.input.rpkm) <- as.character(1:length(xinliRIP.input.rpkm))
quantile.size              <- 4
xinliRIP.groups            <- xinliRIP.input.rpkm[xinliRIP.input.rpkm > 1] %>%
                              na.omit %>% sort(decreasing = TRUE) %>%
                              split( ceiling(seq_along(.)/(length(.)/quantile.size))) %>%
                              map(names)

vsmc.gene.levels           <- rep('Q1',length(vsmc.ann.id))

vsmc.gene.levels[vsmc.ann.id %in% xinliRIP.groups[[1]]] <- 'Q1'
vsmc.gene.levels[vsmc.ann.id %in% xinliRIP.groups[[2]]] <- 'Q2'
vsmc.gene.levels[vsmc.ann.id %in% xinliRIP.groups[[3]]] <- 'Q3'
vsmc.gene.levels[vsmc.ann.id %in% xinliRIP.groups[[4]]] <- 'Q4'

vsmc.gene.levels <- as.factor(vsmc.gene.levels)

thocs.affinity      <- xinliRIP.rpkm[vsmc.ann.id,] %>%
                       add(0.25) %>% as.data.frame %>%
                       mutate( Ratio2 = thocs2/input, 
                               Ratio5 = thocs5/input,
                               Class  = vsmc.genes.class,
                               Groups = vsmc.gene.levels,
                               Symbol = Ann$SYMBOL[vsmc.ann.id])

thocs.affinity %>% filter(Groups == 'Q1' & Class == 'con') %>% nrow

thocs.affinity %>% filter(Groups == 'Q1') %>% wilcox.test(Ratio2 ~ Class, data = .)
thocs.affinity %>% filter(Groups == 'Q1') %>% wilcox.test(Ratio5 ~ Class, data = .)

rip.thoc2.boxplot   <- ggplot(data = thocs.affinity, aes(x = Class, y = Ratio2)) + 
                       geom_violin(trim = T) +
                       geom_jitter(position = position_jitter(.2) ) +
                       stat_summary( fun.data = 'mean_sdl', fun.args = list(mult = 1),
                                     geom = 'pointrange', color = 'red') +
                       facet_wrap( ~ Groups, scales = 'free', ncol = 2) +
                       xlab('Thocs2 paired-wise non-param comparision \n according to gene expression level') +
                       ylab('Binding RPKM ration: treat/input')
(rip.thoc2.boxplot)

rip.thoc5.boxplot   <- ggplot(data = thocs.affinity, aes(x = Class, y = Ratio5)) + 
                       geom_violin(trim = FALSE) +
                       geom_jitter(position = position_jitter(.2) ) +
                       stat_summary( fun.data = 'mean_sdl', fun.args = list(mult = 1),
                                     geom = 'pointrange', color = 'red') +
                       facet_wrap( ~ Groups, scales = 'free', ncol = 2) +
                       xlab('Thocs5 paired-wise non-param comparision \n according to gene expression level') +
                       ylab('Binding RPKM ration: treat/input')
                      
(rip.thoc5.boxplot)



# Covergae figure
#---
setwd(xinliRNP.output.dir)
xinliRIP.bam.rename <- list('input','thocs2','thocs5')
xinliRIP.bai        <- paste(xinliRIP.bam.rename,'.bam', sep = '')
map2( xinliRNP.output.filenames, 
      xinliRIP.bam.rename,
      sortBam, overwrite = TRUE )
map( xinliRIP.bai ,
     indexBam, overwrite = TRUE )
# this is a temporary file folder to
# to store then sorted bam files
# and index file
#---
setwd("D:\\wangli_data\\Rdata\\xinRIP")
# http://www.sthda.com/english/wiki/ggbio-visualize-genomic-data
# http://www.tengfei.name/ggbio/docs/man/tracks.html
# author blog, please see more fine-tune method
#---
mouse.txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
columns(Mus.musculus)
xinli.control <- genes(mouse.txdb, filter = list(gene_id = 11459)) 
wh            <- keepSeqlevels(xinli.control, "chr8")

gene.model   <- autoplot( Mus.musculus, which = wh, 
                          columns = c("GENENAME", "SYMBOL"), 
                          names.expr = "GENENAME::SYMBOL")
thocs5.cov   <- autoplot( 'thocs5.bam', which = wh) + ylim(0,100) + ylab('thocs5')
input.cov    <- autoplot( 'input.bam', which = wh ) + ylim(0,100) + ylab('input')
tracks(gene.model, thocs5.cov, input.cov, heights = c(1,2,2))

session_info()
setwd(xinliRNP.output.dir)
save.image('xinliRNP.Rdata')
quit('no')