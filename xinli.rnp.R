# @author Yisong Zhen
# @since   2017-07-10
# @update  2017-08-22
# nohup R CMD BATCH /home/zhenyisong/biodata/wanglab/wangcode/xinli.rnp.R &
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
           'cluster','factoextra','ggpubr', 'outliers',
           'fitdistrplus', 'VennDiagram', 'vioplot',
           'Rsamtools', 'devtools','parallel',
           'GenomicAlignments', 'BiocParallel',
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
load.lib <- lapply(pkgs, require, character.only = TRUE)

# import data from my X1 desktop
#---
x1.runing.path <- file.path('D:\\wangli_data\\Rdata')
setwd(x1.runing.path)
load('xinliRNP.Rdata')
#load('multiple.xinli.Rdata')

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
xinliRNP.files           <- list.files( path = xinliRNP.rawdata, pattern = 'm_vsmc_.*$', 
                                        all.files = FALSE, full.names  = TRUE, 
                                        recursive = FALSE, ignore.case = FALSE, include.dirs = F)
read.1.files             <- xinliRNP.files[grep('R1',xinliRNP.files)]
read.2.files             <- xinliRNP.files[grep('R2',xinliRNP.files)]

xinliRNP.output.filenames<- basename(read.1.files) %>% 
                            sub(pattern = '_R1_001.fastq.gz', replacement = '') %>%
                            paste0(xinliRNP.output.dir,'/', . ,'.bam')
xinliRNP.sample.names    <- basename(read.1.files) %>% 
                            sub(pattern = '_R1_001.fastq.gz', replacement = '')

sampleFile      <- tempfile(pattern = 'zhen3temp', tmpdir = tempdir())
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
                             aligner = 'Rbowtie',
                             maxHits = 1,
                             paired  = 'fr',
                             splicedAlignment = FALSE,
                             snpFile = NULL,
                             bisulfite = 'no',
                             alignmentParameter = NULL,
                             projectName = 'qProject',
                             alignmentsDir = xinliRNP.QC.dir,
                             lib.loc  = NULL,
                             cacheDir = NULL,
                             clObj = cluster,
                             checkOnly = F)

qQCReport( xinliRNP.qPorject, pdfFilename = 'xinliRNP.QC.pdf', 
           useSampleNames = TRUE, clObj = cluster)
# add this one at 2017-07-18
stopCluster(cluster)

setwd(rsubread.index.lib)
base.string          <-  'mm10'
align( index          = base.string, 
       readfile1      = read.1.files, 
       readfile2      = read.2.files, 
       input_format   = "gzFASTQ", 
       type           = 'rna',
       output_file    = xinliRNP.output.filenames, 
       output_format  = 'BAM',
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
                                       annot.inbuilt = 'mm10', allowMultiOverlap = TRUE)

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


# now read the xinli 772 & FBS data from start
#---

xinli.rawdata        <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/xinli/multiple')
xinli.output.dir     <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/xinli/multiple/rsubread')
xinli.files          <- list.files( path = xinli.rawdata, pattern = '*.fastq$', 
                                    all.files = FALSE, full.names  = TRUE, 
                                    recursive = FALSE, ignore.case = FALSE, include.dirs = F) %>% 
                        {.[c(5:8,16:21)]}

xinli.output.filenames <- basename(xinli.files) %>% sub(pattern = '_R1_001.fastq', replacement = '') %>%
                             paste0(xinli.output.dir ,'/', . ,'.bam') 

setwd(rsubread.index.lib)
base.string           <- 'mm10'
align( index          = base.string, 
       readfile1      = xinli.files, 
       input_format   = 'FASTQ', 
       type           = 'rna',
       output_file    = xinli.output.filenames, 
       output_format  = 'BAM',
       PE_orientation = 'fr', 
       nthreads       = 20, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

xinli.genes     <- featureCounts( xinli.output.filenames, 
                                  useMetaFeatures = TRUE,
                                  countMultiMappingReads = FALSE,
                                  strandSpecific         = 0, 
                                  isPairedEnd            = FALSE,
                                  autosort               = TRUE,
                                  nthreads               = 20,
                                  annot.inbuilt          = 'mm10', 
                                  GTF.featureType        = 'exon',
                                  GTF.attrType           = 'gene_id',
                                  allowMultiOverlap      = TRUE)

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
xinliRIP.rpkm    <- xinliRNP.genes %$% counts %>% DGEList(genes = xinliRNP.genes$annotation) %>%
                    rpkm(normalized.lib.sizes = TRUE, log = FALSE)
colnames(xinliRIP.rpkm) <- c('input','thocs2','thocs5')
xinliRNP.symbols <- mapIds( org.Mm.eg.db, keys = as.character(xinliRNP.genes$annotation$GeneID), 
                            column = 'SYMBOL', keytype = 'ENTREZID', multiVals = 'first') %>% 
                    unlist %>%
                    make.names(unique = T)
xinliRIP.rpkm.df <- xinliRIP.rpkm %>% as.data.frame %>% 
                    mutate(GeneSymbol = xinliRNP.symbols, EntrezID = xinliRNP.genes$annotation$GeneID) %>%
                    dplyr::select(GeneSymbol, EntrezID, input, thocs2, thocs5)


setwd("C:\\Users\\Yisong\\Desktop")
write.xlsx(xinliRIP.rpkm.df , file = 'rip.rpkm.xlsx', colNames = TRUE, borders = 'columns')

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
control.groups.ratio <- control.groups.df %>% mutate( ratio2 = xinliRIP.rpkm.df[control.match.id,'thocs2']/
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
# these two data set were deprecated
#
"
vsmc.contractile.df     <- read_tsv('VSMC_contract_geneTSS1k.bed',  col_names = FALSE)
vsmc.secretion.df       <- read_tsv('VSMC_synthetic_geneTSS1k.bed', col_names = FALSE)
"




"
# second dataset from P1~P10
#---
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

"
vsmc.pdgfDD.sec.df     <- read_tsv('PDGFDD_contract_mus',  col_names = FALSE)
vsmc.pdgfDD.con.df     <- read_tsv('PDGFDD_synthetic_mus', col_names = FALSE)
vsmc.contractile.id    <- unique(vsmc.pdgfDD.con.df$X1) %>%
                          match(Ann$SYMBOL) %>% na.omit
vsmc.secretion.id      <- unique(vsmc.pdgfDD.sec.df$X1) %>%
                          match(Ann$SYMBOL) %>% na.omit




#
vsmc.genes.class    <- c( rep('con',length(vsmc.contractile.id)),
                          rep('sec',length(vsmc.secretion.id)) ) %>%
                       as.factor
# please notice that this is not Entrez GeneID
# this the row index
#---
vsmc.ann.id                <- c(vsmc.contractile.id, vsmc.secretion.id)

xinliRIP.input.rpkm        <- xinliRIP.rpkm[,1]
names(xinliRIP.input.rpkm) <- as.character(1:length(xinliRIP.input.rpkm))
quantile.size              <- 6
vsmc.gene.Q.index          <- paste('Q', 1:quantile.size, sep = '')
xinliRIP.groups            <- xinliRIP.input.rpkm[xinliRIP.input.rpkm > 1] %>%
                              na.omit %>% sort(decreasing = TRUE) %>%
                              split( ceiling(seq_along(.)/(length(.)/quantile.size))) %>%
                              map(names)

vsmc.gene.levels           <- rep('Q0',length(vsmc.ann.id))

for( i in 1:quantile.size) {
    vsmc.gene.levels[vsmc.ann.id %in% as.integer(xinliRIP.groups[[i]]) ] <- vsmc.gene.Q.index[i]
}



"
vsmc.mapping.func          <- function(x, levels, groups, index) {
                                  vsmc.gene.levels  <- levels
                                  xinliRIP.groups   <- groups
                                  vsmc.gene.Q.index <- index
                                  vsmc.gene.levels[vsmc.ann.id %in% xinliRIP.groups[[x]]] <- vsmc.gene.Q.index[x]
                                  return(vsmc.gene.levels)
                              }


map(1:10, vsmc.mapping.func, levels = vsmc.gene.levels,  groups = xinliRIP.groups, index = vsmc.gene.Q.index )
"


vsmc.gene.levels <- as.factor(vsmc.gene.levels)

thocs.affinity      <- xinliRIP.rpkm[vsmc.ann.id,] %>%
                       add(0.25) %>% as_tibble %>%
                       mutate( Ratio2 = thocs2/input, 
                               Ratio5 = thocs5/input,
                               Class  = vsmc.genes.class,
                               Groups = vsmc.gene.levels,
                               Symbol = Ann$SYMBOL[vsmc.ann.id]) %>%
                       filter(!(Groups == 'Q0'))





thocs.table        <- thocs.affinity[sample(1:nrow(thocs.affinity), 10),] %>%
                      tableGrob(rows = NULL)
grid.newpage()
grid.draw(thocs.table)

thocs.affinity %>% filter(Groups == 'Q1' & Class == 'con') %>% nrow

pvalue.matrix <- matrix( data = NA, nrow = quantile.size, ncol = 2)

for( i in 1:quantile.size) {
   pvalue.matrix[i,1] <-  thocs.affinity %>% 
                          filter(Groups == vsmc.gene.Q.index[i]) %>% 
                          wilcox.test(Ratio2 ~ Class, data = .)  %$% p.value
   pvalue.matrix[i,2] <-  thocs.affinity %>% 
                          filter(Groups == vsmc.gene.Q.index[i]) %>% 
                          wilcox.test(Ratio5 ~ Class, data = .)  %$% p.value
}
#dimnames(pvalue.matrix) <- data.frame(row.names = vsmc.gene.Q.index, names = c('Thoc2','Thoc5'))
rownames(pvalue.matrix) <- vsmc.gene.Q.index
colnames(pvalue.matrix) <- c('Thoc2','Thoc5')

thocs.pvalue.table        <- pvalue.matrix %>%
                             {tableGrob(round(., digit = 4), rows = rownames(pvalue.matrix))}
grid.newpage()
grid.draw(thocs.pvalue.table)

rip.thoc2.boxplot   <- ggplot(data = thocs.affinity, aes(x = Class, y = Ratio2)) + 
                       geom_violin(trim = T) +
                       geom_jitter(position = position_jitter(.2) ) +
                       stat_summary( fun.data = 'mean_sdl', fun.args = list(mult = 1),
                                     geom = 'pointrange', color = 'red') +
                       facet_wrap( ~ Groups, scales = 'free', ncol = 2) +
                       xlab('Thoc2 paired-wise non-param comparision \n according to gene expression level') +
                       ylab('Binding RPKM ration: treat/input')
(rip.thoc2.boxplot)

rip.thoc5.boxplot   <- ggplot(data = thocs.affinity, aes(x = Class, y = Ratio5)) + 
                       geom_violin(trim = FALSE) +
                       geom_jitter(position = position_jitter(.2) ) +
                       stat_summary( fun.data = 'mean_sdl', fun.args = list(mult = 1),
                                     geom = 'pointrange', color = 'red') +
                       facet_wrap( ~ Groups, scales = 'free', ncol = 2) +
                       xlab('Thoc5 paired-wise non-param comparision \n according to gene expression level') +
                       ylab('Binding RPKM ration: treat/input')
                      
(rip.thoc5.boxplot)

plot_grid(rip.thoc2.boxplot, rip.thoc5.boxplot, ncol = 1, labels = c('A','B'))




# second trial and error.
#---

xinliRIP.input.rpkm        <- xinliRIP.rpkm[,1]
names(xinliRIP.input.rpkm) <- as.character(1:length(xinliRIP.input.rpkm))
quantile.size              <- 4
vsmc.gene.Q.index          <- paste('Q', 1:quantile.size, sep = '')
xinliRIP.groups            <- xinliRIP.input.rpkm[xinliRIP.input.rpkm > 1] %>%
                              na.omit %>% sort(decreasing = TRUE) %>%
                              split( ceiling(seq_along(.)/(length(.)/quantile.size))) %>%
                              map(names)

vsmc.extra.class                      <- rep('non', nrow(xinliRIP.rpkm))
vsmc.extra.class[vsmc.contractile.id] <- 'con'
vsmc.extra.class[vsmc.secretion.id]   <- 'sec'
vsmc.extra.levels                     <- rep('Q0', nrow(xinliRIP.rpkm))
for( i in 1:quantile.size) {
    vsmc.extra.levels[ as.integer(xinliRIP.groups[[i]]) ] <- vsmc.gene.Q.index[i]
}


thocs.extra      <- xinliRIP.rpkm %>%
                       add(0.25) %>% as_tibble %>%
                       mutate( Ratio2 = thocs2/input, 
                               Ratio5 = thocs5/input,
                               Class  = vsmc.extra.class,
                               Groups = vsmc.extra.levels,
                               Symbol = Ann$SYMBOL) %>%
                       filter(!(Groups == 'Q0'))

pvalue.matrix.non <- matrix( data = NA, nrow = quantile.size, ncol = 2)

for( i in 1:quantile.size) {
   pvalue.matrix.non[i,1] <-  thocs.extra %>% 
                          filter(Class != 'non') %>%
                          filter(Groups == vsmc.gene.Q.index[i]) %>%
                          mutate(Class = as.factor(Class)) %>%
                          wilcox.test(Ratio2 ~ Class, data = .)  %$% p.value
   pvalue.matrix.non[i,2] <-  thocs.extra %>% 
                          filter(Class != 'non' ) %>%
                          filter(Groups == vsmc.gene.Q.index[i]) %>% 
                          mutate(Class = as.factor(Class)) %>%
                          wilcox.test(Ratio5 ~ Class, data = .)  %$% p.value
}
#dimnames(pvalue.matrix) <- data.frame(row.names = vsmc.gene.Q.index, names = c('Thoc2','Thoc5'))
rownames(pvalue.matrix.non) <- vsmc.gene.Q.index
colnames(pvalue.matrix.non) <- c('Thoc2','Thoc5')

thoc.non.pvalue.table        <- pvalue.matrix.non  %>%
                                {tableGrob(round(., digit = 4), rows = rownames(pvalue.matrix.non))}
grid.newpage()
grid.draw(thoc.non.pvalue.table)

rip.thoc2.boxplot   <- thocs.extra %>% 
                       filter((!Class == 'non')) %>%
                       ggplot(data = ., aes(x = Class, y = Ratio2)) + 
                       geom_violin(trim = T) +
                       geom_jitter(position = position_jitter(.2) ) +
                       stat_summary( fun.data = 'mean_sdl', fun.args = list(mult = 1),
                                     geom = 'pointrange', color = 'red') +
                       facet_wrap( ~ Groups, scales = 'free', ncol = 2) +
                       xlab('Thoc2 paired-wise non-param comparision \n according to gene expression level') +
                       ylab('Binding RPKM ration: treat/input')


# third trial and error
# the third set
#---



setwd("E:\\FuWai\\wangli.lab\\vsmc_analysis")
pdgfBB.df <- read_tsv('PDGF_DEG.txt')


# the Venns digram
#---

# see 186
head(xinliRIP.rpkm.df)
percentage                  <- 0.2
xinliRIP.rpkm.sorted.ratio5 <- xinliRIP.rpkm.df %>% mutate(input = input + 0.001) %>%
                               mutate(Ratio5 = thocs5/input) %>% arrange(desc(Ratio5))
gene.num                    <- ( percentage * nrow(xinliRIP.rpkm.sorted.ratio5) ) %>% as.integer

xinliRIP.gene.names         <- xinliRIP.rpkm.sorted.ratio5 %$% GeneSymbol[1:gene.num]

# you have to load the multiple.xinli.Rdata?
# or must use the 772 data set to proceed 
# the data
# direction
# thoc5 - control
# be careful, you must reload the Ann
#
#---

thoc5.772.symbols           <- mapIds( org.Mm.eg.db, keys = as.character(xinli.genes$annotation$GeneID), 
                                       column = 'SYMBOL', keytype = 'ENTREZID', multiVals = 'first') %>% 
                               unlist %>%
                               make.names(unique = T)
cell.lines                  <- factor(rep(c(1,2), 2), levels = 1:2, labels = c('C1','C2'))
groups                      <- factor(c(1,1,2,2), levels = 1:2, labels = c('Control','Thoc5'))
design                      <- model.matrix(~ 0 + groups + cell.lines);
colnames(design)            <- c('Control','Thoc5','Batch')
contrast.matrix             <- makeContrasts(Thoc5 - Control, levels = design)
gene.thoc5.772              <- xinli.genes %$% counts[,5:8] %>%
                               DGEList() %>% calcNormFactors() %>%
                               voom(design = design) %>%
                               lmFit() %>%
                               contrasts.fit(contrast.matrix) %>%
                               eBayes() %>%
                               topTable(number = Inf, adjust.method = 'BH', sort.by = 'none') %>%
                               mutate(symbol = thoc5.772.symbols)
summary(gene.thoc5.772)
filter.lfc                  <- -0.58
filter.pval                 <- 0.05
gene.names.772              <- gene.thoc5.772 %>% 
                               filter( adj.P.Val < filter.pval & logFC <  filter.lfc) %>%
                               dplyr::select(symbol) %>% unlist

# no recylce
# https://stackoverflow.com/questions/29957893/r-convert-vector-to-matrix-without-recycling
#---
commmon.RIP <- intersect(xinliRIP.gene.names, gene.names.772) %>% na.omit
length(commmon.RIP)   <- prod(dim(matrix(commmon.RIP, ncol = 4)))

result.table <- grid.newpage() %>% {matrix(commmon.RIP, ncol= 4, byrow = TRUE)} %>%
                tableGrob(rows = NULL) %>%
                grid.draw()

area.1     <- length(xinliRIP.gene.names %>% na.omit)
area.2     <- length(gene.names.772 %>% na.omit)
cross.area <- length(commmon.RIP)
grid.newpage() %>%  {draw.pairwise.venn( area.1, area.2, cross.area, 
                     alpha = rep(0.5, 2), fill = c('cyan','navyblue'),
                     col = rep('gray97',2), lwd = rep(0,2), 
                     euler.d = T, scaled = F,
                     cex = c(1.5,1.5,1.5),
                     category = c('XinliRIP','722-DGE'),
                     cat.cex = 1.2)}


#
# wangli,assignment
# email title:
# Re: Fw:Fw: PDGFBB_DEG_avgFC,PDGF_DEG
# 08-17,2017
#---
xinliRIP.rpkm <- xinliRNP.genes %$% counts %>% DGEList(genes = xinliRNP.genes$annotation) %>%
                 rpkm(normalized.lib.sizes = TRUE, log = FALSE)
colnames(xinliRIP.rpkm) <- c('input','thoc2','thoc5')

xinliRIP.rpkm.symbols <- mapIds( org.Mm.eg.db, keys= as.character(xinliRNP.genes$annotation$GeneID), 
                                 keytype = "ENTREZID", column = 'SYMBOL') %>% make.names(unique = T)
xinliRIP.rpkm.df <- xinliRIP.rpkm %>% as.data.frame %>% 
                    mutate( GeneSymbol = xinliRIP.rpkm.symbols, 
                            EntrezID   = xinliRNP.genes$annotation$GeneID) %>%
                    dplyr::select(GeneSymbol, EntrezID, input, thoc2, thoc5)


thoc5.symbols               <- mapIds( org.Mm.eg.db, keys= as.character(xinli.genes$annotation$GeneID), 
                                       keytype = "ENTREZID", column = 'SYMBOL') %>% make.names(unique = T)
cell.lines                  <- factor(rep(c(1,2), 2), levels = 1:2, labels = c('C1','C2'))
groups                      <- factor(c(1,1,2,2), levels = 1:2, labels = c('Control','Thoc5'))
design                      <- model.matrix(~ 0 + groups + cell.lines);
colnames(design)            <- c('Control','Thoc5','Batch')
contrast.matrix             <- makeContrasts(Thoc5 - Control, levels = design)
gene.thoc5.772              <- xinli.genes %$% counts[,5:8] %>%
                               DGEList() %>% calcNormFactors() %>%
                               voom(design = design) %>%
                               lmFit() %>%
                               contrasts.fit(contrast.matrix) %>%
                               eBayes() %>%
                               topTable(number = Inf, adjust.method = 'BH', sort.by = 'none') %>%
                               mutate(symbol = thoc5.symbols)

rip.exprs.tidy    <- inner_join( xinliRIP.rpkm.df, gene.thoc5.772, 
                                 by = c('GeneSymbol' = 'symbol')) %>%
                     mutate(Ratio5 = thoc5/(input + 0.0001))

rip.exprs.up      <- rip.exprs.tidy %>% filter((P.Value < 0.05 & logFC > 0)) %>%
                     dplyr::select(Ratio5) %>% mutate(Class = 1)

rip.exprs.down    <- rip.exprs.tidy %>% filter((P.Value < 0.05 & logFC < 0)) %>%
                     dplyr::select(Ratio5) %>% mutate(Class = 2)
rip.exprs.whole   <- rip.exprs.tidy  %>%
                     dplyr::select(Ratio5) %>% 
                     mutate(Class = 3)


# outlier.value     <- outlier(rip.exprs.whole %$% Ratio5 %>% unlist %>% {log(. + 1)})


rip.exprs.class.df <- rbind(rip.exprs.up, rip.exprs.down) %>%
                      rbind(rip.exprs.whole) %>% mutate(Class = as.factor(Class))  


"
comparision.1 <- rip.exprs.class.df %>% filter(Class != 2) %>% 
                 wilcox.test(Ratio5 ~ Class, data = .)  %$% p.value %>%
                 formatC(format = 'g', digits = 2)
 
comparision.2 <- rip.exprs.class.df %>% filter(Class != 1) %>% 
                 wilcox.test(Ratio5 ~ Class, data = .)  %$% p.value %>%
                 formatC(format = 'g', digits = 2)
comparision.3 <- rip.exprs.class.df %>% filter(Class != 3) %>% 
                 wilcox.test(Ratio5 ~ Class, data = .)  %$% p.value %>%
                 formatC(format = 'g', digits = 2)
comparisions  <- p.adjust(c(comparision.1, comparision.2, comparision.3), method = 'bonferroni') %>%
                 formatC(format = 'g', digits = 2)

# find and exclude the extreme outlier
# http://www.itl.nist.gov/div898/handbook/prc/section1/prc16.htm
# maybe this method is argumentable
# without knowing the distribution
# https://stats.stackexchange.com/questions/7155/rigorous-definition-of-an-outlier
#---

IQ.logRatio5        <- rip.exprs.class.df %>%
                       mutate(logRatio5 = log(Ratio5 + 1)) %$%
                       {quantile(logRatio5,3/4) - quantile(logRatio5,1/4)}
                    
Q1.logRatio5        <- rip.exprs.class.df %>%
                       mutate(logRatio5 = log(Ratio5 + 1)) %$%
                       {quantile(logRatio5,1/4)}
Q3.logRatio5        <- rip.exprs.class.df %>%
                       mutate(logRatio5 = log(Ratio5 + 1)) %$%
                       {quantile(logRatio5,3/4)}
upper.outer.fence   <- Q3.logRatio5 + 1.5 * IQ.logRatio5
rip.thoc5.vioplot   <- rip.exprs.class.df %>%
                       mutate(logRatio5 = log(Ratio5 + 1)) %>%
                       filter(logRatio5 < upper.outer.fence) %>%
                       ggplot(aes(x = Class, y = logRatio5)) + 
                       geom_violin(trim = T) +
                       stat_summary( fun.data = 'mean_sdl', fun.args = list(mult = 1),
                                     geom = 'pointrange', color = 'blue') +
                       scale_x_discrete(labels = c('up','down','all')) +
                       xlab('Thoc5 RIP ratio (log-transfromed)\n according to expression change direction') +
                       ylab('Binding affinity(RPKM-log transformed ratio): treat/input') +
                       geom_segment(rip.exprs.class.df, aes(x = 1, y = 4, xend = 2, yend = 4) ) +
                       geom_segment(rip.exprs.class.df, aes(x = 1, y = 3.5, xend = 1, yend = 4) )  +
                       geom_segment(rip.exprs.class.df, aes(x = 2, y = 3.5, xend = 2, yend = 4) )  +
                       geom_text( aes(x = 1.5, y = 4.4), 
                                  label = paste('p','=', comparisions[1], sep = ' '), 
                                  vjust = 'middle')          +
                       geom_segment(rip.exprs.class.df, aes(x = 2, y = 5, xend = 3, yend = 5) ) +
                       geom_segment(rip.exprs.class.df, aes(x = 2, y = 4.5, xend = 2, yend = 5) ) +
                       geom_segment(rip.exprs.class.df, aes(x = 3, y = 4.5, xend = 3, yend = 5) ) +
                       geom_text( aes(x = 2.5, y = 5.4), 
                                  label = paste('p','=', comparisions[2], sep = ' '), 
                                  vjust = 'middle')           +
                       geom_segment(rip.exprs.class.df, aes(x = 1, y = 6, xend = 3, yend = 6) ) +
                       geom_segment(rip.exprs.class.df, aes(x = 1, y = 5.5, xend = 1, yend = 6) ) +
                       geom_segment(rip.exprs.class.df, aes(x = 3, y = 5.5, xend = 3, yend = 6) ) +
                       geom_text( aes(x = 2, y = 6.4), 
                                  label = paste('p','=', comparisions[3], sep = ' '), 
                                  vjust = 'middle')  +
                       theme_classic()



rip.thoc5.boxplot   <- rip.exprs.class.df %>%
                       mutate(logRatio5 = log(Ratio5 + 1)) %>%
                       filter(logRatio5 < upper.outer.fence) %>%
                       ggplot(aes(x = Class, y = logRatio5)) + 
                       geom_boxplot() +
                       scale_x_discrete(labels = c('up','down','all')) +
                       xlab('Thoc5 RIP ratio (log-transfromed)\n according to expression change direction') +
                       ylab('Binding affinity(RPKM-log transformed ratio): treat/input') +
                       geom_segment(rip.exprs.class.df, aes(x = 1, y = 3, xend = 2, yend = 3) ) +
                       geom_segment(rip.exprs.class.df, aes(x = 1, y = 2.8, xend = 1, yend = 3) )  +
                       geom_segment(rip.exprs.class.df, aes(x = 2, y = 2.8, xend = 2, yend = 3) )  +
                       geom_text( aes(x = 1.5, y = 3.4), 
                                  label = paste('p','=', comparisions[1], sep = ' '), 
                                  vjust = 'middle')          +
                       geom_segment(rip.exprs.class.df, aes(x = 2, y = 2.5, xend = 3, yend = 2.5) ) +
                       geom_segment(rip.exprs.class.df, aes(x = 2, y = 2.3, xend = 2, yend = 2.5) ) +
                       geom_segment(rip.exprs.class.df, aes(x = 3, y = 2.3, xend = 3, yend = 2.5) ) +
                       geom_text( aes(x = 2.5, y = 3), 
                                  label = paste('p','=', comparisions[2], sep = ' '), 
                                  vjust = 'middle')           +
                       geom_segment(rip.exprs.class.df, aes(x = 1, y = 4, xend = 3, yend = 4) ) +
                       geom_segment(rip.exprs.class.df, aes(x = 1, y = 3.8, xend = 1, yend = 4) ) +
                       geom_segment(rip.exprs.class.df, aes(x = 3, y = 3.8, xend = 3, yend = 4) ) +
                       geom_text( aes(x = 2, y = 4.4), 
                                  label = paste('p','=', comparisions[3], sep = ' '), 
                                  vjust = 'middle')  +
                       theme_classic()
"                   

# using the vioplot and only compare down vs. up
# 2017-09-04

comparision.final   <- rip.exprs.class.df %>% filter(Class != 3) %>% 
                       wilcox.test(Ratio5 ~ Class, data = .)  %$% p.value %>%
                       formatC(format = 'g', digits = 2)

IQ.logRatio5        <- rip.exprs.class.df %>%
                       filter(Class != 3) %>%
                       mutate(logRatio5 = log(Ratio5 + 1)) %$%
                       {quantile(logRatio5,3/4) - quantile(logRatio5,1/4)}
                    
Q1.logRatio5        <- rip.exprs.class.df %>%
                       filter(Class != 3) %>%
                       mutate(logRatio5 = log(Ratio5 + 1)) %$%
                       {quantile(logRatio5,1/4)}
Q3.logRatio5        <- rip.exprs.class.df %>%
                       filter(Class != 3) %>%
                       mutate(logRatio5 = log(Ratio5 + 1)) %$%
                       {quantile(logRatio5,3/4)}
upper.outer.fence   <- Q3.logRatio5 + 1.5 * IQ.logRatio5

rip.thoc5.vioplot   <- rip.exprs.class.df %>%
                       mutate(logRatio5 = log(Ratio5 + 1)) %>%
                       filter(logRatio5 < upper.outer.fence) %>%
                       filter(Class != 3 ) %>%
                       ggplot(aes(x = Class, y = logRatio5)) + 
                       geom_violin(trim = T, aes(fill = Class), show.legend = F, scale = 'count' ) +
                       stat_summary( fun.data = 'mean_sdl', fun.args = list(mult = 1),
                                     geom = 'pointrange', color = 'blue') +
                       scale_x_discrete(labels = c('up','down')) +
                       xlab('Thoc5 RIP ratio (log-transfromed)\n according to expression change direction') +
                       ylab('Binding affinity(RPKM-log transformed ratio): treat/input') +
                       geom_segment(rip.exprs.class.df, aes(x = 1, y = 2.5, xend = 2, yend = 2.5) ) +
                       geom_segment(rip.exprs.class.df, aes(x = 1, y = 2.3, xend = 1, yend = 2.5) )  +
                       geom_segment(rip.exprs.class.df, aes(x = 2, y = 2.3, xend = 2, yend = 2.5) )  +
                       geom_text( aes(x = 1.5, y = 2.7), 
                                  label = paste('p','=', comparision.final, sep = ' '), 
                                  vjust = 'middle')          +
                       theme_classic()

# traditional methods, deprecated!!!
# boxplot(logRatio5 ~ Class, data = ., outline = F)
#---
"
rip.boxplot   <- rip.exprs.class.df %>%
                 mutate(logRatio5 = log(Ratio5 + 1)) %>%
                 boxplot(logRatio5 ~ Class, data = ., outline = F)
up.vioplot    <- rip.exprs.class.df %>%
                 mutate(logRatio5 = log(Ratio5 + 1)) %>%
                 filter(Class == 1) %>% dplyr::select(logRatio5) %>%
                 unlist() %>% as.vector
down.vioplot  <- rip.exprs.class.df %>%
                 mutate(logRatio5 = log(Ratio5 + 1)) %>%
                 filter(Class == 2) %>% dplyr::select(logRatio5) %>%
                 unlist() %>% as.vector

whole.vioplot <- rip.exprs.class.df %>%
                 mutate(logRatio5 = log(Ratio5 + 1)) %>%
                 filter(Class == 3) %>% dplyr::select(logRatio5) %>%
                 unlist() %>% as.vector

vioplot(up.vioplot, down.vioplot, whole.vioplot)
"

"
1）分析36个合并基因，是否在FBS数据定义的分泌型和收缩型基因中富集；

这个似乎是这样一个模型：
在FBS数据定义的全基因中，随机挑选36个基因，而这36个基因是差异基因（分泌或收缩型基因）的可能？
我们目前得挑法是RIP高的基因和772-DEG的交集，得到的36个基因；
我们的想法是：
这种挑基因的方式，富集了超出预期之外的（分泌、收搜）基因数目，比如，随机挑，可能只有12个基因，预期。

我再想想。
"

# I move the snnipet up
#---

"
xinli.rawdata        <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/xinli/multiple')
xinli.output.dir     <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/xinli/multiple/rsubread')
xinli.files          <- list.files( path = xinli.rawdata, pattern = '*.fastq$', 
                                    all.files = FALSE, full.names  = TRUE, 
                                    recursive = FALSE, ignore.case = FALSE, include.dirs = F) %>% 
                        {.[c(5:8,16:21)]}

xinli.output.filenames <- basename(xinli.files) %>% sub(pattern = '_R1_001.fastq', replacement = '') %>%
                             paste0(xinli.output.dir ,'/', . ,'.bam') 

setwd(rsubread.index.lib)
base.string           <- 'mm10'
align( index          = base.string, 
       readfile1      = xinli.files, 
       input_format   = 'FASTQ', 
       type           = 'rna',
       output_file    = xinli.output.filenames, 
       output_format  = 'BAM',
       PE_orientation = 'fr', 
       nthreads       = 20, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

xinli.genes     <- featureCounts( xinli.output.filenames, 
                                  useMetaFeatures = TRUE,
                                  countMultiMappingReads = FALSE,
                                  strandSpecific         = 0, 
                                  isPairedEnd            = FALSE,
                                  autosort               = TRUE,
                                  nthreads               = 20,
                                  annot.inbuilt          = 'mm10', 
                                  GTF.featureType        = 'exon',
                                  GTF.attrType           = 'gene_id',
                                  allowMultiOverlap      = TRUE)
"

batch.effect          <- factor(c(1,2,1,2), levels = 1:2, labels = c('cell1','cell2'))
groups                <- factor(c(1,1,2,2), levels = 1:2, labels = c('Control','FBS'))
design                <- model.matrix(~ 0 + groups + batch.effect);
colnames(design)      <- c('Control','FBS','Batch')
contrast.matrix       <- makeContrasts(FBS - Control, levels = design)

xinli.FBS             <- xinli.genes %$% counts[,1:4] %>%
                         DGEList() %>% calcNormFactors() %>%
                         voom(design = design) %>%
                         lmFit() %>%
                         contrasts.fit(contrast.matrix) %>%
                         eBayes() %>%
                         topTable(number = Inf, adjust.method = 'BH', sort.by = 'none')


core.36.GeneIDs         <- mapIds( org.Mm.eg.db, keys = commmon.RIP, 
                                 column = 'ENTREZID', keytype = 'SYMBOL', multiVals = 'first') %>% 
                           unlist %>% unique() %>% na.omit
FBS.logFC               <- xinli.FBS$logFC
names(FBS.logFC)        <- xinli.genes$annotation$GeneID
FBS.logFC               <- sort(FBS.logFC, decreasing = TRUE)

core2genes              <- data.frame( diseaseId = rep('coreRIP',length(core.36.GeneIDs)), 
                                       geneId    = core.36.GeneIDs, check.names = TRUE)

core2.GSEA              <- GSEA( FBS.logFC, TERM2GENE = core2genes, 
                                 nPerm = 10^5, maxGSSize = 5000, pvalueCutoff = 1)

result.GSEA.table       <- grid.newpage() %>% {summary(core2.GSEA)} %>% as.data.frame %>%
                           dplyr::select(-(c(core_enrichment, leading_edge, ID))) %>%
                           tableGrob(rows = NULL) %>%
                           grid.draw()

gseaplot(core2.GSEA, geneSetID = 'coreRIP')


# Covergae figure
#---
setwd(xinliRNP.output.dir)
xinliRIP.bam.rename <- list('input','thocs2','thocs5')
xinliRIP.bai        <- paste(xinliRIP.bam.rename,'.bam', sep = '')
cluster             <- makeCluster(4)
system.time( mcmapply( sortBam, 
                       file        = xinliRNP.output.filenames, 
                       destination = xinliRIP.bam.rename,
                       overwrite   = TRUE) )

system.time( mclapply( xinliRIP.bai, indexBam, 
                       overwrite   = TRUE) )

stopCluster(cluster)
# BiocParallel::bpparam()
# restrict the param
# library(BiocParallel)
# register(MulticoreParam(worker = 2))
# MulticoreParam() not supported on Windows, use SnowParam()
#---
fuck <- readGAlignmentPairs(xinliRIP.bai[1])
#---
# deprecated,instead
# I used the parrelle processing procedure
#
#---
#system.time ( map2( xinliRNP.output.filenames, 
#                    xinliRIP.bam.rename,
#                    sortBam, overwrite = TRUE ))
#system.time( map( xinliRIP.bai ,
#                indexBam, overwrite = TRUE ) )
# this is a temporary file folder to
# to store then sorted bam files
# and index file
#---
setwd("D:\\wangli_data\\Rdata\\xinRIP")
# http://www.sthda.com/english/wiki/ggbio-visualize-genomic-data
# http://www.tengfei.name/ggbio/docs/man/tracks.html
# author blog, please see more fine-tune method
#---
mouse.txdb    <- TxDb.Mmusculus.UCSC.mm10.knownGene
columns(Mus.musculus)
xinli.control <- genes(mouse.txdb, filter = list(gene_id = 11459)) 
wh            <- keepSeqlevels(xinli.control, "chr8")

mouse.genes   <- genes(mouse.txdb)

xinli.test.gene <- readGAlignments('thocs5.bam', param = ScanBamParam(which = wh), use.names = T)
autoplot(xinli.test.gene, geom = 'polygon', stat = 'coverage',coverage.col = 'green', fill = 'green', alpha = .2)
ggplot(xinli.test.gene, geom = 'polygon', stat = 'coverage',coverage.col = 'green', fill = 'green', alpha = .2)
gene.model   <- autoplot( Mus.musculus, which = wh, 
                          columns = c("GENENAME", "SYMBOL"), 
                          names.expr = "GENENAME::SYMBOL")
thocs5.cov   <- autoplot( 'thocs5.bam', which = wh) + ylim(0,100) + ylab('thocs5')
input.cov    <- autoplot( 'input.bam', which = wh ) + ylim(0,100) + ylab('input')
tracks(gene.model, thocs5.cov, input.cov, heights = c(1,2,2))

program.environ <- session_info()
setwd(xinliRNP.output.dir)
save.image('xinliRNP.Rdata')
quit('no')