# @author Yisong Zhen
# @since 2017-06-20
# @parent 
#    see fuwai: vsmc_sample_cor.R
#    see fuwai:   microarry2Exprs.R
#    see wanglab: multiple.xinli.R
# @processed R image
#   load from X1 dir
#   vsmc.Rdata
#   multiple.xinli.Rdata

# nohup R CMD BATCH /home/zhenyisong/biodata/wanglab/wangcode/xinli.sampleCorrelation.R &




# smart R programming
# now !

pkgs <- c( 'tidyverse','Rsubread','org.Hs.eg.db',
           'edgeR', 'org.Mm.eg.db', 'affy',
           'annotate','oligo',
           'rae230a.db', 'rat2302.db', 'hgu133a.db',
           'hgu133a2.db', 'hgu133plus2.db', 'mouse4302.db',
           'pd.mogene.1.0.st.v1', 'pd.mogene.2.0.st',
           'pd.ragene.1.0.st.v1','pd.huex.1.0.st.v2',
           'pd.hugene.1.0.st.v1','mogene20sttranscriptcluster.db',
           'mogene11sttranscriptcluster.db','ragene10sttranscriptcluster.db',
           'huex10sttranscriptcluster.db','hugene10sttranscriptcluster.db',
           'mogene10sttranscriptcluster.db',
           'mogene10sttranscriptcluster.db',
           'limma', 'DESeq2', 'genefilter','grid',
           'openxlsx','pheatmap','gridExtra','ggrepel',
           'QuasR','annotate','clusterProfiler',
           'cluster','factoextra',"ggpubr",
           'RColorBrewer')
# install.package(pkgs)
# this will install all necessary packages
# required in this project
#---
#install.lib <- lapply(pkgs, library, character.only = TRUE)
install.lib <- lapply(pkgs, require, character.only = TRUE)

"
x1.running.path <- file.path('D:\\wangli_data\\Rdata')
setwd(x1.running.path)
load('vsmc.Rdata')
load('multiple.xinli.Rdata')
"

# RNA-seq rawdata reproessing
#---

# GSE35664 SRP010854 
# GSE38056 SRP013262
# GSE44461 SRP018779
# GSE51878 SRP032363
# GSE60642 SRP045701
# GSE60641 SRP045702
# GSE65354 SRP052879

# where I depopsit all the indexed genomes in wanglab
#---
rsubread.index.lib <- file.path('/mnt/date/igenomes/rsubread')
# all reference genomes deposit path
#---
hg38.genome.file   <- file.path('/mnt/date/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa')
mm10.genome.file   <- file.path('/mnt/date/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa')
rn6.genome.file    <- file.path('/mnt/date/igenomes/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa')
rn6.GTF.file       <- file.path('/mnt/date/igenomes/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf')


# project output dir

Rdata.output.dir   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/Rdata')
setwd(Rdata.output.dir)
load('xinli.sampleCorrelation.Rdata')


setwd(rsubread.index.lib)
# sessionInfo()
# Rsubread_1.24.2
# the following code implement the indexing the reference genomes by 
# rsubread
#---
"
base.strings     <- list('hg38','mm10','rn6')
reference.sets   <- list( hg38.genome.file, mm10.genome.file, rn6.genome.file)
map2(base.strings, reference.sets, buildindex)
"

#---
# remapping the genome reads
#---

# GSE35664 
# SE model
# Rat
#---

"
GSE35664.rawdata     <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP010854')
GSE35664.output.dir  <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP010854/results')
GSE35664.files       <- list.files( path = GSE35664.rawdata, pattern = '.*.fastq$', 
                                        all.files = FALSE, full.names  = TRUE, 
                                        recursive = FALSE, ignore.case = FALSE, include.dirs = F)
GSE35664.output.filenames  <- basename(GSE35664.files) %>% 
                              sub(pattern = '.fastq', replacement = '') %>%
                              paste0(GSE35664.output.dir,'/', . ,'.bam')
setwd(rsubread.index.lib)
base.string           <-  'rn6'
align( index          = base.string, 
       readfile1      = GSE35664.files, 
       input_format   = 'FASTQ', 
       type           = 'rna',
       output_file    = GSE35664.output.filenames, 
       output_format  = 'BAM',
       PE_orientation = 'fr', 
       nthreads       = 20, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

GSE35664.genes  <- featureCounts( GSE35664.output.filenames, 
                                  useMetaFeatures = TRUE,
                                  countMultiMappingReads = FALSE,
                                  strandSpecific         = 0, 
                                  isPairedEnd            = FALSE,
                                  autosort               = TRUE,
                                  nthreads               = 20,
                                  annot.ext              = rn6.GTF.file , 
                                  isGTFAnnotationFile    = TRUE,
                                  GTF.featureType        = 'exon',
                                  GTF.attrType           = 'gene_id',
                                  allowMultiOverlap      = TRUE)



#
# GSE38056
# PE
# Rat
# Question: 
# How to split paired end SRA file into 2 correct fastq files
# https://www.biostars.org/p/222122/
#---

GSE38056.rawdata     <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP013262')
GSE38056.output.dir  <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP013262/results')
GSE38056.files       <- list.files( path = GSE38056.rawdata, pattern = 'SRR49845.*.fastq$', 
                                        all.files = FALSE, full.names  = TRUE, 
                                        recursive = FALSE, ignore.case = FALSE, include.dirs = F)

GSE38056.1.files         <- GSE38056.files[grep('_1',GSE38056.files)]
GSE38056.2.files         <- GSE38056.files[grep('_2',GSE38056.files)]

GSE38056.output.filenames <- basename(GSE38056.1.files) %>% sub(pattern = '_1.fastq', replacement = '') %>%
                             paste0(GSE38056.output.dir,'/', . ,'.bam')
setwd(rsubread.index.lib)
base.string           <- 'rn6'
align( index          = base.string, 
       readfile1      = GSE38056.1.files , 
       readfile2      = GSE38056.2.files, 
       input_format   = 'FASTQ', 
       type           = 'rna',
       output_file    = GSE38056.output.filenames, 
       output_format  = 'BAM',
       PE_orientation = 'fr', 
       nthreads       = 20, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

GSE38056.genes  <- featureCounts( GSE38056.output.filenames, 
                                  useMetaFeatures = TRUE,
                                  countMultiMappingReads = FALSE,
                                  strandSpecific         = 0, 
                                  isPairedEnd            = TRUE,
                                  requireBothEndsMapped  = TRUE,
                                  autosort               = TRUE,
                                  nthreads               = 20,
                                  annot.ext              = rn6.GTF.file , 
                                  isGTFAnnotationFile    = TRUE,
                                  GTF.featureType        = 'exon',
                                  GTF.attrType           = 'gene_id',
                                  allowMultiOverlap      = TRUE)



#
# GSE44461 SRP018779
# Homo
# PE model
#---

GSE44461.rawdata     <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP018779')
GSE44461.output.dir  <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP018779/results')
GSE44461.files       <- list.files( path = GSE44461.rawdata, pattern = 'SRR.*.fastq$', 
                                    all.files = FALSE, full.names  = TRUE, 
                                    recursive = FALSE, ignore.case = FALSE, include.dirs = F)

GSE44461.1.files         <- GSE44461.files[grep('_1',GSE44461.files)]
GSE44461.2.files         <- GSE44461.files[grep('_2',GSE44461.files)]

GSE44461.output.filenames <- basename(GSE44461.1.files) %>% sub(pattern = '_1.fastq', replacement = '') %>%
                             paste0(GSE44461.output.dir,'/', . ,'.bam')

setwd(rsubread.index.lib)
base.string           <- 'hg38'
align( index          = base.string, 
       readfile1      = GSE44461.1.files , 
       readfile2      = GSE44461.2.files, 
       input_format   = 'FASTQ', 
       type           = 'rna',
       output_file    = GSE44461.output.filenames, 
       output_format  = 'BAM',
       PE_orientation = 'fr', 
       nthreads       = 20, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

GSE44461.genes  <- featureCounts( GSE44461.output.filenames, 
                                  useMetaFeatures = TRUE,
                                  countMultiMappingReads = FALSE,
                                  strandSpecific         = 0, 
                                  isPairedEnd            = TRUE,
                                  requireBothEndsMapped  = TRUE,
                                  autosort               = TRUE,
                                  nthreads               = 20,
                                  annot.inbuilt          = 'hg38', 
                                  GTF.featureType        = 'exon',
                                  GTF.attrType           = 'gene_id',
                                  allowMultiOverlap      = TRUE)

# GSE51878 SRP032363
# Homo
# SE model
#---

GSE51878.rawdata     <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP032363')
GSE51878.output.dir  <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP032363/results')
GSE51878.files       <- list.files( path = GSE51878.rawdata, pattern = 'SRR.*.fastq$', 
                                    all.files = FALSE, full.names  = TRUE, 
                                    recursive = FALSE, ignore.case = FALSE, include.dirs = F)

GSE51878.output.filenames <- basename(GSE51878.files) %>% sub(pattern = '\\.fastq', replacement = '') %>%
                             paste0(GSE51878.output.dir,'/', . ,'.bam')

setwd(rsubread.index.lib)
base.string           <- 'hg38'
align( index          = base.string, 
       readfile1      = GSE51878.files, 
       input_format   = 'FASTQ', 
       type           = 'rna',
       output_file    = GSE51878.output.filenames, 
       output_format  = 'BAM',
       PE_orientation = 'fr', 
       nthreads       = 20, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

GSE51878.genes  <- featureCounts( GSE51878.output.filenames, 
                                  useMetaFeatures = TRUE,
                                  countMultiMappingReads = FALSE,
                                  strandSpecific         = 0, 
                                  isPairedEnd            = FALSE,
                                  autosort               = TRUE,
                                  nthreads               = 20,
                                  annot.inbuilt          = 'hg38', 
                                  GTF.featureType        = 'exon',
                                  GTF.attrType           = 'gene_id',
                                  allowMultiOverlap      = TRUE)


# GSE60642 SRP045701
# Mouse
# SE model
#---


GSE60642.rawdata     <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP045701')
GSE60642.output.dir  <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP045701/results')
GSE60642.files       <- list.files( path = GSE60642.rawdata, pattern = 'SRR.*.fastq$', 
                                    all.files = FALSE, full.names  = TRUE, 
                                    recursive = FALSE, ignore.case = FALSE, include.dirs = F)

GSE60642.output.filenames <- basename(GSE60642.files) %>% sub(pattern = '\\.fastq', replacement = '') %>%
                             paste0(GSE60642.output.dir,'/', . ,'.bam')

setwd(rsubread.index.lib)
base.string           <- 'mm10'
align( index          = base.string, 
       readfile1      = GSE60642.files, 
       input_format   = 'FASTQ', 
       type           = 'rna',
       output_file    = GSE60642.output.filenames, 
       output_format  = 'BAM',
       PE_orientation = 'fr', 
       nthreads       = 20, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

GSE60642.genes  <- featureCounts( GSE60642.output.filenames, 
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


# GSE60641 SRP045702
# mouse
# SE
#---



GSE60641.rawdata     <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP045702')
GSE60641.output.dir  <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP045702/results')
GSE60641.files       <- list.files( path = GSE60641.rawdata, pattern = 'SRR.*.fastq$', 
                                    all.files = FALSE, full.names  = TRUE, 
                                    recursive = FALSE, ignore.case = FALSE, include.dirs = F)

GSE60641.output.filenames <- basename(GSE60641.files) %>% sub(pattern = '\\.fastq', replacement = '') %>%
                             paste0(GSE60641.output.dir,'/', . ,'.bam')

setwd(rsubread.index.lib)
base.string           <- 'mm10'
align( index          = base.string, 
       readfile1      = GSE60641.files, 
       input_format   = 'FASTQ', 
       type           = 'rna',
       output_file    = GSE60641.output.filenames, 
       output_format  = 'BAM',
       PE_orientation = 'fr', 
       nthreads       = 20, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

GSE60641.genes  <- featureCounts( GSE60641.output.filenames, 
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


#  GSE65354 SRP052879
# Homo
# PE model
#---


GSE65354.rawdata     <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP052879')
GSE65354.output.dir  <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/nextseq/SRP052879/results')
GSE65354.files       <- list.files( path = GSE65354.rawdata, pattern = 'SRR.*.fastq$', 
                                    all.files = FALSE, full.names  = TRUE, 
                                    recursive = FALSE, ignore.case = FALSE, include.dirs = F)

GSE65354.1.files         <- GSE44461.files[grep('_1',GSE65354.files)]
GSE65354.2.files         <- GSE44461.files[grep('_2',GSE65354.files)]

GSE65354.output.filenames <- basename(GSE65354.1.files) %>% sub(pattern = '_1.fastq', replacement = '') %>%
                             paste0(GSE65354.output.dir,'/', . ,'.bam')

setwd(rsubread.index.lib)
base.string           <- 'hg38'
align( index          = base.string, 
       readfile1      = GSE65354.1.files , 
       readfile2      = GSE65354.2.files, 
       input_format   = 'FASTQ', 
       type           = 'rna',
       output_file    = GSE65354.output.filenames, 
       output_format  = 'BAM',
       PE_orientation = 'fr', 
       nthreads       = 20, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

GSE65354.genes  <- featureCounts( GSE65354.output.filenames, 
                                  useMetaFeatures = TRUE,
                                  countMultiMappingReads = FALSE,
                                  strandSpecific         = 0, 
                                  isPairedEnd            = TRUE,
                                  requireBothEndsMapped  = TRUE,
                                  autosort               = TRUE,
                                  nthreads               = 20,
                                  annot.inbuilt          = 'hg38', 
                                  GTF.featureType        = 'exon',
                                  GTF.attrType           = 'gene_id',
                                  allowMultiOverlap      = TRUE)



# wangli data, xinli,
# FBS ,thoc5
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

"

# GSE29955
# Affymetrix Human Genome U133A 2.0 Array
#--

GSE29955.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE29955')
setwd(GSE29955.rawdata)
GSE29955.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE29955.symbols   <- rownames(GSE29955.exprs) %>%
                      {unlist(mget(., hgu133a2SYMBOL, ifnotfound = NA))}
GSE29955.data      <- cbind(GSE29955.symbols, GSE29955.exprs )


# GSE36487
#  Affymetrix Human Genome U133A Array
#---

GSE36487.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE36487')
setwd(GSE36487.rawdata)
GSE36487.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE36487.symbols   <- rownames(GSE36487.exprs) %>%
                      {unlist(mget(., hgu133aSYMBOL, ifnotfound = NA))}
GSE36487.data      <- cbind(GSE36487.symbols, GSE36487.exprs )


# GSE12261
# Affymetrix Human Genome U133 Plus 2.0 Array
# ---

GSE12261.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE12261')
setwd(GSE12261.rawdata)
GSE12261.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE12261.symbols   <- rownames(GSE12261.exprs) %>%
                      {unlist(mget(., hgu133plus2SYMBOL, ifnotfound = NA))}
GSE12261.data      <- cbind(GSE12261.symbols, GSE12261.exprs )


# GSE17543
# Affymetrix Human Genome U133 Plus 2.0 Array
#---
GSE17543.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE17543')
setwd(GSE17543.rawdata)
GSE17543.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE17543.symbols   <- rownames(GSE17543.exprs) %>%
                      {unlist(mget(., hgu133plus2SYMBOL, ifnotfound = NA))}
GSE17543.data      <- cbind(GSE17543.symbols, GSE17543.exprs )


# GSE13791
# Affymetrix Human Genome U133 Plus 2.0 Array
#---

GSE13791.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE13791')
setwd(GSE13791.rawdata)
GSE13791.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE13791.symbols   <- rownames(GSE13791.exprs) %>%
                      {unlist(mget(., hgu133plus2SYMBOL, ifnotfound = NA))}
GSE13791.data      <- cbind(GSE13791.symbols, GSE13791.exprs )


# GSE11367
# Affymetrix Human Genome U133 Plus 2.0 Array
#---
GSE11367.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE11367')
setwd(GSE11367.rawdata)
GSE11367.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE11367.symbols   <- rownames(GSE11367.exprs) %>%
                      {unlist(mget(., hgu133plus2SYMBOL, ifnotfound = NA))}
GSE11367.data      <- cbind(GSE11367.symbols, GSE11367.exprs )


# GSE47744
# Affymetrix Mouse Genome 430 2.0 Array
#---

GSE47744.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE47744')
setwd(GSE47744.rawdata)
GSE47744.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE47744.symbols   <- rownames(GSE47744.exprs) %>%
                      {unlist(mget(., mouse4302SYMBOL, ifnotfound = NA))}
GSE47744.data      <- cbind(GSE47744.symbols, GSE47744.exprs )


# GSE42813
# Affymetrix Mouse Genome 430 2.0 Array
#---

GSE42813.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE42813')
setwd(GSE42813.rawdata)
GSE42813.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE42813.symbols   <- rownames(GSE42813.exprs) %>%
                      {unlist(mget(., mouse4302SYMBOL, ifnotfound = NA))}
GSE42813.data      <- cbind(GSE42813.symbols, GSE42813.exprs )


# GSE60447
#  Affymetrix Mouse Genome 430 2.0
#---


GSE60447.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE60447')
setwd(GSE60447.rawdata)
GSE60447.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE60447.symbols   <- rownames(GSE60447.exprs) %>%
                      {unlist(mget(., mouse4302SYMBOL, ifnotfound = NA))}
GSE60447.data      <- cbind(GSE60447.symbols, GSE60447.exprs )

# GSE19441
# Affymetrix Human Exon 1.0 ST Array
# please see the zongna.microarray.R
#---

GSE19441.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE19441')
setwd(GSE19441.rawdata)
GSE19441.rma       <- list.files( getwd(), pattern = '\\.CEL|\\.CEL\\.gz', 
                                  full.names = TRUE, ignore.case = TRUE) %>% 
                      read.celfiles(celPath = ., pkgname = 'pd.huex.1.0.st.v2') %>%
                      oligo::rma(target = 'core')
GSE19441.symbols          <- featureNames(GSE19441.rma) %>% 
                             {mapIds( huex10sttranscriptcluster.db, 
                                      column = 'SYMBOL', keys = ., keytype = 'PROBEID',
                                      multiVals = 'first')}
GSE19441.exprs            <- exprs(GSE19441.rma)
GSE19441.data             <- cbind(GSE19441.symbols, GSE19441.exprs)


# GSE56819
# Affymetrix Human Gene 1.0 ST Array
#---

GSE56819.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE56819')
setwd(GSE56819.rawdata)
GSE56819.rma       <- list.files( getwd(), pattern = '\\.CEL|\\.CEL\\.gz', 
                                  full.names = TRUE, ignore.case = TRUE) %>% 
                      read.celfiles(celPath = ., pkgname = 'pd.hugene.1.0.st.v1') %>%
                      oligo::rma(target = 'core')
GSE56819.symbols          <- featureNames(GSE56819.rma) %>% 
                             {mapIds( hugene10sttranscriptcluster.db, 
                                      column = 'SYMBOL', keys = ., keytype = 'PROBEID',
                                      multiVals = 'first')}
GSE56819.exprs            <- exprs(GSE56819.rma)
GSE56819.data             <- cbind(GSE56819.symbols, GSE56819.exprs )

# GSE30004
# Affymetrix Human Gene 1.0 ST Array
#---

GSE30004.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE30004')
setwd(GSE30004.rawdata)
GSE30004.rma       <- list.files( getwd(), pattern = '\\.CEL|\\.CEL\\.gz', 
                                  full.names = TRUE, ignore.case = TRUE) %>% 
                      read.celfiles(celPath = ., pkgname = 'pd.hugene.1.0.st.v1') %>%
                      oligo::rma(target = 'core')
GSE30004.symbols          <- featureNames(GSE30004.rma) %>% 
                             {mapIds( hugene10sttranscriptcluster.db, 
                                      column = 'SYMBOL', keys = ., keytype = 'PROBEID',
                                      multiVals = 'first')}
GSE30004.exprs            <- exprs(GSE30004.rma)
GSE30004.data             <- cbind(GSE30004.symbols, GSE30004.exprs )


# GSE17556
# Affymetrix Human Gene 1.0 ST Array
#---

GSE17556.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE17556')
setwd(GSE17556.rawdata)
GSE17556.rma       <- list.files( getwd(), pattern = '\\.CEL|\\.CEL\\.gz', 
                                  full.names = TRUE, ignore.case = TRUE) %>% 
                      read.celfiles(celPath = ., pkgname = 'pd.hugene.1.0.st.v1') %>%
                      oligo::rma(target = 'core')
GSE17556.symbols          <- featureNames(GSE17556.rma) %>% 
                             {mapIds( hugene10sttranscriptcluster.db, 
                                      column = 'SYMBOL', keys = ., keytype = 'PROBEID',
                                      multiVals = 'first')}
GSE17556.exprs            <- exprs(GSE17556.rma)
GSE17556.data             <- cbind(GSE17556.symbols, GSE17556.exprs )


# GSE63425
# Affymetrix Human Gene 1.0 ST Array
#---


GSE63425.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE63425')
setwd(GSE63425.rawdata)
GSE63425.rma       <- list.files( getwd(), pattern = '\\.CEL|\\.CEL\\.gz', 
                                  full.names = TRUE, ignore.case = TRUE) %>% 
                      read.celfiles(celPath = ., pkgname = 'pd.hugene.1.0.st.v1') %>%
                      oligo::rma(target = 'core')
GSE63425.symbols          <- featureNames(GSE63425.rma) %>% 
                             {mapIds( hugene10sttranscriptcluster.db, 
                                      column = 'SYMBOL', keys = ., keytype = 'PROBEID',
                                      multiVals = 'first')}
GSE63425.exprs            <- exprs(GSE63425.rma)
GSE63425.data             <- cbind(GSE63425.symbols, GSE63425.exprs )


# GSE68021
# Affymetrix Human Gene 1.0 ST Array
# ---

GSE68021.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE68021')
setwd(GSE68021.rawdata)
GSE68021.rma       <- list.files( getwd(), pattern = '\\.CEL|\\.CEL\\.gz', 
                                  full.names = TRUE, ignore.case = TRUE) %>% 
                      read.celfiles(celPath = ., pkgname = 'pd.hugene.1.0.st.v1') %>%
                      oligo::rma(target = 'core')
GSE68021.symbols          <- featureNames(GSE68021.rma) %>% 
                             {mapIds( hugene10sttranscriptcluster.db, 
                                      column = 'SYMBOL', keys = ., keytype = 'PROBEID',
                                      multiVals = 'first')}
GSE68021.exprs            <- exprs(GSE68021.rma)
GSE68021.data             <- cbind(GSE68021.symbols, GSE68021.exprs )


# GSE50251
# Affymetrix Mouse Gene 1.0 ST Array
#---


GSE50251.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE50251')
setwd(GSE50251.rawdata)
GSE50251.rma       <- list.files( getwd(), pattern = '\\.CEL|\\.CEL\\.gz', 
                                  full.names = TRUE, ignore.case = TRUE) %>% 
                      read.celfiles(celPath = ., pkgname = 'pd.hugene.1.0.st.v1') %>%
                      oligo::rma(target = 'core')
GSE50251.symbols          <- featureNames(GSE50251.rma) %>% 
                             {mapIds( hugene10sttranscriptcluster.db, 
                                      column = 'SYMBOL', keys = ., keytype = 'PROBEID',
                                      multiVals = 'first')}
GSE50251.exprs            <- exprs(GSE50251.rma)
GSE50251.data             <- cbind(GSE50251.symbols, GSE50251.exprs )


# GSE66538
# Affymetrix Mouse Gene 1.0 ST Array
#---

GSE66538.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE66538')
setwd(GSE66538.rawdata)
GSE66538.rma       <- list.files( getwd(), pattern = '\\.CEL|\\.CEL\\.gz', 
                                  full.names = TRUE, ignore.case = TRUE) %>% 
                      read.celfiles(celPath = ., pkgname = 'pd.mogene.1.0.st.v1') %>%
                      oligo::rma(target = 'core')
GSE66538.symbols          <- featureNames(GSE66538.rma) %>% 
                             {mapIds( mogene10sttranscriptcluster.db, 
                                      column = 'SYMBOL', keys = ., keytype = 'PROBEID',
                                      multiVals = 'first')}
GSE66538.exprs            <- exprs(GSE66538.rma)
GSE66538.data             <- cbind(GSE66538.symbols, GSE66538.exprs )


# GSE66280
#  Affymetrix Mouse Gene 2.0 ST Array
#---

GSE66280.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE66280')
setwd(GSE66280.rawdata)
GSE66280.rma       <- list.files( getwd(), pattern = '\\.CEL|\\.CEL\\.gz', 
                                  full.names = TRUE, ignore.case = TRUE) %>% 
                      read.celfiles(celPath = ., pkgname = 'pd.mogene.2.0.st') %>%
                      oligo::rma(target = 'core')
GSE66280.symbols          <- featureNames(GSE66280.rma) %>% 
                             {mapIds( mogene20sttranscriptcluster.db, 
                                      column = 'SYMBOL', keys = ., keytype = 'PROBEID',
                                      multiVals = 'first')}
GSE66280.exprs            <- exprs(GSE66280.rma)
GSE66280.data             <- cbind(GSE66280.symbols, GSE66280.exprs )

# GSE15713
#  Affymetrix Rat Expression 230A Array
#---


GSE15713.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE15713')
setwd(GSE15713.rawdata)
GSE15713.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE15713.symbols   <- rownames(GSE15713.exprs) %>%
                      {unlist(mget(., rae230aSYMBOL, ifnotfound = NA))}
GSE15713.data      <- cbind(GSE15713.symbols, GSE15713.exprs )


# GSE66624
#  Affymetrix Rat Expression 230A Array
#---

GSE66624.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE66624')
setwd(GSE66624.rawdata)
GSE66624.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE66624.symbols   <- rownames(GSE66624.exprs) %>%
                      {unlist(mget(., rae230aSYMBOL, ifnotfound = NA))}
GSE66624.data      <- cbind(GSE66624.symbols, GSE66624.exprs )


# GSE13594
# Affymetrix Rat Genome 230 2.0 Array
#---

GSE13594.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE13594')
setwd(GSE13594.rawdata)
GSE13594.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE13594.symbols   <- rownames(GSE13594.exprs) %>%
                      {unlist(mget(., rat2302SYMBOL, ifnotfound = NA))}
GSE13594.data      <- cbind(GSE13594.symbols, GSE13594.exprs )


# GSE19106
# Affymetrix Rat Genome 230 2.0 Array
#---
GSE19106.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE19106')
setwd(GSE19106.rawdata)
GSE19106.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE19106.symbols   <- rownames(GSE19106.exprs) %>%
                      {unlist(mget(., rat2302SYMBOL, ifnotfound = NA))}
GSE19106.data      <- cbind(GSE19106.symbols, GSE19106.exprs )


# GSE21573
# Affymetrix Rat Genome 230 2.0 Array
#---

GSE21573.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE21573')
setwd(GSE21573.rawdata)
GSE21573.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE21573.symbols   <- rownames(GSE21573.exprs) %>%
                      {unlist(mget(., rat2302SYMBOL, ifnotfound = NA))}
GSE21573.data      <- cbind(GSE21573.symbols, GSE21573.exprs )

# GSE31080
# Affymetrix Rat Genome 230 2.0 Array
#---

GSE31080.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE31080')
setwd(GSE31080.rawdata)
GSE31080.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE31080.symbols   <- rownames(GSE31080.exprs) %>%
                      {unlist(mget(., rat2302SYMBOL, ifnotfound = NA))}
GSE31080.data      <- cbind(GSE31080.symbols, GSE31080.exprs )


# 
# GSE15841
# Affymetrix Rat Genome 230 2.0 Array
#---

GSE15841.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE15841')
setwd(GSE15841.rawdata)
GSE15841.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE15841.symbols   <- rownames(GSE15841.exprs) %>%
                      {unlist(mget(., rat2302SYMBOL, ifnotfound = NA))}
GSE15841.data      <- cbind(GSE15841.symbols, GSE15841.exprs )


# GSE19909
#  Affymetrix Rat Genome 230 2.0 Array
#---

GSE19909.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE19909')
setwd(GSE19909.rawdata)
GSE19909.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE19909.symbols   <- rownames(GSE19909.exprs) %>%
                      {unlist(mget(., rat2302SYMBOL, ifnotfound = NA))}
GSE19909.data      <- cbind(GSE19909.symbols, GSE19909.exprs )


# GSE21403
# Affymetrix Rat Genome 230 2.0 Array
#---
GSE21403.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/vsmc/affytex/GSE21403')
setwd(GSE21403.rawdata)
GSE21403.exprs     <- ReadAffy() %>% affy::rma() %>% exprs() 
GSE21403.symbols   <- rownames(GSE21403.exprs) %>%
                      {unlist(mget(., rat2302SYMBOL, ifnotfound = NA))}
GSE21403.data      <- cbind(GSE21403.symbols, GSE21403.exprs )



# now start the correlation analysis
#---
GSE29955_1


#---
# GSE35664 SRP010854 
# GSE38056 SRP013262
# GSE44461 SRP018779
# GSE51878 SRP032363
# GSE60642 SRP045701
# GSE60641 SRP045702
# GSE65354 SRP052879
#---
rnaseq.list      <- list( GSE35664.genes, GSE38056.genes, GSE44461.genes,
                          GSE51878.genes, GSE60642.genes, GSE60641.genes,
                          GSE65354.genes )

annot.list       <- map(rnaseq.list, . %$% annotation)
rnaseq.rlog.func <- function(genes,annot) {
                        genes %$% counts %>% 
                        DGEList(genes = annot) %>% 
                        calcNormFactors() %>% rpkm() %>%
                        apply(c(1,2), as.integer) %>% rlog()
                    }

rnaseq.matrix.list <- map2(rnaseq.list, annot.list, rnaseq.rlog.func)


# microarray data analysis
#---



GSE29955.1 <- as.matrix(GSE29955.data[,2:ncol(GSE29955.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(4:6)], 1, median) - 
               apply(.[,c(1:3)], 1, median))}
GSE29955.2 <- as.matrix(GSE29955.data[,2:ncol(GSE29955.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(10:12)], 1, median) - 
               apply(.[,c(7:9)], 1, median))}

GSE29955.3 <- as.matrix(GSE29955.data[,2:ncol(GSE29955.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(16:18)], 1, median) - 
               apply(.[,c(13:15)], 1, median))}

GSE36487.1 <- as.matrix(GSE36487.data[,2:ncol(GSE36487.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,2], 1, median) - 
               apply(.[,1], 1, median))}
GSE36487.2 <- as.matrix(GSE36487.data[,2:ncol(GSE36487.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,3], 1, median) - 
               apply(.[,1], 1, median))}
GSE12261.1 <- as.matrix(GSE12261.data[,2:ncol(GSE12261.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,4:6], 1, median) - 
               apply(.[,1:3], 1, median))}
GSE12261.2 <- as.matrix(GSE12261.data[,2:ncol(GSE12261.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,10:12], 1, median) - 
               apply(.[,7:9], 1, median))}
GSE13791.1 <- as.matrix(GSE13791.data[,2:ncol(GSE13791.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,17:18], 1, median) - 
                apply(.[,19:20], 1, median))}
GSE13791.2 <- as.matrix(GSE13791.data[,2:ncol(GSE13791.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,21:23], 1, median) - 
                apply(.[,24:26], 1, median))} 

GSE11367   <- as.matrix(GSE11367.data[,2:ncol(GSE11367.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(2,4,6)], 1, median) - 
                apply(.[,c(1,3,5)], 1, median))}  

GSE47744.1 <- as.matrix(GSE47744.data[,2:ncol(GSE47744.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(3,4)], 1, median) - 
                apply(.[,c(1,2)], 1, median))}  
GSE47744.2 <- as.matrix(GSE47744.data[,2:ncol(GSE47744.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(7,8)], 1, median) - 
                apply(.[,c(5,6)], 1, median))}  
GSE47744.3 <- as.matrix(GSE47744.data[,2:ncol(GSE47744.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(9,10)], 1, median) - 
                apply(.[,c(5,6)], 1, median))}  
GSE42813.1 <- as.matrix(GSE42813.data[,2:ncol(GSE42813.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(1,2)], 1, median) - 
                apply(.[,c(3,4)], 1, median))} 
 
GSE42813.2 <- as.matrix(GSE42813.data[,2:ncol(GSE42813.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(5,6)], 1, median) - 
                apply(.[,c(3,4)], 1, median))} 
GSE60447.1 <- as.matrix(GSE60447.data[,2:ncol(GSE60447.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(1,2)], 1, median) - 
                apply(.[,c(3,4,5)], 1, median))} 
GSE60447.2 <- as.matrix(GSE60447.data[,2:ncol(GSE60447.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(9:11)], 1, median) - 
                apply(.[,c(6,7,8)], 1, median))} 
GSE19441   <- as.matrix(GSE19441.data[,2:ncol(GSE19441.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(2,4,6)], 1, median) - 
                apply(.[,c(1,3,5)], 1, median))} 
GSE56819.1 <- as.matrix(GSE56819.data[,2:ncol(GSE56819.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(4,5,6)], 1, median) - 
                apply(.[,c(1,2,3)], 1, median))}
GSE56819.2 <- as.matrix(GSE56819.data[,2:ncol(GSE56819.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(7,8,9)], 1, median) - 
                apply(.[,c(1,2,3)], 1, median))} 

GSE17556.1 <- as.matrix(GSE17556.data[,2:ncol(GSE17556.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(4,5,6)], 1, median) - 
                apply(.[,c(1,2,3)], 1, median))} 
GSE17556.2 <- as.matrix(GSE17556.data[,2:ncol(GSE17556.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(7,8,9)], 1, median) - 
                apply(.[,c(1,2,3)], 1, median))} 
GSE17556.3 <- as.matrix(GSE17556.data[,2:ncol(GSE17556.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(13,14,15)], 1, median) - 
                apply(.[,c(10,11,12)], 1, median))} 
GSE17556.4 <- as.matrix(GSE17556.data[,2:ncol(GSE17556.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(13,14,15)], 1, median) - 
                apply(.[,c(16,17,18)], 1, median))} 
GSE63425.1 <- as.matrix(GSE63425.data[,2:ncol(GSE63425.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(1,2,8,12)], 1, median) - 
                apply(.[,c(3,4,7,11)], 1, median))} 
GSE63425.2 <- as.matrix(GSE63425.data[,2:ncol(GSE63425.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(5,6,9,10)], 1, median) - 
                apply(.[,c(3,4,7,11)], 1, median))} 

GSE68021.1 <- as.matrix(GSE68021.data[,2:ncol(GSE68021.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(4,5,6)], 1, median) - 
                apply(.[,c(1,2,3)], 1, median))} 

GSE68021.2 <- as.matrix(GSE68021.data[,2:ncol(GSE68021.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(7,8,9)], 1, median) - 
                apply(.[,c(1,2,3)], 1, median))} 
GSE68021.3 <- as.matrix(GSE68021.data[,2:ncol(GSE68021.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(10,11,12)], 1, median) - 
                apply(.[,c(1,2,3)], 1, median))} 
GSE68021.4 <- as.matrix(GSE68021.data[,2:ncol(GSE68021.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(13,14,15)], 1, median) - 
                apply(.[,c(1,2,3)], 1, median))} 
GSE68021.5 <- as.matrix(GSE68021.data[,2:ncol(GSE68021.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(16,17,18)], 1, median) - 
                apply(.[,c(1,2,3)], 1, median))} 
GSE68021.6 <- as.matrix(GSE68021.data[,2:ncol(GSE68021.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(19,20,21)], 1, median) - 
                apply(.[,c(1,2,3)], 1, median))} 
GSE50251  <- as.matrix(GSE50251.data[,2:ncol(GSE50251.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(1,2,3)], 1, median) - 
                apply(.[,c(4,5,6)], 1, median))} 
GSE66538  <- as.matrix(GSE66538.data[,2:ncol(GSE66538.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,c(5,6,7,8)], 1, median) - 
                apply(.[,c(1,2,3,4)], 1, median))} 
GSE66280  <- as.matrix(GSE66280.data[,2:ncol(GSE66280.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,1:6], 1, median) - 
                apply(.[,7:12], 1, median))} 
GSE15713  <- as.matrix(GSE15713.data[,2:ncol(GSE15713.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,1], 1, median) - 
                apply(.[,2], 1, median))} 
GSE66624  <- as.matrix(GSE66624.data[,2:ncol(GSE66624.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,11:12], 1, median) - 
                apply(.[,5:6], 1, median))} 

GSE13594  <- as.matrix(GSE13594.data[,2:ncol(GSE13594.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,3:4], 1, median))} - 
                apply(.[,1:2], 1, median))} 
GSE19106  <- as.matrix(GSE19106.data[,2:ncol(GSE19106.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,3:4], 1, median))} - 
                apply(.[,1:2], 1, median))} 
GSE21573.1 <- as.matrix(GSE21573.data[,2:ncol(GSE21573.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,3:4], 1, median))} - 
                apply(.[,1:2], 1, median))} 
GSE21573.2 <- as.matrix(GSE21573.data[,2:ncol(GSE21573.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,5:6], 1, median))} - 
                apply(.[,1:2], 1, median))} 
GSE21573.3 <- as.matrix(GSE21573.data[,2:ncol(GSE21573.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,7:8], 1, median))} - 
                apply(.[,1:2], 1, median))} 
GSE21573.4<- as.matrix(GSE21573.data[,2:ncol(GSE21573.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,9:10], 1, median))} - 
                apply(.[,1:2], 1, median))} 
GSE31080.1  <- as.matrix(GSE31080.data[,2:ncol(GSE31080.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,2], 1, median))} - 
                apply(.[,1], 1, median))} 
GSE31080.2  <- as.matrix(GSE31080.data[,2:ncol(GSE31080.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,3], 1, median))} - 
                apply(.[,4], 1, median))} 
GSE15841.1  <- as.matrix(GSE15841.data[,2:ncol(GSE15841.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,3:4], 1, median))} - 
                apply(.[,1:2], 1, median))} 
GSE15841.2  <- as.matrix(GSE15841.data[,2:ncol(GSE15841.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,5:6], 1, median))} - 
                apply(.[,1:2], 1, median))} 
GSE15841.3  <- as.matrix(GSE15841.data[,2:ncol(GSE15841.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,9:10], 1, median))} - 
                apply(.[,7:8], 1, median))} 
GSE15841.4  <- as.matrix(GSE15841.data[,2:ncol(GSE15841.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,11:12], 1, median))} - 
                apply(.[,7:8], 1, median))} 
GSE19909    <- as.matrix(GSE19909.data[,2:ncol(GSE19909.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,4:6], 1, median))} - 
                apply(.[,1:3], 1, median))} 
GSE21403   <- as.matrix(GSE21403.data[,2:ncol(GSE21403.data)]) %>%
              apply(c(1,2),as.numeric) %>% 
              {(apply(.[,3:4], 1, median))} - 
                apply(.[,1:2], 1, median))} 



setwd(Rdata.output.dir)
save.image('xinli.sampleCorrelation.Rdata')
quit('no')






lfc.cutoff <- 0.58
# this code is partially extracted from 
# multiple.xinli.R
#---
gene         <- gene.mouse
gene.counts  <- gene$counts[,12:17]
gene.ids     <- gene$annotation$GeneID


columns     <- c("ENTREZID","SYMBOL","MGI","GENENAME");
GeneInfo    <- select( org.Mm.eg.db, keys= as.character(gene.ids), 
                   keytype = "ENTREZID", columns = columns);
m           <- match(gene$annotation$GeneID, GeneInfo$ENTREZID);
Ann         <- cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                          GeneInfo[m, c("SYMBOL", "MGI","GENENAME")]);

Ann$Chr     <- unlist( lapply(strsplit(Ann$Chr, ";"), 
                    function(x) paste(unique(x), collapse = "|")))
Ann$Chr     <- gsub("chr", "", Ann$Chr)


gene.exprs.772   <- gene.counts[,c(1:4)] %>% DGEList(genes = Ann) %>% calcNormFactors
design           <- model.matrix(~ 0 + groups + cell.lines);
colnames(design) <- c('Control','Thoc5','Batch')
contrast.matrix  <- makeContrasts(Thoc5 - Control, levels = design)

gene.result.772  <- gene.exprs.772 %>% voom(design = design) %>% lmFit() %>%
                    contrasts.fit(contrast.matrix) %>% 
                    eBayes() %>%
                    topTable(  number = Inf)
# deprecated
#gene.names.xinli <- gene.result.772 %$% SYMBOL %>% na.omit %>% toupper
#
sample.info     <- data.frame( treat  = c('Control','Control','FBS','FBS') )
xinli.dds       <- gene.exprs.772 %$% t(t(counts) * samples$norm.factors) %>% 
                   apply(2, as.integer) %>% DESeqDataSetFromMatrix(colData   = sample.info, design = ~ treat) %>%
                   varianceStabilizingTransformation %>%
                   assay
xinli.log       <- gene.exprs.772 %$% t(t(counts) * samples$norm.factors) %>% 
                   apply(2, as.integer) %>% add(1) %>% log

xinli.data.filter      <- (apply(xinli.log[,1:2], 1, median) - apply(xinli.log[,3:4], 1, median))
xinli.names            <- Ann %$% SYMBOL %>% subset(abs(xinli.data.filter) > lfc.cutoff) %>%
                          na.omit %>% toupper

RNA.seq.commonnames.xinli   <- xinli.names %>% union(GSE35664_1.name) %>% union(GSE35664_2.name) %>%
                               union(GSE35664_3.name) %>% union(GSE38056.name) %>% union(GSE44461.name) %>%
                               union(GSE51878_1.name) %>% union(GSE51878_2.name) %>% union(GSE60642_1.name) %>%
                               union(GSE60642_2.name) %>% union(GSE60641.name) %>% union(GSE65354_1.name) %>%
                               union(GSE65354_2.name) %>% na.omit

names(xinli.data.filter)    <- make.names(Ann %$% SYMBOL %>% toupper, unique = TRUE)

common.names.xinli          <- intersect(rownames(affy.kd.matrix), RNA.seq.commonnames.xinli)


# I used the vector() and matrix()
# to create the empty matrix /vector
# but it is not applied in this case
#---
xinli.whole.matrix          <- NULL
xinli.rownames              <- NULL


for (gene.name in common.names.xinli) {
    if ( gene.name %in% toupper(Ann$SYMBOL) &
             gene.name %in% toupper(rownames(rna.matrix)) ) {
        row.buffer         <- c(affy.kd.matrix[gene.name,], rna.matrix[gene.name,], xinli.data.filter[gene.name])
        xinli.whole.matrix <- rbind(xinli.whole.matrix, row.buffer)
        xinli.rownames     <- append(xinli.rownames, gene.name)
    }
}

xinli.dist.final <- cor(xinli.whole.matrix ,method = 'spearman')

heatmap.2(xinli.dist.final, scale  = 'none', 
						   Rowv = F,Colv = F, density.info = 'none',
                           key  = TRUE, trace='none', symm = T,symkey = F,symbreaks = T,
						   margins = c(5.5,5.5),dendrogram = 'none',
                           cexRow  = 0.8, cexCol = 0.8, srtCol = 30,
                           labRow  = rownames(xinli.whole.matrix)
						  );