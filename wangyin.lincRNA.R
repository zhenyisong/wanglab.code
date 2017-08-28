# bcl2fastq.sh  WangYin-Total_RNAseq-20180801.csv wangyinLincR
# see lnRNA_next500.sh how to demultiplex samples
#
# nohup R CMD BATCH /home/zhenyisong/biodata/wanglab/wangcode/wangyin.lincRNA.R &
# raw data demultiplexing protocol
# the original sample sheet was sent by wangyin by his recent email
# Wang Yin <wangyin_fuwai@163.com>
# Re:Re: WY-20170801-total
#---
#

# please refer the bcl2fastq manual for parameter meaning
# 
#bcl2fastq --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions \
#--no-lane-splitting  \
#--runfolder-dir /mnt/date/Sequencing/RawData/170801_NB501912_0028_AHF53CBGX3 \
#--sample-sheet /mnt/date/Sequencing/FastQ/sampleSheet/170801_WY_total.csv  \
#--output-dir /mnt/date/Sequencing/FastQ/170801wangyinLincR

#
# I did not specify the threads, see manual page 17, using system default;
#  --processing-threads
#  --demultiplexingthreads
#--- 

# library loading
#---
pkgs         <- c( 'tidyverse','Rsubread', 'tools', 'QuasR','DiagrammeR',
                   'BSgenome.Hsapiens.UCSC.hg38')
install.lib  <- lapply(pkgs, require, character.only = TRUE)
#---end


# rankdir = LR ??
graph <-
    create_graph() %>%
    set_graph_name("miRNA quantification pipeline") %>%
    set_global_graph_attrs("graph", "layout", "dot") %>%
    set_global_graph_attrs("graph", "rankdir", "LR") %>%
    set_global_graph_attrs("node", "color", "white") %>%
    set_global_graph_attrs("node", "fontname", "Arial") %>%
    add_n_nodes(11) %>%
    select_nodes_by_id(c(1:11)) %>% 
    set_node_attrs_ws("shape", "retangle") %>%
    set_node_attrs_ws("style", "filled") %>%
    clear_selection %>%
    add_edges_w_string(
      '1->2 2->3 2->6 2->9 3->4 4->5 6->7 7->8 9->10 10->11', 'black') %>%
    set_node_attrs('label',c( 'H1.d1',  'mesoderm.d3', 
                              'CM.d5',  'CM.d9',    'CM.12d',
                              'EC.d5',  'EC.d9',    'EC.12d',
                              'VSMC.5d','VSMC.9d',  'VSMC.12d')) %>%
    set_node_attrs('fontsize',10) %>%
    set_edge_attrs('arrowsize', 1)
render_graph(graph)
#---
# predifined data and annotation files
#---

hg38.genome.file   <- file.path('/mnt/date/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa')
mm10.genome.file   <- file.path('/mnt/date/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa')
rn6.genome.file    <- file.path('/mnt/date/igenomes/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa')
rn6.GTF.file       <- file.path('/mnt/date/igenomes/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf')
Rdata.output.dir   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/Rdata')

setwd(Rdata.output.dir)
load('wangyin170801.Rdata')

"
x1.running.path <- file.path('D:\\wangli_data\\Rdata')
setwd(x1.running.path)
load('wangyin170801.Rdata')
"

# where I depopsit all the indexed genomes in wanglab
#---
rsubread.index.lib <- file.path('/mnt/date/igenomes/rsubread')

# cd /mnt/date/genomelib/annotation
# NONCODE2016_human_hg38_lncRNA.gtf was dowloaded from noncode database, Zhao, Yi;
# cp NONCODE2016_human_hg38_lncRNA.gtf hg38_UCSC_lncRNA_mRNA.gtf
# cat /mnt/date/igenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf >> hg38_UCSC_lncRNA_mRNA.gtf
# md5sum  hg38_UCSC_lncRNA_mRNA.gtf
# e5098f82ebb9cf4b8c99d357ea43ffb9  hg38_UCSC_lncRNA_mRNA.gtf
#---


lncRmRNA.hg38.GTF.file <- file.path('/mnt/date/genomelib/annotation/g38_UCSC_lncRNA_mRNA.gtf')

rawdata.path    <- '/mnt/date/Sequencing/FastQ/170801wangyinLincR'
wang.output.dir <- '/home/zhenyisong/biodata/wanglab/wangdata/wangyin/170801resuls/rsubread'
wang.QC.dir     <- '/home/zhenyisong/biodata/wanglab/wangdata/wangyin/170801resuls/quasR'
dir.create(file.path(wang.output.dir), showWarnings = FALSE)
dir.create(file.path(wang.QC.dir), showWarnings = FALSE)
setwd(rawdata.path)
wangyin.lincR.files       <- list.files( path = rawdata.path, pattern = '.fastq.gz$', 
                                         all.files = FALSE, full.names  = TRUE, 
                                         recursive = FALSE, ignore.case = FALSE, include.dirs = F) %>%
                             {.[-c(33:34)]}

# first.check  <- md5sum(wangyin.lincR.files)
# second.check <- md5sum(wangyin.lincR.files)
# all.equal(first.check,second.check)
# identical(first.check,second.check)
# save the finger.prints
#---

wang.1.files          <- wangyin.lincR.files[grep('R1_001.fastq.gz',wangyin.lincR.files)]
wang.2.files          <- wangyin.lincR.files[grep('R2_001.fastq.gz',wangyin.lincR.files)]



sample.names.quasR    <- basename(wang.1.files) %>% 
                         sub(pattern = '_R1_001.fastq.gz', replacement = '')
wang.output.files     <- sample.names.quasR %>%
                         paste0(wang.output.dir,'/', . ,'.bam')


# start QC module
# from the liangp.reload.R
#---

sampleFile      <- tempfile(pattern = 'zhen3temp', tmpdir = tempdir())
sample.file     <- data.frame( FileName1  = unlist(wang.1.files),
                               FileName2  = unlist(wang.2.files), 
                               SampleName = unlist(sample.names.quasR))
write_tsv(sample.file, path = sampleFile)
genome          <- 'BSgenome.Hsapiens.UCSC.hg38'
cluster         <- makeCluster(18)

setwd(file.path(wang.QC.dir))

wang.qPorject <- qAlign(   sampleFile,
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
                           alignmentsDir = wang.QC.dir,
                           lib.loc  = NULL,
                           cacheDir = NULL,
                           clObj = cluster,
                           checkOnly = F) %>%
                  qQCReport( pdfFilename = 'wangyin.QC.pdf', 
                             useSampleNames = TRUE, clObj = cluster)
                           
stopCluster(cluster)

setwd(rsubread.index.lib)
base.string           <- 'hg38'
align( index          = base.string, 
       readfile1      = wang.1.files , 
       readfile2      = wang.2.files, 
       input_format   = 'gzFASTQ', 
       type           = 'rna',
       output_file    = wang.output.files, 
       output_format  = 'BAM',
       PE_orientation = 'fr', 
       nthreads       = 20, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

# I checked the strandness
# the annotation data was downloaded from 
# https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/
# hg38_RefSeq.bed.gz
#infer_experiment.py -r /mnt/date/genomelib/annotation/hg38_RefSeq.bed \
#-i /home/zhenyisong/biodata/wanglab/wangdata/wangyin/170801resuls/rsubread/CM-12d-1_S5.bam

wangyin.genes  <- featureCounts(  wang.output.files, 
                                  useMetaFeatures = TRUE,
                                  countMultiMappingReads = FALSE,
                                  strandSpecific         = 2, 
                                  isPairedEnd            = TRUE,
                                  requireBothEndsMapped  = TRUE,
                                  autosort               = TRUE,
                                  nthreads               = 20,
                                  annot.inbuilt          = 'hg38', 
                                  GTF.featureType        = 'exon',
                                  GTF.attrType           = 'gene_id',
                                  allowMultiOverlap      = TRUE)

setwd(Rdata.output.dir)
save.image('wangyin170801.Rdata')
quit('no')


