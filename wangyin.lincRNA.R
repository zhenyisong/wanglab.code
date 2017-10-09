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
                   'edgeR', 'limma', 'DESeq2', 'fastqcr','bamsignals',
                   'Rsamtools', 'org.Hs.eg.db', 'openxlsx', 'gridExtra',
                   'grid', 'rtracklayer', 'huex10sttranscriptcluster.db',
                   'oligo','annotate', 'pd.huex.1.0.st.v2',
                   'BSgenome.Hsapiens.UCSC.hg38')
load.lib     <- lapply(pkgs, require, character.only = TRUE)


# image file 
Rdata.output.dir   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/Rdata')
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
    set_node_attrs('label',c( 'H1.d1',  'Mesoderm.d3', 
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


"
setwd(Rdata.output.dir)
load('wangyin170801.Rdata')

"
x1.running.path <- file.path('D:\\wangli_data\\Rdata')
setwd(x1.running.path)
load('wangyin170801.Rdata')


#---
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

# this is from NONCODE database
#---
lncRmRNA.hg38.GTF.file <- file.path('/mnt/date/genomelib/annotation/hg38_UCSC_lncRNA_mRNA.gtf')
lncRmRNA.mm10.GTF.file <- file.path('/mnt/date/genomelib/annotation/mm10_lnc.protein.all.gtf')

# this is from gencode database
#  md5sum  gencode.vM15.annotation.gtf.gz
# 97a8e76a1e02f96b421bda741efb4a13  gencode.vM15.annotation.gtf.gz
#---
lncRmRNA.GRCh38.GTF.file  <- file.path('/mnt/date/genomelib/annotation/gencode/GRCh38.human.GTF')
lncRmRNA.GRCm38.GTF.file  <- file.path('/mnt/date/genomelib/annotation/gencode/GRCm38.mouse.GTF')


rawdata.path    <- '/mnt/date/Sequencing/FastQ/170801wangyinLincR'
wang.output.dir <- '/home/zhenyisong/biodata/wanglab/wangdata/wangyin/170801resuls/rsubread'
#---
# the procedure of fastqcR
#---
wang.fastQC.dir     <- '/home/zhenyisong/biodata/wanglab/wangdata/wangyin/170801resuls/fastqcR'
dir.create(file.path(wang.fastQC.dir), showWarnings = FALSE)
fastqc(fq.dir = rawdata.path, qc.dir = wang.fastQC.dir, threads = 8)
fastqc.results <- qc_aggregate(wang.fastQC.dir)



"

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

wangyin.lncR.mRNA.noncode  <- featureCounts(  wang.output.files, 
                                      useMetaFeatures        = TRUE,
                                      countMultiMappingReads = FALSE,
                                      strandSpecific         = 2, 
                                      isPairedEnd            = TRUE,
                                      requireBothEndsMapped  = TRUE,
                                      autosort               = TRUE,
                                      nthreads               = 20,
                                      annot.ext              = lncRmRNA.hg38.GTF.file,
                                      isGTFAnnotationFile    = TRUE, 
                                      GTF.featureType        = 'exon',
                                      allowMultiOverlap      = TRUE)


wangyin.lncR.mRNA.gencode  <- featureCounts(  wang.output.files, 
                                      useMetaFeatures        = TRUE,
                                      countMultiMappingReads = FALSE,
                                      strandSpecific         = 1, 
                                      isPairedEnd            = TRUE,
                                      requireBothEndsMapped  = TRUE,
                                      autosort               = TRUE,
                                      nthreads               = 20,
                                      annot.ext              = lncRmRNA.GRCh38.GTF.file,
                                      isGTFAnnotationFile    = TRUE, 
                                      GTF.featureType        = 'exon',
                                      allowMultiOverlap      = TRUE)
"
#---
# postQC, cross-check with public dataset
# GSE52313
# mouse
# PE,
# strandness?
# Ounzain S, Micheletti R, Beckmann T, Schroen B et al. 
# Genome-wide profiling of the cardiac transcriptome after myocardial 
# infarction identifies novel heart-specific long non-coding RNAs. 
# Eur Heart J 2015 Feb 7;36(6):353-68a. PMID: 24786300
#---


public.rawdata.path <- '/home/zhenyisong/biodata/wanglab/wangdata/wangyin/publicdata'
public.output.dir   <- '/home/zhenyisong/biodata/wanglab/wangdata/wangyin/publicdata/rsubread'


dir.create(file.path(public.output.dir), showWarnings = FALSE)
setwd(public.rawdata.path)
pedrazzini.lincR.files  <- list.files( path = public.rawdata.path, pattern = '.fastq$', 
                                         all.files = FALSE, full.names  = TRUE, 
                                         recursive = FALSE, ignore.case = FALSE, include.dirs = F) %>%
                            {.[1:16]}
pedrazzini.1.files           <- pedrazzini.lincR.files[grep('_1.fastq',pedrazzini.lincR.files)]
pedrazzini.2.files           <- pedrazzini.lincR.files[grep('_2.fastq',pedrazzini.lincR.files)]


pedrazzini.output.files     <- basename(pedrazzini.1.files) %>% 
                               sub(pattern = '_1.fastq', replacement = '') %>%
                               paste0(public.output.dir, '/', . ,'.bam')

setwd(rsubread.index.lib)
base.string           <- 'mm10'
align( index          = base.string, 
       readfile1      = pedrazzini.1.files , 
       readfile2      = pedrazzini.2.files, 
       input_format   = 'FASTQ', 
       type           = 'rna',
       output_file    = pedrazzini.output.files, 
       output_format  = 'BAM',
       PE_orientation = 'fr', 
       nthreads       = 8, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

#infer_experiment.py -r /mnt/date/genomelib/annotation/mm10_RefSeq.bed \
#-i /home/zhenyisong/biodata/wanglab/wangdata/wangyin/publicdata/rsubread/SRR1028904.bam

# the above strandness checking suggested that pedrazzini data was not stranded
#--- 
pedrazzini.genes <- featureCounts(  pedrazzini.output.files, 
                                    useMetaFeatures        = TRUE,
                                    countMultiMappingReads = FALSE,
                                    strandSpecific         = 0, 
                                    isPairedEnd            = TRUE,
                                    requireBothEndsMapped  = TRUE,
                                    autosort               = TRUE,
                                    nthreads               = 20,
                                    annot.ext              = lncRmRNA.mm10.GTF.file,
                                    isGTFAnnotationFile    = TRUE, 
                                    GTF.featureType        = 'exon',
                                    allowMultiOverlap      = TRUE)

pedrazzini.genes.gencode <- featureCounts(  pedrazzini.output.files, 
                                            useMetaFeatures        = TRUE,
                                            countMultiMappingReads = FALSE,
                                            strandSpecific         = 0, 
                                            isPairedEnd            = TRUE,
                                            requireBothEndsMapped  = TRUE,
                                            autosort               = TRUE,
                                            nthreads               = 20,
                                            annot.ext              = lncRmRNA.GRCm38.GTF.file,
                                            isGTFAnnotationFile    = TRUE, 
                                            GTF.featureType        = 'exon',
                                            allowMultiOverlap      = TRUE)



# 
# rsubread 
# https://groups.google.com/forum/#!topic/subread/JPPw9lVfgpw
# stat explaination
# I will do extra QC for wangyin data
#---



setwd(public.rawdata.path)
encode.lincR.files  <- list.files( path = public.rawdata.path, pattern = '.fastq$', 
                                         all.files = FALSE, full.names  = TRUE, 
                                         recursive = FALSE, ignore.case = FALSE, include.dirs = F) %>%
                       {.[17:30]}
encode.1.files           <- encode.lincR.files[grep('_1.fastq',encode.lincR.files)]
encode.2.files           <- encode.lincR.files[grep('_2.fastq',encode.lincR.files)]


encode.output.files     <- basename(encode.1.files) %>% 
                           sub(pattern = '_1.fastq', replacement = '') %>%
                           paste0(public.output.dir, '/', . ,'.bam')

setwd(rsubread.index.lib)
base.string           <- 'hg38'
align( index          = base.string, 
       readfile1      = encode.1.files , 
       readfile2      = encode.2.files, 
       input_format   = 'FASTQ', 
       type           = 'rna',
       output_file    = encode.output.files, 
       output_format  = 'BAM',
       PE_orientation = 'fr', 
       nthreads       = 8, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

# chece strandness
#infer_experiment.py -r /mnt/date/genomelib/annotation/hg38_RefSeq.bed \
#-i /home/zhenyisong/biodata/wanglab/wangdata/wangyin/publicdata/rsubread/SRR4421966.bam
#---


#---
# Encode data source
# GSE35585 SRX1603416 
# GSE35585 SRX1603417 
# GSE35585 SRX2243992
# GSE35585 SRX2244046
# GSE35585 SRX2244088
# GSE35585 SRX2244108
# GSE35585 SRX2244255
#---





encode.genes <- featureCounts(  encode.output.files, 
                                useMetaFeatures        = TRUE,
                                countMultiMappingReads = FALSE,
                                strandSpecific         = 2, 
                                isPairedEnd            = TRUE,
                                requireBothEndsMapped  = TRUE,
                                autosort               = TRUE,
                                nthreads               = 20,
                                annot.ext              = lncRmRNA.hg38.GTF.file,
                                isGTFAnnotationFile    = TRUE, 
                                GTF.featureType        = 'exon',
                                allowMultiOverlap      = TRUE)

encode.genes.gencode <- featureCounts(  encode.output.files, 
                                        useMetaFeatures        = TRUE,
                                        countMultiMappingReads = FALSE,
                                        strandSpecific         = 2, 
                                        isPairedEnd            = TRUE,
                                        requireBothEndsMapped  = TRUE,
                                        autosort               = TRUE,
                                        nthreads               = 20,
                                        annot.ext              = lncRmRNA.GRCh38.GTF.file,
                                        isGTFAnnotationFile    = TRUE, 
                                        GTF.featureType        = 'exon',
                                        allowMultiOverlap      = TRUE)

# obtain all processed raw data for x1 analysis
#---
gencode.GRCh38.human    <- import(lncRmRNA.GRCh38.GTF.file)twd(Rdata.output.dir)


# to check the data pattern in wang lncRNA design
# using the washingU ,murry data set the golden standard?
# esc,  mesoderm( cardiac mesoderm?), cardiac progenitor stage  , defferentiation stages?
# wang yin data lack this stage, I am not sure why he design such such way.
#---

# to compare and validate the wangyin data with
# murry data set which was published in CELL
# data source:
# GSE19090 
# [HuEx-1_0-st] Affymetrix Human Exon 1.0 ST Array
#
#---


GSE19090.rawdata   <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/wangyin/publicdata/GSE19090')
setwd(GSE19090.rawdata)
GSE19090.rma       <- list.files( getwd(), pattern = '\\.CEL|\\.CEL\\.gz', 
                                  full.names = TRUE, ignore.case = TRUE) %>% 
                      read.celfiles(celPath = ., pkgname = 'pd.huex.1.0.st.v2') %>%
                      oligo::rma(target = 'core')
GSE19090.symbols          <- featureNames(GSE19090.rma) %>% 
                             {mapIds( huex10sttranscriptcluster.db, 
                                      column = 'SYMBOL', keys = ., keytype = 'PROBEID',
                                      multiVals = 'first')} %>%
                             make.names(unique = T)
GSE19090.exprs            <- exprs(GSE19090.rma)
colnames(GSE19090.exprs)  <- c( 'H7.T5.1','H7.T5.2','H7.T9.1','H7.T9.2','H7.T14.1',
                                'H7.T14.2', 'H7.T2.1','H7.T2.2','H7.T0.1','H7.T0.2')
rownames(GSE19090.exprs)  <- GSE19090.symbols


# GSE58363
# PMID: 25296024
# Early patterning and specification of cardiac 
# progenitors in gastrulating mesoderm. Elife 2014 
# mouse
# publication used the NONCODE annotation
# PE
# ---

setwd(Rdata.output.dir)
save.image('wangyin170801.Rdata')
quit('no')


#---
# data from Encode
# human heart tissue
# to check the QC and seqeunceing depth
#---


"
bf <- BamFile(wang.output.files[1])
sorted.bam.files <- sortBam(bf,  destination = wang.output.dir)
indexed.bam.files <- indexBam(sorted.bam.files, destination = wang.output.dir)
countBam(sorted.bam.files)
"

wangyin.stat <- wangyin.lncR.mRNA.gencode$stat
colnames(wangyin.stat) <- c( 'status','mesoderm.d3.1','mesoderm.d3.2','CM.d12.1',
                             'CM.d12.2', 'CM.d5.1', 'CM.d5.2', 'CM.d9.1','CMd9.2',
                             'EC.d12.1', 'EC.d12.2', 'EC.d6.1', 'EC.d6.2', 'EC.d9.1',
                             'EC.d9.2',  'H1.1','H1.2','VSMC.d12.1', 'VSMC.d12.2',
                             'VSMC.d6.1', 'VSMC.d6.2', 'VSMC.d9.1', 'VSMC.d9.2' )

encode.stat           <- encode.genes.gencode$stat
colnames(encode.stat) <- c( 'status', 'SRR3192433', 'SRR3192434', 'SRR4421689',
                            'SRR4421747', 'SRR4421789', 'SRR4421817', 'SRR4421966')

pedrazzini.stat      <- pedrazzini.genes.gencode$stat
colnames(pedrazzini.stat) <- c( 'status', 'SRR1028904', 'SRR1028905', 'SRR1028906',
                                'SRR1028907', 'SRR1028908', 'SRR1028909', 
                                'SRR1028910', 'SRR1028911')




"
setwd('D:\\yisong.data\\')
load('lncRNA.LR.Rdata')
wang.ucsf.stat           <- lncAll.LR.genes.gencode$stat
colnames(wang.ucsf.stat) <- c( 'status', 'SRR4044044', 'SRR4044045', 'SRR4044046',
                               'SRR4044047', 'SRR4044048', 'SRR4044049', 
                               'SRR4044050', 'SRR4044051', 'SRR4044052','SRR4044053',
                               'SRR4044054', 'SRR4044055', 'SRR4044056', 'SRR4044057',
                               'SRR4044058', 'SRR4044059', 'SRR4044060')

table.seq.depth <- grid.newpage() %>% 
                   {tableGrob( wangyin.stat )} %>%
                   grid.draw()
"

#---
# post QC fro wangyin data using gencode annotation
#---
rnaseq.rlog.func   <- function(genes,annot) {
                          genes %$% counts %>% 
                          DGEList(genes = annot) %>% 
                          calcNormFactors() %>% rpkm(log = T)
                      }
wangyin.lincR.rpkm <- rnaseq.rlog.func( wangyin.lncR.mRNA.gencode, 
                                        wangyin.lncR.mRNA.gencode$annotation)

wangyin.colnames   <- c( 'mesoderm.d3.1','mesoderm.d3.2','CM.d12.1',
                         'CM.d12.2', 'CM.d5.1', 'CM.d5.2', 'CM.d9.1','CM.d9.2',
                         'EC.d12.1', 'EC.d12.2', 'EC.d6.1', 'EC.d6.2', 'EC.d9.1',
                         'EC.d9.2',  'H1.1','H1.2','VSMC.d12.1', 'VSMC.d12.2',
                         'VSMC.d6.1', 'VSMC.d6.2', 'VSMC.d9.1', 'VSMC.d9.2' )
colnames(wangyin.lincR.rpkm) <- wangyin.colnames
apply.mean.func    <- function(df.matrix, col_nums) {
                          apply(df.matrix[,col_nums],1, mean)
                      }
comparision.list   <- list( H1.d1 = 15:16, Mesoderm.d3 = 1:2,
                            CM.d5 = 5:6, CM.d9 = 7:8, CM.d12 = 3:4,
                            EC.d5 = 11:12, EC.d9 = 13:14, EC.d12 = 9:10,
                            VSMC.d5 = 19:20, VSMC.d9 = 21:22, VSMC.d12 = 17:18 )
GRCh38.genenames <- wangyin.lncR.mRNA.gencode %$% annotation %$% 
                    GeneID %>% sub('.\\d+$','',., perl = F) %>%
                    mapIds( org.Hs.eg.db, keys = ., column = 'SYMBOL', 
                            keytype = 'ENSEMBL', multiVals = 'first') %>%
                    make.names( unique = TRUE)

GRCh38.extra.annot      <- gencode.GRCh38.human %$% 
                           {data.frame( ensembl.id = gene_id, 
                                          gene       = type, 
                                          gene_type  = gene_type)}%>%
                           filter(gene == 'gene') %>% unique()

wangyin.extra.annot     <- wangyin.lncR.mRNA.gencode %$% 
                           annotation %$% 
                           GeneID %>% 
                           data.frame(GeneID = .) %>%
                           inner_join( GRCh38.extra.annot, 
                                       by = c('GeneID' = 'ensembl.id'))
time.point.rpkm  <- map(comparision.list, apply.mean.func, df.matrix = wangyin.lincR.rpkm) %>%
                    as.data.frame() %>% mutate(GeneName = GRCh38.genenames) %>%
                    mutate(EnsemblID = wangyin.lncR.mRNA.gencode$annotation$GeneID) %>%
                    cbind(wangyin.extra.annot)


# PCA analysis for wangyin, data replication check
#---

cell.groups <- wangyin.colnames %>%
               sub('.\\d$', '', perl = F, .) %>% as.factor() %>%
               data.frame(cellGroups = .)
wangyin.vsd.data <- wangyin.lncR.mRNA.gencode %$% counts %>% 
                    DESeqDataSetFromMatrix( countData = ., 
                                           colData = cell.groups, 
                                           design = ~ cellGroups) %>%
                    DESeq() %>% getVarianceStabilizedData()
colnames(wangyin.vsd.data) <- wangyin.colnames

wangyin.pca <- wangyin.vsd.data %>% t() %>% prcomp()

# quick and dirty map
# deprecated
#---
"
plot(wangyin.pca$x, col = 'white', main = 'wangyin LncRNA project')
text(wangyin.pca$x[,1], wangyin.pca$x[,2], labels = colnames(wangyin.vsd.data), cex = 0.7)
"

# please see the original code
# from cardioMaturation.R
#---

pca.no         <- dim(wangyin.pca$x)[2]

scree.plot.ggplot <- {wangyin.pca$sdev^2} %>%
                     {./sum(.)} %>% 
                     {data.frame(variance = ., pca = c(1:pca.no)) } %>%
                     ggplot(data = .) +
                     xlab('Wangyin lincRNA Principle Components') +
                     ylab('Proportion of Variance Explained') +
                     scale_x_continuous( breaks = c(1:pca.no), 
                                         labels = as.character(c(1:pca.no), 
                                         limits = as.character(c(1:pca.no)))) +
                     geom_point(aes(x = pca, y = variance), size = 3) +
                     geom_line(aes(x = pca, y = variance), size = 0.8) +
                     scale_linetype_discrete() +
                     theme(legend.position = 'none') +
                     theme_classic()
scree.plot.ggplot

cumsum.plot.ggplot <- {wangyin.pca$sdev^2} %>%
                      {./sum(.)} %>% 
                      {data.frame(variance = cumsum(.), pca = c(1:pca.no))}%>%
                      ggplot(data = .) +
                      xlab('Wangyin lncR Principle Components') +
                      ylab('Culmulative Proportion of Variance Explained') +
                      scale_x_continuous( breaks = c(1:pca.no), labels = as.character(c(1:pca.no), 
                                          limits = as.character(c(1:pca.no)))) +
                      geom_point(aes(x = pca, y = variance), size = 3) +
                      geom_line(aes(x = pca, y = variance),size = 0.8) +
                      scale_linetype_discrete() +
                      theme(legend.position = 'none') +
                      theme_classic()
cumsum.plot.ggplot
                       


wangyin.final.ggplot <- as.data.frame(wangyin.pca$x) %>% 
                        mutate(color.choice = { rownames(.) %>% 
                                                sub('.\\d$', '', perl = F, .) %>%
                                                as.factor}) %>%
                        mutate(color.choice = factor( color.choice, 
                                                      levels(color.choice)[c(7,8,2,3,1,5,6,4,10,11,9)]) )%>%
                        ggplot(data = . ) + 
                        geom_point(aes(x = PC1, y = PC2, color = color.choice), size = 3) + 
                        scale_colour_manual( name   = 'lncRNA sample stages',
                                             values = c( 'gray8','blueviolet', 'brown1','brown2',
                                                         'brown3','green1','green2',
                                                         'green4','darkgoldenrod1',
                                                         'darkgoldenrod2','darkgoldenrod3'),
                                             labels = c( 'H1','mesoderm.d3', 'CM.5', 'CM.9','CM.12',
                                                         'EC.6','EC.9','EC.12','VSMC.6',
                                                         'VSMC.9','VSMC.12')) +
                        theme_classic()

source.pca.data <-  as.data.frame(wangyin.pca$x) %>% 
                    mutate(color.choice = { rownames(.) %>% 
                                            sub('.\\d$', '', perl = F, .) %>%
                                            as.factor}) %>%
                    mutate(color.choice = factor( color.choice, 
                                                  levels(color.choice)[c(7,8,2,3,1,5,6,4,10,11,9)]) ) %>%
                    dplyr::select(PC1,PC2, color.choice) %>% rename(developmental.stages = color.choice)

source.pca.table <- grid.newpage() %>% 
                   {tableGrob( source.pca.data)} %>%
                   grid.draw()
wangyin.final.ggplot
           

#time.point.rpkm

# who is the driver to turn h1 to mesoderm, mRNA?
# who is the guidian to turn the mesoderm to three lineage
# H1, cardiac mesoderm, other three lineage ChIP-seq data
# where to donwload?
# gene.QC.check <- function(genesymbol = gene) 
#---

markers.list    <- list( mesoderm   = c('T','EOMES','MESP1','MESP2','KDR','PDGFRA'), 
                         endothlium = c('PECAM1','FLT1','VWF','TEK','CDH5','KDR','ENG'),
                         VSMC       = c('MYH11','TAGLN','ACTA2','CNN1'),
                         myocardium = c( 'TNNT2','TNNI3','GATA4','NKX2-5','MYH7','MYH6','ACTN1',
                                         'GJA1', 'MYL2', 'GJA5','GJC1','MEF2C','FENDRR', 'TBX20',
                                         'HAND1','HAND2','PITX2','MEF2C','ACTN2'),
                         hESC       = c('B3GALT5','POU5F1'),
                         A.process  = c('MSX2','PCGF2','KDM6B'))

markers.df      <- time.point.rpkm %>% 
                   filter(GeneName %in% unlist(markers.list))


QC.genename         <- 'KDM6B'
gene.pattern.ggplot <- time.point.rpkm %>% 
                       dplyr::select(H1.d1:VSMC.d12,GeneName) %>% 
                       filter(GeneName ==  QC.genename) %>%
                       dplyr::select(H1.d1:VSMC.d12) %>% 
                       dplyr::select(c(1:5)) %>% t() %>%
                       as.data.frame() %>%
                       ggplot( data = ., 
                               aes( x = as.factor(rownames(.)), 
                                    y = V1, group = 1)) +
                        geom_line( linetype = 'dashed', 
                                   color = 'steelblue') + 
                        geom_point(color = 'steelblue') +
                        scale_x_discrete( limits = c('H1.d1','Mesoderm.d3',
                                                     'CM.d5', 'CM.d9', 'CM.d12')) +
                        xlab('developmental stages') +
                        ylab('gene expression log(rpkm)') +
                        ggtitle( paste('gene symbol', QC.genename, 'QC checke', sep = ' '))
             

gene.pattern.ggplot




               


setwd('C:\\Users\\Yisong\\Desktop')
wangyin.wb <- createWorkbook('wangyin')
addWorksheet(wangyin.wb, 'log_whole_data_rpkm')
addWorksheet(wangyin.wb, 'log_with_replication')
addWorksheet(wangyin.wb, 'selected_markers_in_stages')
addWorksheet(wangyin.wb, 'sequence_depth')
addWorksheet(wangyin.wb, 'sequence_depth_encode')
addWorksheet(wangyin.wb, 'sequence_depth_pedrazzini')
addWorksheet(wangyin.wb, 'sequence_depth_wangUCSF')

writeData(wangyin.wb, sheet = 1, time.point.rpkm )
writeData(wangyin.wb, sheet = 2, wangyin.lincR.rpkm)
writeData(wangyin.wb, sheet = 3, markers.df )
writeData(wangyin.wb, sheet = 4, wangyin.stat )
writeData(wangyin.wb, sheet = 5, encode.stat )
writeData(wangyin.wb, sheet = 6, pedrazzini.stat )
writeData(wangyin.wb, sheet = 7, wang.ucsf.stat )

saveWorkbook(wangyin.wb, 'wangyin-2017-09-28.xlsx', overwrite = TRUE)


# now parse the GTF informantion and 
# determine which is the lncRNA gene and its 
# gene ensembl ID
#---



"
setwd(Rdata.output.dir)
save.image('wangyin170801.Rdata')
quit('no')
"
