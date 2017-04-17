#
#
#  fastqc -t 8 -f fastq -o lncRNA_QC /mnt/date/Sequencing/FastQ/totalRNA_2017_03_14/A*.fastq


library(gplots)
library(xlsx)
library(Rsubread)
library(edgeR)
library(limma)
library(org.Hs.eg.db)
library(DESeq2)
library(gplots)
library(genefilter)
library(RColorBrewer)
library(org.Rn.eg.db)
library(cluster)
library(factoextra)
library(clusterProfiler)
library(pathview)
library(sva)
library(systemPipeR)
library(rtracklayer)
library(stringr)
library(GenomicFeatures)

setwd('/home/zhenyisong/biodata/wanglab/wangdata/lncRNAClinic/rsubread')
setwd('D:\\wangli_data\\Rdata')
load('lincRNA_SE.Rdata')

human.genome_ref.path   <- "/mnt/date/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
human.hg19.lncGene.GTF  <- '/mnt/date/genomelib/annotation/human_hg19_GenePusLncRNA.gtf'
human.hg19.lncRNA.GTF   <- '/mnt/date/genomelib/annotation/NONCODE2016.hg19.gtf'
human.hg38.lncRNA.fpkm  <- '/mnt/date/biodata/wanglab/wangdata/lncRNAClinic/annotation/humanLncGene.fpkm'
setwd('/mnt/date/Sequencing/FastQ/totalRNA_2017_03_14')
reads.files.names       <- list.files(pattern = '*.fastq$')[1:16]
raw.data.path           <- '/mnt/date/Sequencing/FastQ/totalRNA_2017_03_14'
setwd('/home/zhenyisong/biodata/wanglab/wangdata/lncRNAClinic')
#unlink('rsubread', force = TRUE, recursive = TRUE)
#dir.create('rsubread')
output.path             <- '/home/zhenyisong/biodata/wanglab/wangdata/lncRNAClinic/rsubread'
setwd('/home/zhenyisong/biodata/wanglab/wangdata/lncRNAClinic/rsubread')

human.base              <- 'h19_index'

read.path.1             <- reads.files.names[grep("R1",reads.files.names )]
lincRNA.outputs.files   <- paste0(output.path,'/', read.path.1,'.bam')
read.path.1             <- paste0(raw.data.path, '/',read.path.1)
read.path.2             <- reads.files.names[grep("R2",reads.files.names )]
read.path.2             <- paste0(raw.data.path, '/',read.path.2)

# generate QC report
name.vector.args        <- c(read.path.1, read.path.2)
sample.names            <- sub('_001.fastq','',basename(name.vector.args))
names(name.vector.args) <- sample.names

#rawdata.list           <- seeFastq(fastq = name.vector.args, batchsize = 10000, klength = 8)
#pdf("fastqReport.pdf", height = 18, width = 4 * length(rawdata.list))
#seeFastqPlot(rawdata.list )
#dev.off()
#
#buildindex( basename = human.base, reference = human.genome_ref.path )
#align( index          = human.base, 
#       readfile1      = read.path.1, 
#       readfile2      = read.path.2, 
#       input_format   = "FASTQ", 
#       type           = 'rna',
#       output_file    = lincRNA.outputs.files, 
#       output_format  = "BAM",
#       PE_orientation = 'fr', 
#       nthreads       = 8, 
#       indels         = 1,
#       maxMismatches  = 3,
#       phredOffset    = 33,
#       unique         = T )
#  

align( index         = human.base, 
       readfile1     = read.path.1, 
       input_format  = "FASTQ", 
       type          = 'rna',
       output_file   = lincRNA.outputs.files, 
       output_format = "BAM", 
       nthreads      = 15, 
       indels        = 1,
       maxMismatches = 3,
       phredOffset   = 33,
       unique        = T )


hg19.lncGene        <- featureCounts( lincRNA.outputs.files, useMetaFeatures = TRUE, 
                                      annot.ext  = human.hg19.lncGene.GTF, isGTFAnnotationFile = TRUE,
                                      nthreads   = 15, strandSpecific = 2, allowMultiOverlap = TRUE)

gene.counts        <- hg19.lncGene$counts
gene.exprs         <- DGEList(counts = gene.counts)
gene.exprs         <- calcNormFactors(gene.exprs)

colnames(gene.counts) <- c('70OM_A','70PM_A','7OM_A','7PM_A','70OM_B','70PM_B','7OM_B','7PM_B')
write.csv(gene.counts, file = 'lncRNA_PEmodel.csv')
dge.tmm               <- t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
log.counts            <- log(dge.tmm + 1.0)
colnames(log.counts)  <- c('70OM_A','70PM_A','7OM_A','7PM_A','70OM_B','70PM_B','7OM_B','7PM_B')
write.csv(log.counts, file = 'lncRNA_PEmodel.logTrans.csv')

log.genes             <- cbind( log.counts[,1] - log.counts[,3], log.counts[,1]- log.counts[,2],
                                log.counts[,3] - log.counts[,4], log.counts[,2] - log.counts[,4],
                                log.counts[,5] - log.counts[,7], log.counts[,5] - log.counts[,6],
                                log.counts[,7] - log.counts[,8], log.counts[,6] - log.counts[,8])

colnames(log.genes) <- c('70OM-TOM/A','70OM-70PM/A','TOM-TPM/A','70PM-TPM/A','70OM-TOM/B','70OM-70PM/B','TOM-TPM/B','70PM-TPM/B')

write.csv(log.genes, file = 'lncRNA_PEmodel.LFC.csv')


# post-processing QC
#---

##hg19.lnRNA.GTF.transcript <- import.gff( human.hg19.lncRNA.GTF, format = "gtf", genome = "hg19",
##                              feature.type = "transcript", sequenceRegionsAsSeqinfo = TRUE)
##hg19.lncRNA.gene.db       <-  makeTxDbFromGFF(human.hg19.lncRNA.GTF, format = 'gtf')
##gene.len      <- width( ranges(hg19.lnRNA.GTF) )

lncRNA.index         <- grepl( 'NONHSAG', rownames(gene.counts), ignore.case = FALSE)
hg19.lncRNA.names    <- sub("\\.\\d+$","", rownames(gene.counts)[lncRNA.index])
hg19.lncRNA.len      <- hg19.lncGene$annotation$Length[lncRNA.index]
hg19.lncRNA.exprs    <- gene.counts[lncRNA.index,]
hg19.lncRNA.rpkm     <- rpkm(hg19.lncRNA.exprs, gene.length = hg19.lncRNA.len )
hg38.lncRNA.pub.fpkm <- read.table(human.hg38.lncRNA.fpkm, header = F)



# this is tedous long time 
# should use the do.call instead of 
#----

results                   <- hg38.lncRNA.pub.fpkm$V3
hg38.lncRNA.fpkm.matrix   <- NULL

for( i in (1:length(results))) {
   hg38.lncRNA.fpkm.matrix  <- rbind(hg38.lncRNA.fpkm.matrix , as.double(str_split(results[i],",")[[1]]) )
}

rownames(hg38.lncRNA.fpkm.matrix) <- hg38.lncRNA.pub.fpkm$V1
match.index                       <- match(hg19.lncRNA.names,hg38.lncRNA.pub.fpkm$V1)
hg19.lncRNA.allTissue.fpkm        <- hg38.lncRNA.fpkm.matrix[match.index,]
setwd('/home/zhenyisong/biodata/wanglab/wangdata/lncRNAClinic/annotation')
tissues.colnames     <- read.table(file = 'lncRNA_tissue.txt', header = F,  comment.char = "#")

samplePlusAll.fpkm   <- cbind(hg19.lncRNA.allTissue.fpkm, hg19.lncRNA.rpkm)
mean(is.na(hg19.lncRNA.allTissue.fpkm))
# 0.005390896
# http://stackoverflow.com/questions/3798998/cor-shows-only-na-or-1-for-correlations-why
cor(log(samplePlusAll.fpkm + 1), method = 'spearman', use = "complete.obs")



#setwd('/home/zhenyisong/biodata/wanglab/wangdata/lncRNAClinic/rsubread')
#save.image(file = 'lincRNA_SE.Rdata')
#quit('no')




#
# branch -- see the code at
setwd('E:\\FuWai\\wangli.lab\\lncRNAClinic')
hisat2.counts <- read.csv( file = 'transcript_count_matrix.csv',
                           header = TRUE, row.names = 1)
#group.A.counts <- hisat2.counts[,1:4] + 1
#group.B.counts <- hisat2.counts[,5:8] + 1
gene.exprs     <- DGEList(counts = hisat2.counts) 
gene.exprs     <- calcNormFactors(gene.exprs)

norm.counts    <- t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
norm.counts    <- norm.counts + 1

keytypes(org.Hs.eg.db)

gene.refseq   <- rownames(hisat2.counts)

columns  <- c("ENTREZID","SYMBOL", "REFSEQ", "GENENAME");
GeneInfo <- select( org.Hs.eg.db, keys = as.character(gene.refseq), 
                    keytype = "REFSEQ", columns = columns);

norm.counts.ann <- cbind(norm.counts,GeneInfo  )
norm.counts.log.ann <- cbind(log(norm.counts), GeneInfo)
setwd("C:\\Users\\Yisong\\Desktop")
write.csv( norm.counts.ann, file = 'lncRNAClinic.csv')
write.csv( norm.counts.ann, file = 'lncRNAClinic.log.txt')


# cannot use the limma to process the data and caculate the p.value
#---