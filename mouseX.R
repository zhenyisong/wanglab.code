# the script is to demultiplexing the result from next-500
# @rawdate path
# @samplesheet edited by Yisong, log record from Zongnna, lib_id, index
#---
#baseCalls='/mnt/date/Sequencing/RawData/0314/170314_NB501912_0005_AH3WNYBGX2/Data/Intensities/BaseCalls'
#sampleSheet='/mnt/date/Sequencing/FastQ/sampleSheet/2017_03_14_MouseX.csv'
#dataOutput='/mnt/date/Sequencing/FastQ/mouseXRNA_2017_03_14'
#runFolder='/mnt/date/Sequencing/RawData/0314/170314_NB501912_0005_AH3WNYBGX2'
#bcl2fastq --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions \
#          --no-lane-splitting -R $runFolder -i $baseCalls -r 10 -d 10 -p 10 \
#          --sample-sheet $sampleSheet -o $dataOutput 2>> /dev/null &
#cd /mnt/date/Sequencing/FastQ/mouseXRNA_2017_03_14
#mkdir mouseX_QC
#files=(sample*.gz)
#for filename in ${files[@]};do
#    fastqc -t 15 -q -o mouseX_QC $filename
#done
#multiqc mouseX_QC 

library(gplots)
library(xlsx)
library(Rsubread)
library(edgeR)
library(limma)
library(DESeq2)
library(gplots)
library(genefilter)
library(RColorBrewer)
library(org.Mm.eg.db)
library(cluster)
library(factoextra)
library(clusterProfiler)
library(pathview)
library(sva)
library(systemPipeR)
library(rtracklayer)
library(stringr)
library(GenomicFeatures)

setwd('/home/zhenyisong/biodata/wanglab/wangdata/mouseX/rsubread')
load('mouseX.Rdata')
mouse.genome_ref.path   <- "/mnt/date/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
setwd('/mnt/date/Sequencing/FastQ/mouseXRNA_2017_03_14')
reads.files.names       <- list.files(pattern = "^sample.*\\.fastq\\.gz$")
raw.data.path           <- getwd()
setwd('/home/zhenyisong/biodata/wanglab/wangdata/mouseX')
#unlink('rsubread', force = TRUE, recursive = TRUE)
#dir.create('rsubread')
output.path             <- '/home/zhenyisong/biodata/wanglab/wangdata/mouseX/rsubread'

setwd('/home/zhenyisong/biodata/wanglab/wangdata/mouseX/rsubread')
mouse.base              <- 'mm10_index'

read.path.1             <- reads.files.names[grep("R1",reads.files.names )]
mouseX.outputs.files    <- paste0(output.path,'/', read.path.1,'.bam')
read.path.1             <- paste0(raw.data.path, '/',read.path.1)
read.path.2             <- reads.files.names[grep("R2",reads.files.names )]
read.path.2             <- paste0(raw.data.path, '/',read.path.2)

#buildindex( basename = mouse.base, reference = mouse.genome_ref.path )
#
#align( index          = mouse.base, 
#       readfile1      = read.path.1, 
#       readfile2      = read.path.2, 
#       input_format   = "gzFASTQ", 
#       type           = 'rna',
#       output_file    = mouseX.outputs.files, 
#       output_format  = "BAM",
#       PE_orientation = 'fr', 
#       nthreads       = 15, 
#       indels         = 1,
#       maxMismatches  = 3,
#       phredOffset    = 33,
#       unique         = T )

rsubreadQC <- align( index          = mouse.base, 
                     readfile1      = read.path.1, 
                     input_format   = "gzFASTQ", 
                     type           = 'rna',
                     output_file    = mouseX.outputs.files, 
                     output_format  = "BAM",
                     PE_orientation = 'fr', 
                     nthreads       = 15, 
                     indels         = 1,
                     maxMismatches  = 3,
                     phredOffset    = 33,
                     unique         = T )

mm10.genes       <- featureCounts( mouseX.outputs.files, useMetaFeatures = TRUE, 
                                   annot.inbuilt = "mm10", allowMultiOverlap = TRUE,
                                   nthreads = 15, strandSpecific = 0)

gene.counts <- mm10.genes$counts
gene.ids    <- mm10.genes$annotation$GeneID

columns     <- c("ENTREZID","SYMBOL", "MGI", "GENENAME");
GeneInfo    <- select( org.Mm.eg.db, keys= as.character(gene.ids), 
                       keytype="ENTREZID", columns = columns);
m           <- match(mm10.genes$annotation$GeneID, GeneInfo$ENTREZID);
Ann         <- cbind( mm10.genes$annotation[, c("GeneID", "Chr","Length")],
                      GeneInfo[m, c("SYMBOL", "MGI", "GENENAME")]);
Ann$Chr     <-  unlist( lapply(strsplit(Ann$Chr, ";"), 
                    function(x) paste(unique(x), collapse = "|")))
Ann$Chr     <- gsub("chr", "", Ann$Chr)
gene.counts.export <- cbind(Ann, gene.counts)
write.csv(gene.counts.export, file = 'mouseX.csv')
gene.exprs  <- DGEList(counts = gene.counts, genes = Ann)


gene.exprs  <- calcNormFactors(gene.exprs)
dge.tmm         <- t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
dge.tmm.counts  <- apply(dge.tmm,2, as.integer)

sample.info              <- data.frame( treat  = c('sample_11','sample_12','sample_13',
                                                   'sample_14','sample_15','sample_16') )
dds                      <- DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                                    colData   = sample.info,
                                                    design    = ~ treat)
vsd                      <- varianceStabilizingTransformation(dds)
vsd.exprs                <- assay(vsd)
colnames(vsd.exprs)      <- c('sample_11','sample_12','sample_13',
                              'sample_14','sample_15','sample_16')
mouseX.PCA  <- prcomp(t(vsd.exprs))
pr.var      <- mouseX.PCA$sdev^2
pve         <- pr.var/sum(pr.var)

pve.df <- data.frame(variance = pve, pca = c(1:6))

pve.pdf <- ggplot(pve.df) +
           xlab('Principle Component') +
           ylab('Proportion of Variance Explained') +
           scale_x_continuous( breaks = c(1:6), labels = as.character(c(1:6), 
                               limits = as.character(c(1:6)))) +
           geom_point(aes(x = pca, y = variance), size = 3) +
           geom_line(aes(x = pca, y = variance), size = 0.8) +
           scale_linetype_discrete() +
           theme(legend.position = "none")

mosueX.cord <- as.data.frame(mouseX.PCA$x)

mosueX.cord$cardio.type<- factor(1:6)
ggplot(mosueX.cord) + 
           geom_point(aes(x = PC1, y = PC2, color = cardio.type), size = 3) + 
           scale_colour_manual( name   = 'sample classification',
                                values = c("cadetblue", "cadetblue3", "deeppink","royalblue",
                                           "deeppink1","royalblue1"),
                                labels = c( 'sample_11','sample_12','sample_13',
                                            'sample_14','sample_15','sample_16')) +
           theme( legend.position    = 'bottom',
                  legend.direction   = 'horizontal',
                  legend.title.align = 0.5, 
                  legend.text        = element_text(size = 7)) +
           guides( color = guide_legend(title.position = 'top') )

save.image(file = 'mouseX.Rdata')
quit('no')