# the raw data were generated from next-seq 500
# SE or PE?

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


# cd /home/zhenyisong/data/wanglilab/wangcode
# nohup R CMD BATCH multiple.xinli.R &
genome_ref.path      <- "/home/zhenyisong/data/bringback/igenome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
setwd('/home/zhenyisong/data/wanglilab/projects/xinli')
reads.files.names    <- list.files(pattern = '*.fastq')
reads.paths          <- paste0(getwd(),'/',reads.files.names)

setwd('/home/zhenyisong/data/wanglilab/projects/xinli')
unlink('rsubread')
dir.create('rsubread')
output.path          <- '/home/zhenyisong/data/wanglilab/projects/xinli/rsubread'
setwd('/home/zhenyisong/data/wanglilab/projects/xinli/rsubread')

outputs.files        <- paste0(output.path,'/', reads.files.names,'.bam')
base.string          <- 'mm10_index'

"
use the Rsubread command to generate index file
this index file will be generated and saved at getwd()
you do not need to generate the script
"
setwd("/home/zhenyisong/data/wanglilab/projects/xinli/rsubread")
buildindex( basename = base.string, reference = genome_ref.path )

"
this is the function which is called to align the genome
sequence
"

align( index         = base.string, 
       readfile1     = reads.paths, 
       input_format  = "FASTQ", 
       type          = 'rna',
       output_file   = outputs.files, 
       output_format = "BAM", 
       nthreads      = 8, 
       indels        = 1,
       maxMismatches = 3,
       phredOffset   = 33,
       unique        = T )

# get gene's counts
gene         <-  featureCounts( outputs.files, useMetaFeatures = TRUE, 
                                annot.inbuilt = "mm10", allowMultiOverlap = TRUE)

setwd('/home/zhenyisong/data/cardiodata')
save.image(file = 'multiple.xinli.Rdata')
quit("no")  


gene.counts  <- gene$counts

colnames(gene.counts)
gene.ids     <- gene$annotation$GeneID
colnames(gene.counts) <- c( 'CM_GFPOE_1','CM_Ino80_1_1','CM_Ino80_1_2',
                            'CM_Ino80_3_1','CM_Ino80_3_2','CM_NT_1',
                            'CM_NT_2','CM_Rif1B_1','CM_Rif1C_1',
                            'CM_Rif1C_2','CM_Rif1OE_1','CM_Rif1OE_2');

keytypes(org.Hs.eg.db)

columns  <- c("ENTREZID","SYMBOL", "OMIM", "GENENAME");
GeneInfo <- select( org.Hs.eg.db, keys= as.character(gene.ids), 
                   keytype="ENTREZID", columns = columns);
m        <- match(gene$annotation$GeneID, GeneInfo$ENTREZID);
Ann      <- cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                   GeneInfo[m, c("SYMBOL", "OMIM", "GENENAME")]);

Ann$Chr  <-  unlist( lapply(strsplit(Ann$Chr, ";"), 
                    function(x) paste(unique(x), collapse = "|")))
Ann$Chr  <- gsub("chr", "", Ann$Chr)

gene.exprs      <- DGEList(counts = gene.counts, genes = Ann)
gene.exprs      <- calcNormFactors(gene.exprs)
dge.tmm         <- t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
dge.tmm.counts  <- apply(dge.tmm,2, as.integer)

group           <- factor(c( 'GFP_oe','Ino80_1','Ino80_1','Ino80_3',
                             'Ino80_3','CM_NT','CM_NT','Rif1B',
                             'Rif1C','Rif1C','Rif1OE','Rif1OE'));
sample.info     <- data.frame(treat = group)
dds             <- DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                           colData   = sample.info,
                                           design    = ~ treat)
vsd                     <- varianceStabilizingTransformation(dds, blind = FALSE);
vsd.expr.wholeGroup                 <- assay(vsd)
rownames(vsd.expr.wholeGroup )      <- GeneInfo[m,'SYMBOL']
colnames(vsd.expr.wholeGroup )      <- colnames(gene.counts) 

setwd('/wa/zhenyisong/wanglilab/wangdata') 
mature.im.markers <- 'cardio_manual_maturation_markers.xlsx'
mature.im.table   <- read.xlsx(mature.im.markers, header = TRUE, stringsAsFactors = FALSE, sheetIndex = 1)  
mature.im.names   <- toupper( mature.im.table$geneSymbol )
mature.im.heatmap <- vsd.expr[mature.im.names,] 
my_palette        <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
heatmap.result    <- heatmap.2(mature.im.heatmap, col = my_palette, scale  = 'row', 
						       Rowv = TRUE,Colv = FALSE, density.info = 'none',
                               key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
						       cexRow  = 1, cexCol = 1,srtCol = 30,
                               distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						       hclustfun  = function(d) hclust(d, method = 'complete'),
						       dendrogram = 'row',margins = c(12,6),labRow = mature.im.names,
						      ); 

heart.fibro.markers <- 'cardio_fibro_manual_markers.xlsx'
heart.table         <- read.xlsx(heart.fibro.markers, header = TRUE, stringsAsFactors = FALSE, sheetIndex = 2)
fibro.table         <- read.xlsx(heart.fibro.markers, header = TRUE, stringsAsFactors = FALSE, sheetIndex = 1)    
heart.names         <- toupper( heart.table$symbol )
fibro.names         <- toupper( fibro.table$symbol )
heart.heatmap       <- vsd.expr[heart.names[-14],] 
fibro.heatmap       <- vsd.expr[fibro.names,] 
my_palette          <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
heatmap.result      <- heatmap.2( heart.heatmap, col = my_palette, scale  = 'row', 
						          Rowv = TRUE,Colv = FALSE, density.info = 'none',
                                  key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
						          cexRow  = 1, cexCol = 1,srtCol = 30,
                                  distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						          hclustfun  = function(d) hclust(d, method = 'complete'),
						          dendrogram = 'row',margins = c(12,6),labRow = heart.names[-14],
						        ); 
heatmap.result      <- heatmap.2(fibro.heatmap, col = my_palette, scale  = 'row', 
						       Rowv = TRUE,Colv = FALSE, density.info = 'none',
                               key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
						       cexRow  = 1, cexCol = 1,srtCol = 30,
                               distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						       hclustfun  = function(d) hclust(d, method = 'complete'),
						       dendrogram = 'row',margins = c(12,6),labRow = fibro.names,
						      ); 

"
extract fetal.genes from maturation.Rdata
"
## Ino80
## vsd.expr <- vsd.expr.wholeGroup[,c(2:7)]
## Rif1OE
## vsd.expr <- vsd.expr.wholeGroup[,c(1,11,12)]
## Rif1B
## vsd.expr <- vsd.expr.wholeGroup[,c(6,7,8)]
## Rif1C
## vsd.expr <- vsd.expr.wholeGroup[,c(6,7,9,10)]

## load('tissue-specific.Rdata')
##  see code tissue-specific.R
## fetal.genes <- cardiac.gene.names
#  mature.genes <- fibroblast.gene.names
fetal <- intersect(rownames(vsd.expr),toupper(fetal.genes))  
## discarded, variation is 0
##fetal.heatmap       <- vsd.expr[fetal,] 
##heatmap.result      <- heatmap.2(fetal.heatmap, col = my_palette, scale  = 'row', 
##						       Rowv = TRUE,Colv = FALSE, density.info = 'none',
##                               key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
##						       cexRow  = 1, cexCol = 1,srtCol = 30,
##                               distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
##						       hclustfun  = function(d) hclust(d, method = 'complete'),
##						       dendrogram = 'row',margins = c(12,6),labRow = fetal,
##						      ); 
##
mature <- intersect(rownames(vsd.expr),toupper(mature.genes))
mature.heatmap       <- vsd.expr[mature,]
ind    <- apply(mature.heatmap, 1, sd) == 0
subset <- mature.heatmap[!ind,]
heatmap.result      <- heatmap.2( subset, col = my_palette, scale  = 'row', 
						          Rowv = TRUE,Colv = FALSE, density.info = 'none',
                                  key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
						          cexRow  = 1, cexCol = 1,srtCol = 30,
                                  distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						          hclustfun  = function(d) hclust(d, method = 'complete'),
						          dendrogram = 'row',margins = c(12,6),labRow = NA,
						        ); 

fetal <- intersect(rownames(vsd.expr),toupper(fetal.genes))
fetal.heatmap       <- vsd.expr[fetal,]
ind <- apply(fetal.heatmap, 1, sd) == 0
subset <- fetal.heatmap[!ind,]
heatmap.result      <- heatmap.2(subset, col = my_palette, scale  = 'row', 
						         Rowv = TRUE,Colv = FALSE, density.info = 'none',
                                 key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
						         cexRow  = 1, cexCol = 1,srtCol = 30,
                                 distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						         hclustfun  = function(d) hclust(d, method = 'complete'),
						         dendrogram = 'row',margins = c(12,6),labRow = NA,
						        );

"
extract fibro genes/heart genes from 
fibroblast.Rdata
load('fibroblast.Rdata')
"

setwd('/home/zhenyisong/data/cardiodata')
save.image(file = 'ino80.zn.Rdata')
quit("no")   