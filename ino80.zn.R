# ino80.zn.R
# raw data is here from wangli Excel annotation
# I saved the following inforamtion in the files.txt
# NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-15.single.sanger.fastq CM_Rif1OE_1.fastq
# NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-16.single.sanger.fastq CM_NT_1.fastq
# NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-17.single.sanger.fastq CM_Ino80_1_1.fastq
# NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-18.single.sanger.fastq CM_Ino80_3_1.fastq
# NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-19.single.sanger.fastq CM_Rif1C_1.fastq
# NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-1.single.sanger.fastq  CM_NT_2.fastq
# NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-20.single.sanger.fastq CM_GFPOE_1.fastq
# NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-22.single.sanger.fastq CM_Rif1OE_2.fastq
# NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-2.single.sanger.fastq  CM_Ino80_1_2.fastq
# NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-3.single.sanger.fastq  CM_Ino80_3_2.fastq
# NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-4.single.sanger.fastq  CM_Rif1B_1.fastq
# NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-5.single.sanger.fastq  CM_Rif1C_2.fastq

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

# the processed bam files are needed.
# @parent
#    /home/zhenyisong/data/wanglilab/projects/2016-06-07
#        rsubread_limma.R

#--------------------------------------------------------------------------
#  the following code is to extract gene count information from
#  processed sam files.
#  the code is extracted from parent program
#--------------------------------------------------------------------------


"
this is the targets file which indicate the fastq file path and other experiment
inforamtion regarding the sequence
"
targets.file      = '/home/zhenyisong/data/wanglilab/projects/2016-06-07/targets.txt'
reads.files       = read.table(targets.file,header = F)

"
output path, where the resuls are saved
"
reads.path        = '/home/zhenyisong/data/wanglilab/projects/2016-06-07/'
output.path       = '/home/zhenyisong/data/wanglilab/projects/2016-06-07/rsubread_results/'

"
generate the path vectors
"
reads.paths       = paste0(reads.path,reads.files$V1)
outputs.files     = paste0(output.path,reads.files$V1,'.sam')

# get gene's counts
gene         <-  featureCounts( outputs.files, useMetaFeatures = TRUE, 
                        annot.inbuilt = "hg19", allowMultiOverlap = TRUE)
gene.counts  <- gene$counts

colnames(gene.counts)
gene.ids     <- gene$annotation$GeneID
colnames(gene.counts) <- c( 'CM_GFPOE_1','CM_Ino80_1_1','CM_Ino80_1_2',
                            'CM_Ino80_3_1','CM_Ino80_3_2','CM_NT_1',
                            'CM_NT_2','CM_Rif1B_1','CM_Rif1C_1',
                            'CM_Rif1C_2','CM_Rif1OE_1','CM_Rif1OE_2');
gene.counts   <- gene.counts[,c(2:7)]

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

group           <- factor(c( 'Ino80_1','Ino80_1','Ino80_3',
                             'Ino80_3','CM_NT','CM_NT'));
sample.info     <- data.frame(treat = group)
dds             <- DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                           colData   = sample.info,
                                           design    = ~ treat)
vsd                     <- varianceStabilizingTransformation(dds, blind = FALSE);
vsd.expr                <- assay(vsd)
rownames(vsd.expr)      <- GeneInfo[m,'SYMBOL']
colnames(vsd.expr)      <- colnames(gene.counts) 

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
fetal <- intersect(rownames(vsd.expr),toupper(fetal.genes))
fetal.heatmap       <- vsd.expr[fetal,] 
heatmap.result      <- heatmap.2(fetal.heatmap, col = my_palette, scale  = 'row', 
						       Rowv = TRUE,Colv = FALSE, density.info = 'none',
                               key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
						       cexRow  = 1, cexCol = 1,srtCol = 30,
                               distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						       hclustfun  = function(d) hclust(d, method = 'complete'),
						       dendrogram = 'row',margins = c(12,6),labRow = fetal,
						      ); 

mature <- intersect(rownames(vsd.expr),toupper(mature.genes))
mature.heatmap       <- vsd.expr[mature,]
ind <- apply(mature.heatmap, 1, sd) == 0
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