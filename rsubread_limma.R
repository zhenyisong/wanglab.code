# @author Yisong zhen
# @since  2016-05-26
# @update 2016-08-02
# @name rsubread.R

"
load module 
"
library(Rsubread)
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(DESeq2)
library(gplots)
library(genefilter)
library(xlsx)
library(clusterProfiler)
library(plotrix)
library(pathview)
#---------------------------------------------------------------------------------------
# @parameter
#     input file
#       Fastq raw data:
#            these raw data were orginally from
#            23: VSMC NT-1 day 4
#                NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-23.single.sanger.fastq
#            24: VSMC NT-2 day 4
#                NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-24.single.sanger.fastq
#            25: VSMC H2AZ44-1 day 4 050416
#                NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-25.single.sanger.fastq
#            26:
#               NS50281_160511_NS500489_AHNCGTBGXX.Cnot3-26.single.sanger.fastq
#            index file is saved as 'RNAseq & ChIPseq sample list.xlsx'
#       reference genome sequence
#       target.txt which specify the fastq path
#       output file path
#       base index name
#     output file
#--------------------------------------------------------------------------------------- 


# call the old saved data from image file

setwd('/home/zhenyisong/data/wanglilab/wangcode')
setwd("C:\\Users\\Yisong\\Desktop")
load(file = 'rsubread_limma.Rdata')

"
define the path where the reference genome sequence is deposited
"
genome_ref.path   = "/home/zhenyisong/data/genome/mm10_Bowtie1Index/mm10.fa"

"
this is the targets file which indicate the fastq file path and other experiment
inforamtion regarding the sequence
"
targets.file      = '/home/zhenyisong/data/wanglilab/projects/2016-05-26/targets.txt'
reads.files       = read.table(targets.file,header = F)

"
output path, where the resuls are saved
"
reads.path        = '/home/zhenyisong/data/wanglilab/projects/2016-05-26/'
output.path       = '/home/zhenyisong/data/wanglilab/projects/2016-05-26/results/'


"
generate the path vectors
"
reads.paths       = paste0(reads.path,reads.files$V1)
outputs.files     = paste0(output.path,reads.files$V1,'.sam')

"
the base index name
"
base.string       = 'mm10_ref_index'


"
use the Rsubread command to generate index file
this index file will be generated and saved at getwd()
you do not need to generate the script
"
buildindex( basename = base.string, reference = genome_ref.path )

align( index         = base.string, 
       readfile1     = reads.paths, 
       input_format  = "FASTQ", 
       type          = 'rna',
       output_file   = outputs.files, 
       output_format = "SAM", 
       nthreads      = 8, 
       indels        = 1,
       maxMismatches = 3,
       phredOffset   = 33,
       unique        = T )


# get gene's counts
gene = featureCounts( outputs.files, useMetaFeatures = TRUE, 
                      annot.inbuilt = "mm10", allowMultiOverlap = TRUE)

gene.counts  = gene$counts
gene.ids     = gene$annotation$GeneID
colnames(gene.counts) = c( 'VSMC_H2AZ44_Day_4_1','VSMC_H2AZ44_Day_4_2',
                           'VSMC_NT_Day_4_1','VSMC_NT_Day_4_1');

setwd('/home/zhenyisong/data/wanglilab/wangcode')
save.image(file = 'rsubread_limma.Rdata')

#"
#this code is to transfer data to window system
#and R package needed cannot be installed in Linux HPC
#"
#
#write.table( gene.counts, file = "vsmc.genecount.txt", quote = FALSE, 
#             sep = "\t", row.names = TRUE, col.names = TRUE);
#
## -- start 
#setwd("E:\\FuWai\\王利实验室\\RNA-seq\\2016_05_26")
#gene.df  = read.table( 'vsmc.genecount.txt', header = TRUE, 
#                       sep = "\t", row.names = 1)
#GeneInfo = select( org.Mm.eg.db, keys= as.character(rownames(gene.df)), 
#                   keytype="ENTREZID", columns = c("SYMBOL"));
#m        = match(rownames(gene.df), GeneInfo$ENTREZID);
#Ann      = cbind( rownames(gene.df),
#                  GeneInfo[m, c("SYMBOL")]);
#
#gene.counts = as.matrix(gene.df)
#
## --- you can comment on the above code

keytypes(org.Mm.eg.db)

columns  = c("ENTREZID","SYMBOL", "MGI", "GENENAME");
GeneInfo = select( org.Mm.eg.db, keys= as.character(gene.ids), 
                   keytype="ENTREZID", columns = columns);
m        = match(gene$annotation$GeneID, GeneInfo$ENTREZID);
Ann      = cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                          GeneInfo[m, c("SYMBOL", "MGI", "GENENAME")]);

rownames(gene.counts) = GeneInfo[m,'SYMBOL'];

#write.table( gene.counts, file = "vsmc.counts.txt", quote = FALSE, 
#             sep = "\t", row.names = TRUE, col.names = TRUE);

Ann$Chr  =  unlist( lapply(strsplit(Ann$Chr, ";"), 
                    function(x) paste(unique(x), collapse = "|")))
Ann$Chr  = gsub("chr", "", Ann$Chr)

# if us window, start here
gene.exprs = DGEList(counts = gene.counts, genes = Ann)
A          = rowSums(gene.exprs$counts)
isexpr     = A > 50

hasannot   = rowSums(is.na(gene.exprs$genes)) == 0
gene.exprs = gene.exprs[isexpr & hasannot, , keep.lib.size = FALSE]
gene.exprs = calcNormFactors(gene.exprs)

#---------------------------------------------------------------------
# the following is to generate the heatmap 
# and PCA analysis result
#---------------------------------------------------------------------

dge.tmm                  = t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
#dge.tmm.counts <- round(dge.tmm, digits = 0)
dge.tmm.counts           = apply(dge.tmm,2, as.integer)
rownames(dge.tmm.counts) = gene.exprs$genes$SYMBOL

sample.info              = data.frame( treat  = c('TR','TR','CT','CT') )
dds                      = DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                                   colData   = sample.info,
                                                   design    = ~ treat)
vsd                      = varianceStabilizingTransformation(dds, blind = FALSE);
vsd.expr                 = assay(vsd)
colnames(vsd.expr)       = c('VSMC_H2AZ44_1','VSMC_H2AZ44_2','VSMC_NT_1','VSMC_NT_2')

f1 <- function(x) (IQR(x) > 0.5)
f2 <- function(x) (sd(x)/abs(mean(x)) < 0.1)
ff <- filterfun(f1,f2)
gene.index         <- genefilter(vsd.expr,ff)
heatmap.vsd        <- vsd.expr[gene.index ,]

# Figure 1. group clustering
heatmap( cor(heatmap.vsd),  margins = c(10, 10),
            cexCol = 1, cexRow = 1);


# Figure 2. gene heatmap, gene clustering
# http://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
# The position of each element in the heatmap.2 plot can be 
# controlled using the lmat, lhei and lwid parameters. 
# These are passed by heatmap.2 to the layout command as:
# 
# layout(mat = lmat, widths = lwid, heights = lhei)
# lmat is a matrix describing how the screen is to be broken up. 
# By default, heatmap.2 divides the screen into a four element grid, 
# so lmat is a 2x2 matrix. The number in each element of the matrix 
# describes what order to plot the next four plots in. 
# Heatmap.2 plots its elements in the following order:
# 
#   1. Heatmap,
#   2. Row dendrogram,
#   3. Column dendrogram,
#   4. Key
heatmap.result <- heatmap.2( heatmap.vsd, col = greenred(75),scale  = 'row', 
						     Rowv = TRUE,Colv = FALSE, density.info = 'none',key = TRUE, trace = 'none', 
						     cexCol = 1.5,distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						     hclustfun = function(d) hclust(d, method = 'complete'),
						     dendrogram = 'row',margins = c(12,9),labRow = NA, srtCol = 30,
						     lmat = rbind(c(4,0), c(2,1),c(0,3)), lhei = c(1,3, 0.5), lwid = c(1,4));

# Figure 3. PCA analysis
pr.out <- prcomp(t(heatmap.vsd))
plot( pr.out$x, col = c('blue','blue','brown','brown'), main = 'Principle Component Analysis',
      pch = 18, cex = 2, xlim = c(-25,20))
thigmophobe.labels(pr.out$x[,1],pr.out$x[,2],labels = colnames(heatmap.vsd),cex = 0.7)

pr.var <- pr.out$sdev^2
pve    <- pr.var/sum(pr.var)
plot( seq(1:4),pve, type = 'b', xlab = 'PC number',
      ylab = 'proportion of variance')
#---------------------------------------------------------------------
# PCA and heatmap END
# END
#---------------------------------------------------------------------

d      <- gene.exprs

group  = factor(c('CT','CT','TR','TR'));
design = model.matrix(~ 0 + group);
colnames(design) = c('Control','Treatment')
contrast.matrix  = makeContrasts(Treatment - Control, levels = design)
d.norm          = voom(d, design = design)
fit             = lmFit(d.norm, design)
fit2            = contrasts.fit(fit,contrast.matrix)
fit2            = eBayes(fit2)


gene.result = topTable(  fit2, 
                         number        = Inf, 
                         adjust.method = "BH", 
                         sort.by       = "p");


setwd("E:\\CardioSignal\\publicData")
file.name <- 'transfac.name.final.table'
tf.whole.set <- read.table(file.name, header = FALSE)
tf.exprs.log <- vsd.expr[rownames(vsd.expr) %in% tf.whole.set$V1,]
tf.gene.result <- gene.result[gene.result$SYMBOL %in% tf.whole.set$V1, ]
setwd("C:\\Users\\Yisong\\Desktop")
write.xlsx( gene.result, file = "vsmc_all.xlsx", row.names = TRUE, col.names = TRUE);
write.xlsx( tf.exprs.log, file = "tf_exprs_all.xlsx", row.names = TRUE, col.names = TRUE);
write.xlsx(tf.gene.result, file = 'tf.dge.xlsx',row.names = TRUE, col.names = TRUE);

##gene.result = topTable( fit2, number = Inf, 
##                        adjust.method="BH", sort.by="p",
##                        lfc = 0.58, p.value = 0.05);
##
#write.table( gene.result, file = "vsmc.xls", quote = FALSE, 
#             sep = "\t", row.names = TRUE, col.names = TRUE);

write.xlsx(file = 'vsmc.xls',gene.result)


#-------------------------------------------------------------
# Figure KEGG enrichment analysis
# Figure Go enrichment analysis
#-------------------------------------------------------------

gene.entrez.id = gene.result$GeneID

kegg.table  =  enrichKEGG( gene.entrez.id, organism = "mouse", 
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH", 
                           qvalueCutoff  = 0.1, readable = TRUE)
kegg.result       = summary(kegg.table)
kegg.qvalue       = -log(kegg.result$qvalue)
kegg.pathway.name = kegg.result$Description

#barplot(kegg.qvalue, names.arg = kegg.pathway.name, las = 3.5)

par(mar = c(12,4,1,1), fin = c(4,4))

x = barplot( kegg.qvalue, cex.lab = 0.8,cex.axis= 0.8,
             main = 'KEGG enrichment anlysis', cex.main = 0.8,
             ylab = '-log(q-value of enrichment)')
text( cex = 0.6, x = x - 0.25, y = -1.25, 
      kegg.pathway.name, 
      xpd = TRUE, srt = 50, pos = 2)

# this Figure is from the code 
# http://www.bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
dotplot(kegg.table)
enrichMap(kegg.table)


go.table = enrichGO( gene.entrez.id, organism = "mouse",
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
go.cc.result      = summary(go.table)
go.qvalue         = -log(go.cc.result$qvalue)
go.term.name      = go.cc.result$Description


par(mar = c(6,12,1,1), fin = c(4.5,4.5))
x = barplot( go.qvalue[1:10], cex.lab = 0.8,cex.axis= 0.8,
             names.arg = go.term.name[1:10], horiz = TRUE,
             main = 'GO_MF enrichment anlysis', cex.main = 0.8,
             las  = 2, cex.name = 0.8,
             xlab = '-log(q-value of GO enrichment)')


go.table = enrichGO( gene.entrez.id, organism = "mouse",
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
go.cc.result      = summary(go.table)
go.qvalue         = -log(go.cc.result$qvalue)
go.term.name      = go.cc.result$Description


par(mar = c(6,12,1,1), fin = c(4.5,4.5))
x = barplot( go.qvalue[1:10], cex.lab = 0.8,cex.axis= 0.8,
             names.arg = go.term.name[1:10], horiz = TRUE,
             main = 'GO_CC enrichment anlysis', cex.main = 0.8,
             las  = 2, cex.name = 0.8,
             xlab = '-log(q-value of GO enrichment)')


go.table = enrichGO( gene.entrez.id, organism = "mouse",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
go.cc.result      = summary(go.table)
go.qvalue         = -log(go.cc.result$qvalue)
go.term.name      = go.cc.result$Description


par(mar = c(6,12,1,1), fin = c(4.5,4.5))
x = barplot( go.qvalue[1:10], cex.lab = 0.8,cex.axis= 0.8,
             names.arg = go.term.name[1:10], horiz = TRUE,
             main = 'GO_BP enrichment anlysis', cex.main = 0.8,
             las  = 2, cex.name = 0.75,
             xlab = '-log(q-value of GO enrichment)')

do.table <- enrichDO( gene.entrez.id,
                     ont = "DO",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
pathway.genelist = gene.result$logFC
names(pathway.genelist) = gene.result$GeneID
hippo <- pathview(gene.data  = pathway.genelist,
                     pathway.id = "mmu04390",
                     species    = "mmu",
                     limit      = list(gene = max(abs(pathway.genelist)), cpd=1))

## End(Not run)