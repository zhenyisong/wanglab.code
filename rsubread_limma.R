# @author Yisong zhen
# @since  2016-05-26
# @update 2016-05-26
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

#-------------------------------------------------------
# @parameter
#     input file
#       reference genome sequence
#       target.txt which specify the fastq path
#       output file path
#       base index name
#     output file
#-------------------------------------------------------    

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


keytypes(org.Mm.eg.db)

columns  = c("ENTREZID","SYMBOL", "MGI", "GENENAME");
GeneInfo = select( org.Mm.eg.db, keys= as.character(gene.ids), 
                   keytype="ENTREZID", columns = columns);
m        = match(gene$annotation$GeneID, GeneInfo$ENTREZID);
Ann      = cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                          GeneInfo[m, c("SYMBOL", "MGI", "GENENAME")]);

rownames(gene.counts) = GeneInfo[m,'SYMBOL'];
write.table( gene.counts, file = "vsmc.counts.txt", quote = FALSE, 
             sep = "\t", row.names = TRUE, col.names = TRUE);

Ann$Chr  =  unlist( lapply(strsplit(Ann$Chr, ";"), 
                    function(x) paste(unique(x), collapse = "|")))
Ann$Chr  = gsub("chr", "", Ann$Chr)

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
heatmap(   cor(heatmap.vsd),  margins = c(10, 10),
            cexCol = 1, cexRow = 1);


# Figure 2. gene heatmap, gene clustering
heatmap.result<-heatmap.2(heatmap.vsd, col=greenred(75),scale = 'row', 
						 Rowv = TRUE,Colv=FALSE, density.info='none',key=FALSE, trace='none', 
						 cexRow=0.1,distfun= function(d) as.dist(1-cor(t(d),method = 'pearson')),
						 hclustfun = function(d) hclust(d, method = 'complete'),
						 dendrogram = 'row',margins =c(12,9),labRow = NA,
						 lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 10, 0.25 ));

# Figure 3. PCA analysis
pr = prcomp(t(heatmap.vsd))
plot( pr$x, col = 'white', main = 'PC plot',
      xlim = c(-22,15))
text(pr$x[,1],pr$x[,2],labels = colnames(heatmap.vsd),cex = 0.7)

#---------------------------------------------------------------------
# PCA and heatmap END
# END
#---------------------------------------------------------------------

d          = gene.exprs

group  = factor(c('CT','CT','TR','TR'));
design = model.matrix(~ 0 + group);
colnames(design) = c('Control','Treatment')
contrast.matrix  = makeContrasts(Treatment - Control, levels = design)
d.norm          = voom(d, design = design)
fit             = lmFit(d.norm, design)
fit2            = contrasts.fit(fit,contrast.matrix)
fit2            = eBayes(fit2)


#gene.result = topTable( fit2, coef = ncol(design), 
#                        number = Inf, adjust.method="BH", sort.by="p");
gene.result = topTable( fit2, number = Inf, 
                        adjust.method="BH", sort.by="p",
                        lfc = 0.58, p.value = 0.05);

#write.table( gene.result, file = "vsmc.xls", quote = FALSE, 
#             sep = "\t", row.names = TRUE, col.names = TRUE);

write.xlsx(file = 'vsmc.xls',gene.result)


## End(Not run)