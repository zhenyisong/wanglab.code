library(gplots)
library(xlsx)
library(Rsubread)
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(DESeq2)
library(gplots)
library(genefilter)

setwd('/home/zhenyisong/data/wanglilab/vsmc_db');

"
these are processed high-through-put data, those 
data needed the Perl script to extract the 
information. the data were inputed by read.table method
"
rna.seq.filename = 'final_rna_seq.cos';
non.db.filename  = 'final_nons.cos';
affy.db.filename = 'final_affy.cos';
rna.seq.db = read.table( rna.seq.filename, header = TRUE, sep = "\t",
                         row.names = 1)
non.db     = read.table( non.db.filename , header = TRUE, sep = "\t",
                         row.names = 1)
affy.db    = read.table( affy.db.filename , header = TRUE, sep = "\t",
                         row.names = 1)
"
transform the data.frame into the log format
"
rna.matrix     = as.matrix(rna.seq.db)
rna.log.matrix = log(rna.matrix + 1)
affy.matrix    = as.matrix(affy.db)

"
the following module is to extract the differential expression
gene name from VSMC data analysis.
see the whole script rsubread_limma.R
get the gene names from DEG analysis
"



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
gene.result = topTable( fit2, number  = Inf, 
                        adjust.method = "BH", 
                        sort.by = "p",
                        lfc     = 0.58
                        p.value = 0.05);
"
get the DEG gene names and transform the 
gene name to upper case.
"
deg.names   = rownames(gene.result)
deg.names   = toupper(deg.names)
#
# -- module end
#

"
get the common name
and gene expression matrix
"
common.name          = intersect(rownames(rna.matrix), rownames(affy.db))
common.name          = intersect(common.name, deg.names)
colnames.vector      = c(colnames(rna.matrix),colnames(affy.matrix))
common.matrix        = seq(1:length(colnames.vector))

for( gene in common.name) {
    gene.exprs    = matrix( c( rna.log.matrix[gene,],
                                affy.matrix[gene,]),
                            byrow = F, nrow = 1)
    common.matrix = rbind(common.matrix,gene.exprs)
}
length(common.name)
exprs.matrix            = common.matrix[-1,]
colnames(exprs.matrix)  = colnames.vector

results = cor(exprs.matrix,method = 'spearman')

#rna.cor = cor(rna.log.matrix, method = 'spearman')


#----------------------------------------------
# random shuffle the matrix and generate
# the Spearman p.value distribution
# simulation begin
#----------------------------------------------

pseudo_num = seq(288 * 288)
random.matrix = exprs.matrix
shuffle.times = 100
for( j in 1:shuffle.times) {
    for ( i in 1:dim(exprs.matrix)[1]) {
        random.matrix[1,] = exprs.matrix[1,sample(dim(exprs.matrix)[2])]
    }
    pseudo.cor = cor(random.matrix, method = 'spearman')
    pseudo_num = cbind(pseudo_num,as.vector(pseudo.cor))
}
pseudo_num = pseudo_num[,-1]

pdf("pseudo.pvalue.pdf")
hist(as.vector(pseudo_num))
dev.off()
cutoff   = mean(as.vector(pseudo_num) > 0.8)
fileConn = file("cutoff.txt")
writeLines(c("hello","world",cutoff), fileConn)
close(fileConn)

# ---simulation end

whole.heatmap = heatmap( results,  margins = c(10, 10),
                         cexCol = 0.2, cexRow = 0.2);
partial.map = results[results['SRR01',] > 0.9,results['SRR01',] > 0.9]
partial.map = results[results['SRR01',] > 0.62,results['SRR01',] > 0.62]
write.xlsx(partial.map, file = "cutoff_0.8.xls")
heatmap(  partial.map,  margins = c(10, 10),
           cexCol = 1, cexRow = 1);
rownames(results)[results['SRR01',] > 0.9]
colnames.vector[whole.heatmap$rowInd]
summary(as.vector(rna.cor))
spearman.d = as.vector(rna.cor)
hist(spearman.d, prob = TRUE, n = 200, col = 'grey')
lines(density(spearman.d), col = "blue", lwd = 2) # add a density estimate with defaults
lines(density(spearman.d, adjust=2), lty = "dotted", col = "darkgreen", lwd = 2) 

# Figure
# heatmap of smooth muscle differentiation markers

vcms.markers  = 'SM-markers.xlsx'
vcms.table    = read.xlsx(vcms.markers,header = TRUE, stringsAsFactors = FALSE, 1)
vsmc.genename = vcms.table$GeneSymbol

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
vsd.markers              = vsd.expr[vsmc.genename,]

# this code is extracted from 
# http://stackoverflow.com/questions/17820143/how-to-change-heatmap-2-color-range-in-r
colors      = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))
my_palette  =  colorRampPalette(c("red", "black", "green"))(n = 299)
heatmap.result = heatmap.2(vsd.markers, col = my_palette,scale  = 'row', 
						   Rowv = TRUE,Colv = FALSE, density.info = 'none',
                           key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
						   cexRow = 1.5, cexCol = 1.5,srtCol = 30,
                           distfun= function(d) as.dist(1-cor(t(d),method = 'pearson')),
						   hclustfun  = function(d) hclust(d, method = 'complete'),
						   dendrogram = 'row',margins = c(12,6),labRow = vsmc.genename,
						  );
