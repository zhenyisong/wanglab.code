library(gplots)
setwd('/home/zhenyisong/data/wanglilab/vsmc_db');
rna.seq.filename = 'final_rna_seq.cos';
non.db.filename  = 'final_nons.cos';
affy.db.filename = 'final_affy.cos';
rna.seq.db = read.table( rna.seq.filename, header = TRUE, sep = "\t",
                         row.names = 1)
non.db     = read.table( non.db.filename , header = TRUE, sep = "\t",
                         row.names = 1)
affy.db    = read.table( affy.db.filename , header = TRUE, sep = "\t",
                         row.names = 1)
rna.matrix     = as.matrix(rna.seq.db)
rna.log.matrix = log(rna.matrix + 1)
affy.matrix    = as.matrix(affy.db)

common.name          = intersect(rownames(rna.matrix), rownames(affy.db))
colnames.vector      = c(colnames(rna.matrix),colnames(affy.matrix))
common.matrix = seq(1:length(colnames.vector))

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
rna.cor = cor(rna.log.matrix, method = 'spearman')
whole.heatmap = heatmap( results,  margins = c(10, 10),
                         cexCol = 0.2, cexRow = 0.2);
partial.map = results[results['SRR01',] > 0.9,results['SRR01',] > 0.9]
heatmap(  partial.map,  margins = c(10, 10),
           cexCol = 1, cexRow = 1);
rownames(results)[results['SRR01',] > 0.9]
colnames.vector[whole.heatmap$rowInd]
summary(as.vector(rna.cor))
spearman.d = as.vector(rna.cor)
hist(spearman.d, prob = TRUE, n = 200, col = 'grey')
lines(density(spearman.d), col = "blue", lwd = 2) # add a density estimate with defaults
lines(density(spearman.d, adjust=2), lty = "dotted", col = "darkgreen", lwd = 2) 
