library(affy)
library(annotate)
library(rat2302.db)
library(limma)
library(xlsx)
library(GEOquery)
"
set the dir to read the raw data GSE19106
"
setwd("D:\\wangli_data\\GSE19106")

raw.data    <- ReadAffy();
rma.data    <- rma(raw.data);
exprs.data  <- exprs(rma.data)

"
please see the Limma manual at Page 37
"

design           <- model.matrix(~ 0 + factor(c(1,1,2,2)))
colnames(design) <- c("control", "treat")
fit              <- lmFit(exprs.data, design)
contrast.matrix  <- makeContrasts(treat - control, levels = design)
fit2             <- contrasts.fit(fit, contrast.matrix)
fit2             <- eBayes(fit2)


probes           <- rownames(exprs.data)
symbols          <- unlist(mget(probes, rat2302SYMBOL, ifnotfound = NA))

result.limma     <- topTable( fit2, coef = 1, adjust = "BH",
                              sort.by="B", p.value = 0.05, 
                              lfc = 0.58, genelist = symbols,
                              number = Inf)

write.xlsx(file = 'GSE19106.dzy.xls',result.limma)

"
I downloaded the raw data from NCBI ftp

"
setwd("D:\\wangli_data\\GSE13836")
gse    <- getGEO(filename = 'GSE13836_series_matrix.txt.gz')
cutoff <- 0.58
str(gse)
vsmc.exprs <- exprs(gse)
rownames(vsmc.exprs) <- featureData(gse)@data$GENE_SYMBOL
vmsc.lfc          <- vsmc.exprs[,1] - vsmc.exprs[,2]
vmsc.lfc.final    <- vmsc.lfc[abs(vmsc.lfc) > cutoff] 
vmsc.df           <- as.data.frame(vmsc.lfc.final)
vmsc.df$genename  <- names(vmsc.lfc.final)
write.xlsx(file = 'GSE13836.dzy.xls',vmsc.df)
