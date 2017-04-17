library(gplots)
library(xlsx)
library(Rsubread)
library(edgeR)
library(limma)
library(DESeq2)
library(genefilter)
library(RColorBrewer)
library(ggplot2)
library(cluster)
library(factoextra)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(GEOquery)
library(affy)
library(annotate)
library(mgu74av2.db)
library(DOSE)

# where is the raw data?
# cleaned raw data is deposited in 

# cd /home/zhenyisong/wanglab/wangcode
# nohup R CMD BATCH zj.liangping.R &
genome_ref.path      <- "/bioware/genome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
setwd('/home/zhenyisong/wanglab/wangdata/lianglab/Rawdata')
reads.files.names    <- list.files(pattern = '*.fastq')
read.path.1          <- reads.files.names[grep("R1",reads.files.names)]
read.path.2          <- reads.files.names[grep("R2",reads.files.names)]
reads.paths.1        <- paste0(getwd(),'/',read.path.1)
reads.paths.2        <- paste0(getwd(),'/',read.path.2)
setwd('/home/zhenyisong/wanglab/wangdata/lianglab')
unlink('rsubread')
dir.create('rsubread')
output.path          <- '/home/zhenyisong/wanglab/wangdata/lianglab/rsubread'
setwd('/home/zhenyisong/wanglab/wangdata/lianglab/rsubread')

outputs.files        <- paste0(output.path,'/', read.path.1,'.bam')

base.string          <-  'hg19_index'
buildindex( basename = base.string, reference = genome_ref.path )
setwd(output.path)
align( index          = base.string, 
       readfile1      = reads.paths.1, 
       readfile2      = reads.paths.2, 
       input_format   = "FASTQ", 
       type           = 'rna',
       output_file    = outputs.files, 
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 4, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

liangping.genes      <- featureCounts( outputs.files, useMetaFeatures = TRUE,
                                       countMultiMappingReads = FALSE,
                                       strandSpecific         = 0, 
                                       isPairedEnd            = TRUE,
                                       requireBothEndsMapped  = TRUE,
                                       autosort               = TRUE,
                                       nthreads               = 8,
                                       annot.inbuilt = "hg19", allowMultiOverlap = TRUE)

#
# the temperary file transfer

# setwd("/home/zhenyisong/wanglab/wangdata/lianglab")
# save.image(file = 'liangping.Rdata')
# quit("no")
# setwd("D:\\wangli_data\\Rdata")
# load("liangping.Rdata")
#


gene         <- liangping.genes
gene.counts  <- gene$counts
gene.ids     <- gene$annotation$GeneID

keytypes(org.Hs.eg.db)

columns  <- c("ENTREZID","SYMBOL", "GENENAME");
GeneInfo <- select( org.Hs.eg.db, keys= as.character(gene.ids), 
                   keytype = "ENTREZID", columns = columns);
m        <- match(gene$annotation$GeneID, GeneInfo$ENTREZID);
Ann      <- cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                          GeneInfo[m, c("SYMBOL", "GENENAME")]);

Ann$Chr  <- unlist( lapply(strsplit(Ann$Chr, ";"), 
                    function(x) paste(unique(x), collapse = "|")))
Ann$Chr  <- gsub("chr", "", Ann$Chr)

gene.exprs <- DGEList(counts = gene.counts, genes = Ann)
gene.exprs <- calcNormFactors(gene.exprs)
dge.tmm    <- t(t(gene.exprs$counts) * gene.exprs$samples$norm.factors)
dge.tmm.counts           <- apply(dge.tmm,2, as.integer)
sample.info              <- data.frame( treat  = c('Control','Treat',
                                                   'Control','Treat',
                                                   'Control','Treat') )
dds                      <- DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                                    colData   = sample.info,
                                                    design    = ~ treat)
vsd                      <- varianceStabilizingTransformation(dds);
vsd.expr                 <- assay(vsd)
rownames(vsd.expr)       <- gene.exprs$genes$SYMBOL
# setwd("C:\\Users\\Yisong\\Desktop")
# write.table(vsd.expr, file = 'zhejiang.temp.table')
#rownames(vsd.expr)       <- gene.exprs$genes$SYMBOL
rownames(vsd.expr)       <- NULL
colnames(vsd.expr)       <- c('Control-1','Treat-1',
                              'Control-2','Treat-2',
                              'Control-3','Treat-3')
sds <- rowSds(vsd.expr)
sh  <- shorth(sds)
vsd.filtered.expr <- vsd.expr[sds > 0.3,]

heatmap(cor(vsd.filtered.expr, method = 'spearman'), cexCol = 0.8, cexRow = 0.8)


data.dist  <- as.dist( 1 - cor(vsd.filtered.expr,method = 'spearman'))
cluster.hc <- hclust( d = data.dist, method = 'complete')
groups     <- cutree(cluster.hc, k = 2)
# setwd("D:\\wangli_data\\Rdata")
# save.image('liangping.Rdata')
fviz_cluster( list( data = t(vsd.filtered.expr), cluster = groups),
              palette = c("red", 'blue') ,
              ellipse.type = "convex",
              repel = FALSE,
              show.clust.cent = FALSE,
              ggtheme = theme_minimal()  )

# --- the above result confirmed the expriment design is OK

gene.tmm   <- DGEList(counts = gene.exprs, genes = Ann)
cell.lines <- factor(rep(c('C1','C2','C3'),each = 2), levels = c('C1','C2','C3'))
groups     <- factor( c('Control','Treatment','Control','Treatment','Control','Treatment'), 
                      levels = c('Control','Treatment'));
design     <- model.matrix(~ 0 + groups);
colnames(design) <- c('Control','Treatment')
d.norm           <- voom(gene.tmm, design = design)
cor.fit          <- duplicateCorrelation(d.norm,design,block = cell.lines)
#https://support.bioconductor.org/p/46936/
#the follwoing value should be positive
#cor.fit$consensus
fit <- lmFit(d.norm,design,block = cell.lines,correlation = cor.fit$consensus)

contrast.matrix  <- makeContrasts(Treatment - Control, levels = groups)

fit2             <- contrasts.fit(fit,contrast.matrix)
fit2             <- eBayes(fit2)

gene.result      <- topTable(  fit2, 
                               number        = Inf, 
                               adjust.method = "BH", 
                               sort.by       = "p");
gene.result.dge  <- topTable(  fit2, 
                               number        = Inf, 
                               adjust.method = "BH", 
                               sort.by       = "p",
                               p.value       = 0.05 );


gene.entrez.id <- gene.result.dge$GeneID

kegg.table     <- enrichKEGG( gene.entrez.id, organism = "human", 
                              pvalueCutoff  = 0.05, 
                              pAdjustMethod = "BH", 
                              qvalueCutoff  = 0.1)
kegg.result       <- summary(kegg.table)
kegg.qvalue       <- -log(kegg.result$qvalue)
kegg.pathway.name <- kegg.result$Description



par(mar = c(12,4,1,1), fin = c(4,4))

x = barplot( kegg.qvalue, cex.lab = 0.8,cex.axis= 0.8, width = 0.5,
             main = 'KEGG enrichment anlysis', cex.main = 0.8,
             ylab = '-log(q-value of enrichment)')
text( cex = 0.75, x = x - 0.25, y = -1.25, 
      kegg.pathway.name, 
      xpd = TRUE, srt = 60, pos = 2)


# setwd("C:\\Users\\Yisong\\Desktop")
# write.xlsx(gene.result.dge, file = 'zhejiang.xlsx')
pathway.genelist <- gene.result$logFC
names(pathway.genelist) <-  gene.result$GeneID
akt.pathway <- pathview(gene.data  = pathway.genelist,
                        pathway.id = "hsa04151",
                        species    = "hsa",
                        limit      = list(gene = max(abs(pathway.genelist)), cpd = 1))

#
# to do list
#---

##go.table = enrichGO( gene.entrez.id, organism = "mouse",
##                     ont = "MF",
##                     pAdjustMethod = "BH",
##                     pvalueCutoff  = 0.01,
##                     qvalueCutoff  = 0.05)
##go.cc.result      = summary(go.table)
##go.qvalue         = -log(go.cc.result$qvalue)
##go.term.name      = go.cc.result$Description
##
##
##par(mar = c(6,12,1,1), fin = c(4.5,4.5))
##x = barplot( go.qvalue[1:10], cex.lab = 0.8,cex.axis= 0.8,
##             names.arg = go.term.name[1:10], horiz = TRUE,
##             main = 'GO_MF enrichment anlysis', cex.main = 0.8,
##             las  = 2, cex.name = 0.8,
##             xlab = '-log(q-value of GO enrichment)')
##
##
##go.table = enrichGO( gene.entrez.id, organism = "human",
##                     ont = "CC",
##                     pAdjustMethod = "BH",
##                     pvalueCutoff  = 0.01,
##                     qvalueCutoff  = 0.05)
##go.cc.result      = summary(go.table)
##go.qvalue         = -log(go.cc.result$qvalue)
##go.term.name      = go.cc.result$Description
##
##
##par(mar = c(6,12,1,1), fin = c(4.5,4.5))
##x = barplot( go.qvalue[1:10], cex.lab = 0.8,cex.axis= 0.8,
##             names.arg = go.term.name[1:10], horiz = TRUE,
##             main = 'GO_CC enrichment anlysis', cex.main = 0.8,
##             las  = 2, cex.name = 0.8,
##             xlab = '-log(q-value of GO enrichment)')
##

##
## discard the raw data because the liang'ping side
## use the PI3K pathway to infer the transcription network
## instead of AKT1 signaling pathway
##
##setwd("E:\\FuWai\\wangli.lab\\liangping")
##gse     <- getGEO(filename = "GSE3383_family.soft.gz")
##gsmlist <- GSMList(gse)
##gpl     <- GPLList(gse)
##
##probesets   <- Table(GPLList(gse)[[1]])$ID
##data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
##                                      {tab <- Table(x)
##                                       mymatch <- match(probesets,tab$ID_REF)
##                                       return(tab$VALUE[mymatch])
##                                     }))
##data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
##data.matrix <- log2(data.matrix)
##rownames(data.matrix) <- probesets
##colnames(data.matrix) <- names(gsmlist)
##pdata <- data.frame(samples = names(gsmlist))
##rownames(pdata) <- names(gsmlist)
##pheno <- as(pdata,"AnnotatedDataFrame")
##eset2 <- new('ExpressionSet',exprs = data.matrix,phenoData = pheno)


## this data set contains the gain-of-function
## loss-of-funtion results from PI3K pathway
## 
setwd("E:\\FuWai\\wangli.lab\\liangping\\GSE558_RAW")
raw.data    <- ReadAffy();
rma.data    <- rma(raw.data);
exprs.data  <- exprs(rma.data)
probes      <- rownames(exprs.data)
gene.symbol       <- unlist(mget(probes, mgu74av2SYMBOL, ifnotfound = NA))
rownames(exprs.data) <- gene.symbol
colnames(exprs.data) <- c( 'PI3Kca-1','PI3Kca-2','PI3Kca-2',
                           'PI3Kdn-1','PI3Kdn-2','PI3Kdn-3',
                           'PI3Kntg-1','PI3Kntg-2','PI3Kntg-3')
groups              <- factor( c( rep(1,3),rep(2,3),rep(3,3) ), levels = c(1:3),
                               labels = c("PI3Kca","PI3Kdn","PI3Kntg") )
design              <- model.matrix(~ 0 + groups)
colnames(design)    <- levels(groups)
fit                 <- lmFit(exprs.data, design)
contrast.matrix     <- makeContrasts(PI3Kca - PI3Kntg, PI3Kdn - PI3Kntg, levels = design)
fit2                <- contrasts.fit(fit, contrast.matrix)
fit2                <- eBayes(fit2)
PI3Kca.topTable     <- topTable(  fit2, coef = 1, adjust = "BH", 
                                  number = Inf, sort.by = "P", p.value = 0.05)
PI3Kdn.topTable     <- topTable(  fit2, coef = 2, adjust = "BH", 
                                  number = Inf, sort.by = "P", p.value = 0.05)

PI3K.gene.symbol.set <- na.omit(union(PI3Kca.topTable$ID,PI3Kdn.topTable$ID))
#PI3K.gene.symbol.set <- na.omit(PI3Kdn.topTable$ID)
human.PI3K.set       <- toupper(PI3K.gene.symbol.set)
human.PI3K.id.set    <- select( org.Hs.eg.db, keys= as.character(human.PI3K.set), 
                                keytype = "SYMBOL", columns = 'ENTREZID');
human.PI3K.id.set    <- na.omit(human.PI3K.id.set)
PI3K.data.list       <- data.frame( diseaseId = unlist(rep('PI3K', length(human.PI3K.id.set))), 
                                     geneId   = unlist(human.PI3K.id.set) )
PID_PI3KCI_PATHWAY.MSigDB <- c( 'ADAP1','ARAP3','ARF1','ARF5','ARF6','BLK','BLNK',
                                'BTK','CYTH1','CYTH2','CYTH3','DAPP1','FGR','FOXO3',
                                'FYN','HCK','HRAS','HSP90AA1','INPP5D','INPPL1','ITK',
                                'KRAS','LAT','LCK','LYN','NRAS','PDPK1','PIK3CA','PIK3CB',
                                'PIK3CD','PIK3CG','PIK3R1','PIK3R2','PIK3R3','PIK3R5','PIK3R6',
                                'PLCG1','PLCG2','PLEKHA1','PLEKHA2','PTEN','RAC1','RAP1A','RHOA',
                                'SGK1','SRC','SYK','YES1','ZAP70')
PI3K.MSigDB.id.set           <- select( org.Hs.eg.db, keys= as.character(PID_PI3KCI_PATHWAY.MSigDB), 
                                        keytype = "SYMBOL", columns = 'ENTREZID');
PI3K.MSigDB.id.set           <- data.frame( diseaseId = unlist(rep('PI3K', length(PI3K.MSigDB.id.set))), 
                                             geneId   = unlist(PI3K.MSigDB.id.set) )
PI3K.pathway.genelist        <- gene.result.dge$logFC
names(PI3K.pathway.genelist) <- gene.result.dge$GeneID
PI3K.pathway.genelist        <- sort(PI3K.pathway.genelist, decreasing = TRUE)
gsea.result                  <- GSEA(PI3K.pathway.genelist,TERM2GENE = PI3K.data.list, maxGSSize = 1000 )
gsea.result                  <- GSEA(PI3K.pathway.genelist,TERM2GENE = PI3K.MSigDB.id.set, maxGSSize = 1000 )

## test GSEA
# http://www.disgenet.org/ds/DisGeNET/results/all_gene_disease_associations.tsv.gz
# DisGeNET which contains 381056 associations, between 16666 genes and 13172 diseases 

setwd('C:\\Users\\Yisong\\Desktop\\all_gene_disease_associations.tsv')
gda <- read.delim('all_gene_disease_associations.tsv', header = TRUE,  comment.char = "#")
disease2gene <- gda[, c("diseaseId", "geneId")]
disease2name <- gda[, c("diseaseId", "diseaseName")]
set.seed(123)
zhejiang.GSEA <- GSEA(PI3K.pathway.genelist, TERM2GENE = disease2gene, TERM2NAME = disease2name) 
#colnames(summary(zhejiang.GSEA))
gsea.result <- summary(zhejiang.GSEA)
gsea.result$Des
gseaplot(zhejiang.GSEA,'umls:C0151744')

