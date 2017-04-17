library(affy)
library(annotate)
library(rat2302.db)
library(limma)
library(xlsx)
library(GEOquery)
library(gplots)
library(Rsubread)
library(edgeR)
library(org.Hs.eg.db)
library(DESeq2)
library(gplots)
library(genefilter)

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

result.nofilter.limma <- topTable( fit2, coef = 1, adjust = "BH",
                                   sort.by="B",genelist = symbols,
                                   number = Inf)

write.xlsx(file = 'GSE19106.dzy.xls',result.nofilter.limma)

"
I downloaded the raw data from NCBI ftp

"
setwd("D:\\wangli_data\\GSE13836")
gse    <- getGEO(filename = 'GSE13836_series_matrix.txt.gz')
cutoff <- 0.58
str(gse)
vsmc.exprs <- exprs(gse)
boxplot(vsmc.exprs)
vmsc.lfc             <- vsmc.exprs[,1] - vsmc.exprs[,2]
vsmc.df              <- data.frame( gene.name = featureData(gse)@data$ORF,
                                    log.exprs = vmsc.lfc)
#vmsc.lfc.final    <- vmsc.lfc[abs(vmsc.lfc) > cutoff] 
#vmsc.df           <- as.data.frame(vmsc.lfc.final)
#vmsc.df$genename  <- names(vmsc.lfc.final)
write.xlsx(file = 'GSE13836.dzy.xls',vsmc.df)
#write.table(file = 'GSE13836.dzy.txt',vsmc.exprs)

"
GSE23303
"
setwd("D:\\wangli_data\\GSE23303")

gse    <- getGEO(filename = 'GSE23303_series_matrix.txt.gz')
cutoff <- 0.58
str(gse)
all.exprs <- exprs(gse)
boxplot(all.exprs) # this step is to check whether data is normorlized
vsmc.exprs <- all.exprs[,c('GSM571577','GSM571580','GSM571583',
                           'GSM571579','GSM571582','GSM571585')]
design           <- model.matrix(~ 0 + factor(c(1,1,1,2,2,2)))
colnames(design) <- c("vsmc", "control")
fit              <- lmFit(vsmc.exprs, design)
contrast.matrix  <- makeContrasts(vsmc - control, levels = design)
fit2             <- contrasts.fit(fit, contrast.matrix)
fit2             <- eBayes(fit2)
symbols          <- featureData(gse)@data$ORF

result.nofilter.limma <- topTable( fit2, coef = 1, adjust = "BH",
                                   sort.by="B",genelist = symbols,
                                   number = Inf)
write.xlsx(file = 'GSE23303.dzy.xls',result.nofilter.limma)

"
GSE72696 = SRP063675
mv SRR2378601/*.srr ./
mv SRR2378604/*.sra ./
mv SRR2378603/*.sra ./
mv SRR2378606/*.sra ./
fastq-dump.2.4.5 *.sra
ls *.fastq > targets.txt
mkdir results
"
setwd('/home/zhenyisong/data/cardiodata/SRP063675')
setwd("D:\\wangli_data\\GSE72696")
genome_ref.path   <- "/home/zhenyisong/data/bringback/igenome/igenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
targets.file      <- '/home/zhenyisong/data/cardiodata/SRP063675/targets.txt'
reads.files       <- read.table(targets.file,header = F)
reads.path        <- '/home/zhenyisong/data/cardiodata/SRP063675/'
output.path       <- '/home/zhenyisong/data/cardiodata/SRP063675/results/'
reads.paths       <- paste0(reads.path,reads.files$V1)
outputs.files     <- paste0(output.path,reads.files$V1,'.sam')
base.string       <- 'hg19_index'
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

gene              <-  featureCounts( outputs.files, useMetaFeatures = TRUE, 
                                     annot.inbuilt = "hg19", allowMultiOverlap = TRUE)
gene.counts       <- gene$counts
gene.ids          <- gene$annotation$GeneID



#------------------------------------------------------------
#save.image(file = 'atherosclerosis.Rdata')
# setwd("C:\\Users\\Yisong\\Desktop")
# load("atherosclerosis.Rdata")
#------------------------------------------------------------

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

colnames(dge.tmm.counts)
group           <- factor(c( 'Control','PDGF-BB','Control','PDGF-BB'))
d               <- dge.tmm.counts
rownames(d)     <- gene.exprs$genes$SYMBOL
design          <- model.matrix(~ 0 + group);
colnames(design)<- c('Control','Treatment')
contrast.matrix <- makeContrasts(Treatment - Control, levels = design)
d.norm          <- voom(d, design = design)
fit             <- lmFit(d.norm, design)
fit2            <- contrasts.fit(fit,contrast.matrix)
fit2            <- eBayes(fit2)
gene.result     <- topTable(  fit2, 
                              number        = Inf, 
                              adjust.method = "BH", 
                              genelist      = symbols,
                              sort.by       = "p")
write.xlsx(file = 'GSE72696.dzy.xls',gene.result)