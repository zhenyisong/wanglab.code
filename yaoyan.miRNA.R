#---
# @author Yisong zhen
# @since 2017-03-30
# @updated 2017-10-27
# @parent
#    yaoyan.miRNA.sh
# 
# ageing
# database
# http://senescence.info/
# I downloaded the human ageing related gene
# http://genomics.senescence.info/genes/microarray.php
# unzip and get the 
#---


pkgs <- c( 'tidyverse','Rsubread','org.Mm.eg.db','edgeR',
           'limma', 'DESeq2', 'genefilter',
           'openxlsx','pheatmap','gridExtra','ggrepel',
           'QuasR','annotate','clusterProfiler',
           'GGally','RColorBrewer',
           'cluster','factoextra', 'cowplot',
           'Rsamtools', 'devtools', 'gplots',
           'reshape2', 'stringr','magrittr',
           'multiMiR', 'biomaRt',
           'igraph','org.Hs.eg.db', 'miRNATarget',
           'DiagrammeR','VennDiagram', 'GEOquery',
           'GenomicFeatures', 'rtracklayer',
           'pathview','xlsx'
         )
#--
# deprecated
# install.package(pkgs)
# this will install all necessary packages
# required in this project
#---
#install.lib <- lapply(pkgs, library, character.only = TRUE)
load.lib <- lapply(pkgs, require, character.only = TRUE)

working.env <- 'window'
linux.path  <- file.path('/home/zhenyisong/biodata/wanglab/wangdata/yaoyan/rawdata/miRNA_bwa')
window.path <- file.path('D:\\yisong.data')


#---
#load('yaoyan.encode.Rdata')
#load('yaoyan.Rdata')
#---

switch( working.env, linux  = { setwd(linux.path);
                                load('yaoyan.miRNA.Rdata')},
                     window = { setwd(window.path);
                                load('yaoyan.miRNA.Rdata')} )

figures.path <- file.path('E:\\FuWai\\dead.dir\\wangli.lab\\yaoyan\\yaoyanAgingpaper\\final.figures')
#---
# quantification pipleine digram.
# this is the whole project first step for
# quantifiaction of miR expression level
# pipeline digram, which helps readers to 
# understand the results and inference procedure.
#
# This is the Figure 1A.
#---
graph <-
    create_graph() %>%
    set_graph_name('miRNA quantification pipeline') %>%
    set_global_graph_attrs('graph', 'overlap', 'true') %>%
    set_global_graph_attrs('node', 'color', 'white') %>%
    set_global_graph_attrs('node', 'fontname', 'Arial') %>%
    add_n_nodes(5) %>%
    select_nodes_by_id(c(1:5)) %>% 
    set_node_attrs_ws('shape', 'retangle') %>%
    set_node_attrs_ws('style', 'filled')   %>%
    set_node_attrs_ws('color', 'white')   %>%
    clear_selection %>%
    add_edges_w_string(
      '1->2 2->3 3->4 4->5', 'black') %>%
    set_node_attrs('label',c( 'Reads mapping by the BWA algorithm',
                              'Preprocessing and quality control',
                              'miR reads counting by the HT-Count',
                              'Validation by the public dataset',
                              'Differential expression analysis using the edgeR')) %>%
    set_node_attrs('fontsize',10) %>%
    set_edge_attrs('arrowsize', 1)

render_graph(graph)

#---
# the following code is to test
# the multiplot arrangement
# but failed to save all plots in a figure
#---

"
sampleFile      <- tempfile(pattern = 'zhen3temp', tmpdir = tempdir(), fileext = '.png')
#newsampleFile   <- tempfile(pattern = 'zhen4temp', tmpdir = tempdir(), fileext = '.svg')
pipeline.svg <- export_graph(graph, file_name = sampleFile, file_type = 'png')
#convertPicture(sampleFile, newsampleFile)
#pipeline.img   <-  readPicture(newsampleFile)
pipeline.img   <-  readPNG(sampleFile)

pipeline.ggplot <- rasterGrob(pipeline.img, interpolate = TRUE)
ger.pic <- qplot() + theme_minimal() + annotation_custom(ger)

all.test <- c( coremiR.8.heatmap$gtable, kegg.pathway.plot, 
                  Go.CC.ggplot, ggnet2(miR.targets.go.CC.network))
a <- gtable_add_grob(coremiR.8.heatmap$gtable, rect, 1, 1)
grob.temp <- arrangeGrob( a, kegg.pathway.plot, 
             Go.CC.ggplot, ggnet2(miR.targets.go.CC.network), ncol = 2)
grid.arrange(grobs = grob.temp)
plot2by2 <- plot_grid( coremiR.8.heatmap$gtable, kegg.pathway.plot, 
                       Go.CC.ggplot, ggnet2(miR.targets.go.CC.network),
                       labels=c('A', 'B', 'C', 'D'), ncol = 2)
figure3.ggplot <- plot_grid(aging.miR23.lm.figure3B, plotlist = venns.ggplot2.figure3A, labels=c('A', 'B'), nrow = 2)
grid.arrange(arrangeGrob(aging.miR23.lm.figure3B), venns.ggplot2.figure3A)

"

# Agilen image raw data normlization method
#
# library(LVSmiRNA)

# raw data source from /141219_SN7001347_0224_AHB61HADXX/
# I did not retain the 141217 dir rawdata and did not merge
# them to find the all . discard 1:36
#---


miRNA.count.files  <- list.files(pattern = 'Sample.*\\.txt$')
miRNA.count.matrix <- NULL
for( filename in miRNA.count.files) {
    column <- read.delim(filename, header = F, row.names = 1)
    miRNA.count.matrix <- cbind(miRNA.count.matrix , column$V2 )
}

rownames(miRNA.count.matrix) <- rownames(column)
colnames(miRNA.count.matrix) <- sub('\\.txt','',miRNA.count.files)

groups <- factor(c(rep(1,3),rep(2,3), rep(1,3),rep(2,3)), levels = 1:2, labels = c('Control','Patient'))
ages   <- factor(c(rep(1,6),rep(2,6)), levels = 1:2, labels = c('AF60','AF70'))
AF.miRNA.trans.data <- DGEList(counts = miRNA.count.matrix[,7:18], group = groups ) %>%
                       calcNormFactors(method = 'TMM') %>%
                       estimateCommonDisp() %>%
                       estimateTagwiseDisp(trend = 'movingave')

#
# depracated exactTest for simple analysis
#---
#AF.miRNA.edgeR               <- exactTest(AF.miRNA.trans.data)
#AF.miRNA.edgeR.pvalue        <- AF.miRNA.edgeR$table$PValue
#AF.miRNA.edgeR.p.adj         <- p.adjust(AF.miRNA.edgeR.pvalue, method = 'BH')

# introduce the ages parameter to cancel the batch effect?
#
design               <- model.matrix(~ 0 + groups + ages);
colnames(design)     <- c('Control','Patient','Age')
contrast.matrix      <- makeContrasts(Patient - Control, Age - Control,levels = design)
fit                  <- glmFit(AF.miRNA.trans.data, design)
model.result         <- glmLRT(fit, contrast = contrast.matrix[,1])
AF.miRNA.edgeR.pvalue        <- model.result$table$PValue
AF.miRNA.edgeR.p.adj         <- p.adjust(AF.miRNA.edgeR.pvalue, method = 'BH')


#---
#AF.miRNA.edgeR.sig.table     <- model.result$table[which(AF.miRNA.edgeR.p.adj < 0.05),]
# this step output zero gene below the threshold
# I instead not correct the p.value using BH method
#---
AF.miRNA.edgeR.sig.table           <- model.result$table[which(AF.miRNA.edgeR.pvalue < 0.05),] %>%
                                      arrange(PValue)
rownames(AF.miRNA.edgeR.sig.table) <- rownames(model.result$table)[which(AF.miRNA.edgeR.pvalue < 0.05)]

#---
# deprecated:
#
# require(XLConnect)
#setwd('C:\\Users\\Yisong\\Desktop')
#writeWorksheetToFile('AF-miNRA.xlsx', data = AF.miRNA.edgeR.sig.table, 
#                                 sheet = 'AF-exludedAge')
# thrown out error
# Error: NoSuchMethodError (Java): 
# org.apache.poi.ss.usermodel.Cell.setCellType(Lorg/apache/poi/ss/usermodel/CellType;)V
#---

write.xlsx(x = AF.miRNA.edgeR.sig.table, file = 'AF-miNRA.xlsx', sheet = 'AF-exludedAge')

# now check the age variable influence on 
# samples; we above model using age as the second param to isolate the AF patient, SR patient
# this is a double check.
#---
model.age.result             <- glmLRT(fit, contrast = contrast.matrix[,2])
Aging.miRNA.edgeR.pvalue     <- model.age.result$table$PValue
Aging.miRNA.edgeR.p.adj      <- p.adjust(Aging.miRNA.edgeR.pvalue, method = 'BH')
model.age.result$table$p.adj <- Aging.miRNA.edgeR.p.adj
Aging.miRNA.edgeR.sig.table  <- model.age.result$table[which(Aging.miRNA.edgeR.p.adj < 0.05),] %>%
                                      arrange(PValue)
rownames(Aging.miRNA.edgeR.sig.table) <- rownames(model.age.result$table)[which(Aging.miRNA.edgeR.p.adj < 0.05)]
# check postive control
Aging.miRNA.edgeR.sig.table['hsa-miR-34a-5p',]


#--- AF project end


#---
# branched project
# using the current model to check the miRNA involved in 
# aging process
#---


# follwoing code will generate the figure smiliar to
# NATURE paper regarding the distribution of 
# aging specific miR-34a
# 
sample.SR.info          <- data.frame( treat  = c( 'SR40','SR40','SR40','SR50','SR50', 'SR50',
                                                   'SR60','SR60','SR60','SR70','SR70', 'SR70') )
dge.tmm.miRNA.log       <- miRNA.count.matrix[,c(1:9,13:15)]  %>%
                           DESeqDataSetFromMatrix(colData = sample.SR.info, design = ~ treat) %>%
                           DESeq() %>%
                           counts(normalized = TRUE) %>% add(1) %>% log

rownames(dge.tmm.miRNA.log) <- rownames(column)


# validation or post-QC to check the miR-34a expression

age60.40.epxrs.log    <- dge.tmm.miRNA.log[,c(7:9,1:3)] 

aging.whole.df        <- data.frame( avg   = apply(age60.40.epxrs.log, 1, median), 
                                     ratio = apply(age60.40.epxrs.log[,1:3], 1, median) / apply(age60.40.epxrs.log[,4:6], 1, median))
aging.miR.set         <- aging.whole.df['hsa-miR-34a-5p',]
# remove illagal letter, which can call the function in package
# library(IDPmisc)
# NaRV.omit(df)
# is much better?
#---

#
# this code is from  yaoyan.R
# please re-check the code
#---
"
human.gsea.figure1B    <- gseaplot(aging.1.result, 'aging.pathway')
"
aging.whole.df <- aging.whole.df[complete.cases(aging.whole.df) & !is.infinite(rowSums(aging.whole.df)),]

miR34.valid.figure1C <- ggplot(data = aging.whole.df, aes( x= avg, y = ratio)) +
                        geom_jitter(alpha = 0.5) +
                        geom_point( data = aging.miR.set, aes(x =  avg , y = ratio), 
                                    color = 'red', size = 3.2, shape = 16) +
                        geom_text( data = aging.miR.set, aes(x =  avg , y = ratio), 
                                   label = 'hsa-miR-34a', color = 'blue', 
                                   vjust = -2, hjust = 1) +
                        xlab('average expression counts (log)') + 
                        ylab('ratio old(60)/young(40) (log)') +
                        theme_classic()
#-- figure code end



#---
# post-QC: to check the tissue source
# as this miRNA data set is from human right artium
# using heart specific marker to distinguish the 
# human heart tissue and other tissues
# please refer the reference
# 1: Chistiakov DA, Orekhov AN, Bobryshev YV. Cardiac-specific miRNA in
# cardiogenesis, heart function, and cardiac pathology (with focus on myocardial
# infarction). J Mol Cell Cardiol. 2016 May;94:107-21. doi:
# 10.1016/j.yjmcc.2016.03.015. Epub 2016 Apr 4. Review. PubMed PMID: 27056419.
# I curated the heart specific miRNA from this paper.
#---

# this will generate the figure for cardiac specifc 
# miRs distribution in SR (yaoyan samples)
# cardiac specific miRs were curated from 
#
# PMID: 27056419
# Title:
# Cardiac-specific miRNA in cardiogenesis, 
# heart function, and cardiac pathology 
# (with focus on myocardial infarction)
# author:
# Chistiakov DA1, Orekhov AN2, Bobryshev YV3.
#---

heart.specific.miRNA       <-c( 'miR-499$','miR-133a-5p','miR-208a',
                                'miR-1-3p','miR-1-5p','miR-208b','miR-133a-3p')
sample.SR.info             <- data.frame( treat  = c('SR40','SR40','SR40','SR50','SR50', 'SR50',
                                                     'SR60','SR60','SR60','SR70','SR70', 'SR70') )
SR.samples.exprs.miRNA.log <- miRNA.count.matrix[,c(1:9,13:15)]  %>%
                              DESeqDataSetFromMatrix(colData = sample.SR.info, design = ~ treat) %>%
                              DESeq() %>% counts(normalized = TRUE) %>% add(1) %>% log %>% 
                              apply(1,median) %>% as.data.frame()
                             

# check the matrix manipulation is correct
# colSums(SR.samples.exprs.miRNA)

rownames(SR.samples.exprs.miRNA.log) <- rownames(column)
colnames(SR.samples.exprs.miRNA.log) <- c('SR.median.exprs.log')
SR.samples.exprs.miRNA.filter.log           <- filter(SR.samples.exprs.miRNA.log, SR.median.exprs.log != 0)
rownames(SR.samples.exprs.miRNA.filter.log) <- rownames(column)[SR.samples.exprs.miRNA.log$SR.median.exprs.log != 0]
match.pattern                        <- paste(heart.specific.miRNA,collapse = '|')
heart.ind                            <- grep(match.pattern,rownames(SR.samples.exprs.miRNA.filter.log))
(rownames(SR.samples.exprs.miRNA.filter.log)[heart.ind])

heart.miRNA.exprs                    <- data.frame(exprs = SR.samples.exprs.miRNA.filter.log[heart.ind,1])
heart.miRNA.names                    <- rownames(SR.samples.exprs.miRNA.filter.log)[heart.ind]

heart.specific.jitter  <- jitter(rep(1, dim(heart.miRNA.exprs)[1]), amount = 0.4)
ggplot(data = SR.samples.exprs.miRNA.filter.log, aes(x = 1, y = SR.median.exprs.log)) + 
         geom_point( position = 'jitter' ) + 
         # this label method is deprecated.
         #geom_label_repel( data = heart.miRNA.exprs, mapping = aes(x = 1, y = exprs, label = heart.miRNA.names ))
         geom_point( data = heart.miRNA.exprs, 
                    aes( x = heart.specific.jitter, y = exprs), 
                    color = 'green', size = 2.4) +
         geom_text( data = heart.miRNA.exprs, 
                    aes(x = heart.specific.jitter, y = exprs, label = heart.miRNA.names),
                    color = 'blue', vjust = .5, hjust = -.2) 

#-- tissue specfici QC end




sample.info              <- data.frame( treat  = c('AF40','AF40','AF40','AF50','AF50', 'AF50',
                                                   'AF60','AF60','AF60','AF70','AF70', 'AF70') )
dds.miRNA                <- DESeqDataSetFromMatrix( countData = dge.tmm.miRNA,
                                                    colData   = sample.info,
                                                    design    = ~ treat)
vsd.exprs.miRNA          <- varianceStabilizingTransformation(dds.miRNA, blind = FALSE) %>%
                            assay()

colnames(vsd.exprs.miRNA)   <- c('AF40-1','AF40-2','AF40-3','AF50-1','AF50-2', 'AF50-3',
                                 'AF60-1','AF60-2','AF60-3','AF70-1','AF70-2', 'AF70-3') 

aging.miRNA.sd              <- apply(vsd.exprs.miRNA, 1, sd)
hist(aging.miRNA.sd)
estimator.sd                <- shorth(aging.miRNA.sd)
aging.miRNA.filtered        <- vsd.exprs.miRNA[aging.miRNA.sd > 0.3,]

heatmap.result <- heatmap.2(  aging.miRNA.filtered  , col = greenred(75), scale  = 'row', 
						      Rowv = TRUE,Colv = FALSE, density.info = 'none',
                              key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
                              cexCol = 0.5,srtCol = 30,
                              distfun = function(d) as.dist(1 - cor(t(d), method = 'spearman')),
						      hclustfun  = function(d) hclust(d, method = 'complete'),
						      dendrogram = 'no', labRow = NA);

#aging.dist.miRNA            <- cor(t(aging.miRNA.filtered), method = 'spearman') %>% as.dist()
aging.hclust.miRNA          <- hclust(d = aging.dist.miRNA, method = 'complete')


# updated since 2017-06-29
# now every figure generated depended on 
# the single raw data -- matrix counts
#---


# generate the dge.tmm.miRNA.log for 
# later usage. The new count matrix
# were generated from 
#---
sample.SR.info          <- data.frame( treat  = c('SR40','SR40','SR40','SR50','SR50', 'SR50',
                                                   'SR60','SR60','SR60','SR70','SR70', 'SR70') )
dge.tmm.miRNA.log       <- miRNA.count.matrix[,c(1:9,13:15)]  %>%
                           DESeqDataSetFromMatrix(colData = sample.SR.info, design = ~ treat) %>%
                           DESeq() %>%
                           counts(normalized = TRUE) %>% add(1) %>% log

rownames(dge.tmm.miRNA.log) <- rownames(column)

# miR aging DGE this is the core code
#---
aging.group          <- factor( c(rep(1:4, each = 3)), 
                                levels = 1:4, 
                                labels = c('SR40','SR50','SR60','SR70'))

miRNA.aging.data     <- DGEList( counts = miRNA.count.matrix[,c(1:9,13:15)], group = aging.group ) %>%
                                 calcNormFactors(method = 'TMM') %>%
                                 estimateCommonDisp() %>%
                                 estimateTagwiseDisp(trend = 'movingave')

design               <- model.matrix(~ 0 + aging.group);
colnames(design)     <- levels(aging.group)
contrast.matrix      <- makeContrasts(SR60 - SR40, SR70 - SR50, levels = design)
fit                  <- glmFit(miRNA.aging.data, design)
model.result.1       <- glmLRT(fit, contrast = contrast.matrix[,1])
model.result.2       <- glmLRT(fit, contrast = contrast.matrix[,2])
ageing.miRNA.edgeR.pvalue                <- model.result.1$table$PValue
ageing.miRNA.edgeR.sig.table.1           <- model.result.1$table[which(ageing.miRNA.edgeR.pvalue  < 0.05),] %>%
                                            arrange(PValue)
rownames(ageing.miRNA.edgeR.sig.table.1) <- rownames(model.result.1$table)[which(ageing.miRNA.edgeR.pvalue  < 0.05)]

ageing.miRNA.edgeR.pvalue                <- model.result.2$table$PValue
ageing.miRNA.edgeR.sig.table.2           <- model.result.2$table[which(ageing.miRNA.edgeR.pvalue  < 0.05),] %>%
                                            arrange(PValue)
rownames(ageing.miRNA.edgeR.sig.table.2) <- rownames(model.result.2$table)[which(ageing.miRNA.edgeR.pvalue  < 0.05)]

aging.names <- intersect(rownames(ageing.miRNA.edgeR.sig.table.1), rownames(ageing.miRNA.edgeR.sig.table.2))

aging.names.plus            <- append(aging.names,'hsa-miR-495')
aging.specific.miRNA        <- dge.tmm.miRNA.log[aging.names,]
extra.miRNA                 <- dge.tmm.miRNA.log['hsa-miR-495-3p',] + dge.tmm.miRNA.log['hsa-miR-495-5p',]
aging.specific.miRNA.plus   <- rbind(aging.specific.miRNA, extra.miRNA)
rownames(aging.specific.miRNA.plus) <- aging.names.plus
colnames(aging.specific.miRNA.plus) <- c('SR40-1','SR40-2','SR40-3','SR50-1','SR50-2', 'SR50-3',
                                         'SR60-1','SR60-2','SR60-3','SR70-1','SR70-2', 'SR70-3') 
# see guo howto draw table in ggplot
# this is the Figure S2.
# in corresponding to yao's suggestion, I
# change the miR-34a-5p to mi
#---

table.rownames    <- rownames(aging.specific.miRNA.plus )
table.rownames[3] <- 'hsa-miR-34a'

figureS2          <- grid.newpage() %>% 
                     {tableGrob( round(aging.specific.miRNA.plus, digits = 2), 
                                 rows = table.rownames)} %>%
                     grid.draw()

setwd(figures.path)

grid.newpage() %>% 
{tableGrob( round(aging.specific.miRNA.plus, digits = 2), 
            rows = table.rownames)} %>%
ggsave( 'FigureS2.jpeg', plot = ., dpi = 600, 
        width = 270, height = 80, units = 'mm')


# deprecated, since 2017-06-29
#---

'
g <- tableGrob(rownames(ageing.miRNA.edgeR.sig.table.1), rows = NULL)
grid.newpage()
grid.draw(g)
'

#---
# heatmap of core8 aging miRNA
# core 8, miR list overlap for show
# Figure 2A need ordering
#---

draw_colnames_45 <- function (coln, gaps, ...) {
    coord <- pheatmap:::find_coordinates( length(coln), gaps)
    x     <- coord$coord - 0.5 * coord$size
    res   <- textGrob( coln, x = x, y = unit(1, 'npc') - unit(3,'bigpts'), 
                       vjust   = 0.5, hjust = 1, rot = 45, gp = gpar(...))
    return(res)
}  

assignInNamespace(x = 'draw_colnames', value = 'draw_colnames_45',
                  ns = asNamespace('pheatmap'))

color.bar <- colorRampPalette(c('midnightblue', 'grey', 'mediumvioletred'))(100)
#---
# in corresponding yao's request
# I change the miR name
#---

aging.matrix <- aging.specific.miRNA.plus
rownames(aging.matrix)[3] <- 'hsa-miR-34a'
coremiR.8.figure2A <- pheatmap( aging.matrix, cluster_rows = TRUE, cluster_cols = FALSE,
                                clustering_distance_rows = 'correlation', color = color.bar,
                                scale = 'row', fontsize_col = 10, cellwidth = 23, 
                                cellheight = 23, cuttree_rows = 2)

#--- end figure 2A

#---
# miRNA and lm() graph
#---

aging.correlation <- data.frame( miRNA.exprs = unlist(aging.specific.miRNA['hsa-miR-23a-5p',]), 
                                 age = factor(rep(1:4, each = 3), levels = 1:4, 
                                 labels = c('Age40','Age50','Age60','Age70'),
                                 ordered = TRUE) )


#
# this code is extracted from 
# http://stackoverflow.com/questions/7549694/adding-regression-line-equation-and-r2-on-graph
#---

lm.equation <- function(df){
    m <- lm(miRNA.exprs ~ age, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*','~~italic(r)^2~'='~r2, 
         list( a = format(coef(m)[1], digits = 2), 
               b = format(coef(m)[2], digits = 2), 
              r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

# regressin line using miR-23a
# mimic
# Figure 3?
#---

aging.miR23.lm.figure3B <- ggplot(data = aging.correlation, aes(x = as.numeric(age), y = miRNA.exprs)) +
                           geom_smooth(method = lm, se = TRUE) +
                           geom_point(aes(color = age), show.legend = F) + 
                           xlim('age.40','age.50', 'age.60','age.70') +
                           xlab('human ages') +
                           ylab('miR-23a expression level(log), y')  +
                           geom_text( aes(x = 3, y = 0.4), 
                                      label = lm.equation(aging.correlation), 
                                      parse = TRUE) +
                           theme_classic()


#
# overlapped miRNA genes AF40.60 ~ AF50.70
# prediction their mRNA target
# construct miRNA network concerned with aging process
#---

data(TBL2_HS)

miRNA.map    <- dimnames(TBL2)
ref.seq.homo <- miRNA.map[[1]][miRNA.map[[2]] %in% aging.names]

data(HS_refseq_to_entrezgene)
miRNA.target.geneID <- unique(id_conv[id_conv$refseq_mrna_ncrna %in% ref.seq.homo,2])

kegg.table     <- enrichKEGG( miRNA.target.geneID, organism = 'human', 
                              pvalueCutoff  = 0.05, 
                              pAdjustMethod = 'BH', 
                              qvalueCutoff  = 0.1)
kegg.result       <- summary(kegg.table)
kegg.qvalue       <- -log(kegg.result$qvalue)
kegg.pathway.name <- kegg.result$Description

miR.targets.go.CC <- enrichGO( gene  = miRNA.target.geneID,
                               OrgDb = org.Hs.eg.db,
                               ont   = 'CC',
                               pAdjustMethod = 'BH',
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05)
summary(miR.targets.go.CC)
Go.CC.ggplot      <- summary(miR.targets.go.CC) %>% 
                     as.data.frame %>% mutate( Golabels  = paste(ID, Description, sep = ' '),
                                               LogQvalue = - log(qvalue)) %>%
                     ggplot(aes(x = reorder(Golabels, LogQvalue), y = LogQvalue)) + 
                     geom_bar( stat = 'identity', width = 0.4, 
                               position = position_dodge(width = 0.1), size = 20) +
                     theme(axis.text.x = element_text(angle = 60,hjust = 1, size = 9),
                           axis.text.y = element_text(hjust = 1, size = 9)) +
                     ylab('-log(qvalue)') + xlab('GO CC term enrichment by clusterProfiler') + 
                     coord_flip()


miR.targets.go.CC.network <- enrichMap( miR.targets.go.CC, vertex.label.cex = 1.2, 
                                        layout = igraph::layout.kamada.kawai)

# deprecated miR targets method, I chose the 
# new database for miR targets, multiMiR
#---

"
data(HS_refseq_to_hgnc_symbol)
# these are Core8 miR targets
miRNA.target.Symbol <- unique(id_conv[id_conv$refseq_mrna_ncrna %in% ref.seq.homo,2])
sinkplot()
miRNA.target.Symbol
sinkplot('plot', col = 'black')

miR.target.DGEOverlap.1 <- intersect(miRNA.target.Symbol,aging.result.1.names$SYMBOL)
miR.target.DGEOverlap.2 <- intersect(miRNA.target.Symbol,aging.result.2.names$SYMBOL)
miR.common.targets      <- intersect(miR.target.DGEOverlap.1, miR.target.DGEOverlap.2)
miR.common.targets      <- miR.target.DGEOverlap.1
aging.rpkm.mRNA[miR.common.targets,]


pheatmap( aging.rpkm.mRNA[miR.common.targets,], cluster_rows = TRUE, cluster_cols = FALSE,
          clustering_distance_rows = 'correlation', color = color.bar,
          scale = 'row', fontsize_col = 10, cellwidth = 23, cellheight = 23,
          labels_col = c( 'SR40-1','SR40-2','SR40-3','SR50-1','SR50-2','SR50-3',
                          'SR60-1','SR60-2','SR60-3','SR70-2','SR70-3','SR70-1'))

"

# new miR targets code
#
#---

# these data were from yaoyan.R, please see
# you have to load the yaoyan.Rdata and 
# get the aging.rpkm.mRNA
# Figure 2B
# co-8miRNA expression genes (mRNAs)
#---

aging.core8.targets <- get.multimir( mirna = rownames(aging.specific.miRNA.plus),
                                     org = 'hsa', table = 'predicted', predicted.cutoff = 1,
                                     predicted.cutoff.type = 'p') %$%
                       unique(predicted$target_symbol) %>% as.character

miR.target.40.60   <- intersect(aging.core8.targets,aging.result.1$SYMBOL) %>% na.omit
color.bar          <- colorRampPalette(c('midnightblue', 'grey', 'mediumvioletred'))(100)
aging.core8.targets.figure2B <- miR.target.40.60 %>% sample(50) %>% 
                                {rownames(aging.rpkm.mRNA) %in% . } %>%
                                aging.rpkm.mRNA[.,] %>%
                                pheatmap(cluster_rows = FALSE, cluster_cols = FALSE,
                                         clustering_distance_rows = 'correlation', color = color.bar,
                                         scale = 'row', fontsize_col = 4.5, fontsize_row = 5,
                                         cellwidth = 8, cellheight = 5,
                                         labels_col = c( 'SR40-1','SR40-2','SR40-3',
                                                         'SR50-1','SR50-2','SR50-3',
                                                         'SR60-1','SR60-2','SR60-3',
                                                         'SR70-2','SR70-3','SR70-1'))
#---
keytypes(org.Hs.eg.db)
miR.kegg.df <- mapIds( org.Hs.eg.db, keys = miR.target.40.60, 
                                   column = 'ENTREZID', keytype = 'SYMBOL', 
                                   multiVals = 'CharacterList') %>% 
                         unlist %>% as.integer %>%
                         enrichKEGG( organism = 'human', 
                                     pAdjustMethod = 'none', 
                                     qvalueCutoff = 1) %>%
                         summary() %$%
                         {data.frame( kegg.pvalue  = -log(pvalue),
                                      kegg.pathway = Description )}

kegg.ggplot.figure2C      <- miR.kegg.df[1:13,] %>% 
                             ggplot( aes( x = reorder(kegg.pathway, kegg.pvalue), 
                                          y = kegg.pvalue)) + 
                             geom_bar( stat = 'identity', width = 0.8, 
                                       position = position_dodge(width = 0.1), size = 10) +
                             theme( text        = element_text(size = 9),
                                    axis.text.x = element_text(angle = 60,hjust = 1, size = 9),
                                    axis.text.y = element_text(hjust = 1, size = 9)) +
                             ylab('-log(pvalue)') + 
                             xlab('KEGG pathway GSEA analysis of miRNAs targets') + 
                             coord_flip()


miR.targets.go.CC <- mapIds( org.Hs.eg.db, keys = miR.target.40.60, 
                                   column = 'ENTREZID', keytype = 'SYMBOL', 
                                   multiVals = 'CharacterList') %>% 
                     unlist %>% as.integer %>% 
                     enrichGO( OrgDb = org.Hs.eg.db,
                               ont   = 'BP',
                               pAdjustMethod = 'none',
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 1)
summary(miR.targets.go.CC)
Go.CC.ggplot.figure2D      <- summary(miR.targets.go.CC) %>% 
                              as.data.frame %>% mutate( Golabels  = paste(ID, Description, sep = ' '),
                                                        LogPvalue = - log(pvalue)) %>%
                              arrange(pvalue) %>% filter(pvalue < 0.005) %>%
                              ggplot(aes(x = reorder(Golabels, LogPvalue), y = LogPvalue)) + 
                              geom_bar( stat = 'identity', width = 0.4, 
                                        position = position_dodge(width = 0.05), size = 4) +
                              theme( text        = element_text(size = 9),
                                     axis.text.x = element_text(angle = 60,hjust = 1, size = 9),
                                     axis.text.y = element_text(hjust = 1, size = 9)) +
                              ylab('-log(pvalue)') + xlab('GO CC term enrichment analysis') + 
                              coord_flip()


miR.network.figure2E      <- enrichMap( miR.targets.go.CC, vertex.label.cex = 0.8, 
                                        layout = igraph::layout.kamada.kawai)

# deprecated code!
#

'
par(mar = c(12,4,1,1), fin = c(4,4))

x = barplot( kegg.qvalue, cex.lab = 0.8,cex.axis= 0.8, width = 0.5,
             main = 'KEGG enrichment anlysis', cex.main = 0.8,
             ylab = '-log(q-value of enrichment)')
text( cex = 0.75, x = x - 0.25, y = -1.25, 
      kegg.pathway.name, 
      xpd = TRUE, srt = 60, pos = 2)
'

# the above is the oeverall aging pathway enrichment
#---

# get the correlation between miRNA core 8 with mRNA

aging.mRNA.sd       <- apply(aging.rpkm.mRNA, 1, sd)
shorth(aging.mRNA.sd)
aging.mRNA.filtered <- aging.rpkm.mRNA[aging.mRNA.sd > 0.6,]
rownames(aging.mRNA.filtered)  <- rownames(aging.rpkm.mRNA)[aging.mRNA.sd > 0.6]
aging.mRNA.filtered.reCol      <- aging.mRNA.filtered[,c(1:9,12,11,10)]

aging.total.data               <- rbind(aging.mRNA.filtered.reCol, aging.specific.miRNA.plus)

cor.all.result <- cor(t(aging.total.data), method = 'spearman')
miRNA.cor.excerpt <- cor.all.result[,(dim(cor.all.result)[1] - 7):dim(cor.all.result)[1]]
summary(as.vector(miRNA.cor.excerpt))
hist(as.vector(miRNA.cor.excerpt), breaks = 100)
spearman.data <- data.frame(cor.data = as.vector(miRNA.cor.excerpt))
figureS3A <- ggplot(data = spearman.data, aes(x = cor.data)) +
             geom_histogram( aes(y  = ..density..), binwidth = .01,
                             colour = 'black', fill = 'white') +
             geom_density(alpha = .2, fill = 'lawngreen') +
             xlab('spearman correlation coefficiency') +
             theme_classic() +
             theme(aspect.ratio = 1)

qqPlot(spearman.data$cor.data, ylab = 'Spearman Correlation Coef')

figureS3B <- ggplot(spearman.data, aes(sample = cor.data)) + 
             stat_qq() +
             theme_classic() +
             theme(aspect.ratio = 1)

figureS3 <- plot_grid( figureS3A, figureS3B, 
                       labels = c('A', 'B'), ncol = 2)
setwd(figures.path)
ggsave( 'FigureS3.jpeg', plot = figureS3, 
        width = 250, height = 180, units = 'mm',dpi = 600)


names.vector <- rownames(miRNA.cor.excerpt)
miR.list <- lapply(as.data.frame(miRNA.cor.excerpt), function(x) names.vector[abs(x) > 0.8])

miR.cor.names <- unique(unlist(miR.list))

geneid.2.table <- toTable(org.Hs.egSYMBOL)

miR.cor.geneID <- unique(geneid.2.table$gene_id[geneid.2.table$symbol %in% miR.cor.names])

miR.cor.kegg.df       <- enrichKEGG( miR.cor.geneID , organism = 'human', 
                                      pAdjustMethod = 'none',    
                                      qvalueCutoff  = 1) %>%
                        summary() %$%
                        {data.frame( kegg.pvalue  = -log(pvalue),
                                     kegg.pathway = Description)} %>%
                        filter(kegg.pvalue > 4)

#---
# in 2017-10-30
# I think this is the new Figure 2D
#---

kegg.pathway.figure2F <- miR.cor.kegg.df[1:12,]  %>% ggplot( aes( x = reorder(kegg.pathway, kegg.pvalue), 
                                                       y = kegg.pvalue)) + 
                         geom_bar( stat = 'identity', width = 0.4, 
                                   position = position_dodge(width = 0.1), size = 20) +
                         theme(axis.text.x = element_text(angle = 60,hjust = 1, size = 8),
                               axis.text.y = element_text(hjust = 1, size = 10)) +
                         ylab('-log(pvalue)') + 
                         xlab('KEGG pathway enrichment analysis of miR co-expressed mRNA') + 
                         coord_flip()
     



# test the igraph result
rownames(aging.specific.miRNA.plus)
dim(cor.all.result)
matrix.names <- names(cor.all.result[,14089])[cor.all.result[,14089] > 0.8]

adj.matrix   <- matrix(0, 220, 220)
adj.matrix[220,] <- 1
colnames(adj.matrix) <- matrix.names
rownames(adj.matrix) <- matrix.names
i.graph <- graph_from_adjacency_matrix( adj.matrix, mode = 'directed',
                                        add.colnames = NULL, 
                                        add.rownames = NA)
#plot(i.graph, layout = layout.kamada.kawai, vertex.size = 3)
plot.igraph( i.graph, layout = layout_with_fr, vertex.size = 5,
      vertex.label.cex = 0.8, edge.arrow.size = 0.2)
# pathway mao
setwd('C:\\Users\\Yisong\\Desktop')
aging.pathway.genelist <- aging.result.1$logFC
names(aging.pathway.genelist) <-  aging.result.1$GeneID
aging.pathway <- pathview(gene.data  = aging.pathway.genelist,
                          pathway.id = 'hsa05202',
                          species    = 'hsa',
                          limit      = list(gene = max(abs(aging.pathway.genelist)), cpd = 1))

aging.pathway.2 <- pathview(gene.data  = aging.pathway.genelist,
                          pathway.id = 'hsa05322',
                          species    = 'hsa',
                          limit      = list(gene = max(abs(aging.pathway.genelist)), cpd = 1))

aging.pathway.3 <- pathview(gene.data  = aging.pathway.genelist,
                          pathway.id = 'hsa05034',
                          species    = 'hsa',
                          limit      = list(gene = max(abs(aging.pathway.genelist)), cpd = 1))

#---
# QC data aging
# GSE43556
#
# Boon RA, Iekushi K, Lechner S, Seeger T et al. 
# MicroRNA-34a regulates cardiac ageing and function. 
# Nature 2013 Mar 7;495(7439):107-10. PMID: 23426265
# this QC data will generate the positve set for
# miRNAs correlated with aging process
# and whole aging miRNA dataset from DGE analysis
# DGE (differential gene expression)
#
#---
setwd('E:\\FuWai\\PaperPublished\\CardiacAgeing\\GSE43556')
gse       <- getGEO(filename = 'GSE43556_family.soft.gz')
gsmlist   <- GSMList(gse)
gpl       <- GPLList(gse)

probesets   <- Table(GPLList(gse)[[2]])$ID
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
                                      {tab <- Table(x)
                                       mymatch <- match(probesets,tab$ID_REF)
                                       return(tab$VALUE[mymatch])
                                     }))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
boxplot(data.matrix)
anno.df     <- dataTable(gpl$GPL7732)@table
gene.ind    <- match(probesets, anno.df$ID)
gene.symbol <- anno.df$'miRNA_ID'[gene.ind]
gene.symbol <- make.names(gene.symbol, unique = TRUE)
rownames(data.matrix) <- NULL
rownames(data.matrix) <- gene.symbol

aging.miRNA.pub   <- na.omit(data.matrix[,-c(1:8)]) %>% log
group             <- factor(rep(1:2, each = 4), levels = 1:2, labels = c('Old','Young'))
design            <- model.matrix(~ 0 + group)
colnames(design)  <- levels(group)
contrast.matrix   <- makeContrasts( Old - Young, levels = design)
fit               <- lmFit(aging.miRNA.pub, design)
fit2              <- contrasts.fit(fit,contrast.matrix)
fit2.miRNA        <- eBayes(fit2)

miRNA.result       <- topTable(  fit2.miRNA, 
                                number        = Inf, 
                                adjust.method = 'BH', 
                                sort.by       = 'p',
                                p             = 0.05);

miRNA.result['mmu.miR.34a',]
#                logFC  AveExpr        t     P.Value  adj.P.Val          B
# mmu.miR.34a 0.3961234 6.922052 4.698403 0.001157251 0.01327598 -0.9671097
dim(miRNA.result)

hsa.miRNA.aging.matchnames   <- sub('mmu\\.','hsa-', rownames(miRNA.result))
hsa.miRNA.aging.matchnames   <- gsub('\\.','-',hsa.miRNA.aging.matchnames)
hsa.miRNA.aging.gsub         <- gsub('(hsa-miR-\\w+)-.*','\\1',hsa.miRNA.aging.matchnames,perl = TRUE)

# next is GSEA enrichment analysis

ageing.miRNA.edgeR.1.exprs                <- model.result.1$table$logFC
names(ageing.miRNA.edgeR.1.exprs)         <- seq( 1:length(ageing.miRNA.edgeR.1.exprs ))
ageing.miRNA.edgeR.1.sorted.exprs         <- sort(ageing.miRNA.edgeR.1.exprs, decreasing = TRUE)
aging.index.id                            <- names(ageing.miRNA.edgeR.1.exprs)[rownames(model.result.1$table) %in% hsa.miRNA.aging.matchnames  ]

aging2gene  <- data.frame( diseaseId = unlist(as.character(rep('aging.miRNA',length(aging.index.id )))), 
                                        geneId    = unlist(as.integer(aging.index.id)), check.names = TRUE)

aging.GSEA         <- GSEA( ageing.miRNA.edgeR.1.sorted.exprs, 
                            TERM2GENE = aging2gene, maxGSSize = 5000, pvalueCutoff = 1)
aging.microarray.plot <- gseaplot(aging.GSEA, 'aging.miRNA')

#---
# venns's digram
# AF60 ~ AF40
# AF70 ~ AF50
# Old ~ Young Mouse from Nature miRNA data set
#
#---
homo.AF60.AF40 <- gsub('(hsa-miR-\\w+)-.*','\\1', rownames(ageing.miRNA.edgeR.sig.table.1), perl = TRUE)
homo.AF70.AF50 <- gsub('(hsa-miR-\\w+)-.*','\\1', rownames(ageing.miRNA.edgeR.sig.table.2), perl = TRUE)
mouse.aging    <- hsa.miRNA.aging.gsub 

area.1           <- length(homo.AF60.AF40)
area.2           <- length(homo.AF70.AF50)
area.3           <- length(mouse.aging)
n.12             <- length( intersect(homo.AF60.AF40,homo.AF70.AF50) )
n.13             <- length( intersect(homo.AF60.AF40,mouse.aging) )
n.23             <- length( intersect(homo.AF70.AF50,mouse.aging) )
n.123            <- length( intersect(intersect(homo.AF70.AF50,mouse.aging), intersect(homo.AF60.AF40, mouse.aging)) )
dev.off()
venns.ggplot2.figure3A    <- grid.newpage() %>% 
                             {draw.triple.venn( area.1, area.2, area.3, n.12, n.23, n.13, n.123,
                                                alpha = rep(0.5, 3), fill = c('cyan','limegreen','navyblue'),
                                                col = rep('gray97',3), lwd = rep(0,3), euler.d = T, scaled = F,
                                                cex = c(1.5,1.5,1.5,1.5,1.5,1.5,1.5),
                                                category = c('homo.SR60.SR40','homo.SR70.SR50','mouse.aging'),
                                                cat.cex = 1.5, cat.dist = c(0.08,0.08,0.08)) }


#---
# please check the Nature figure 1.D
# PMID: 23426265
# QC
# this is the QC to check if the data analysis 
# is corrected and evedenced by miNRA analysi by me
#
#---

exprs.df        <- data.frame( avg   = apply(aging.miRNA.pub, 1, median), 
                               ratio = apply(aging.miRNA.pub[,1:4], 1, median) / apply(aging.miRNA.pub[,5:8], 1, median))
exprs.df['mmu.miR.34a',]
ggplot(data = exprs.df, aes( x= avg, y = ratio)) +
geom_point(alpha = 0.02) +
geom_point(aes(x =  exprs.df['mmu.miR.34a', 1] , y = exprs.df['mmu.miR.34a',2]), color = 'red')


#----
# QC for miRNA pipeline
# RNA
# data was download from ENCODE
# please check yaoyan.miRNA.sh
# this preocessed data was save 
# as yaoyan.encode.Rdata
# for tissue specific selection
# 
#---

setwd('/home/zhenyisong/biodata/cardiodata/encodeHeart/humanHeart/miRNA/miRNA_bwa')

miRNA.count.files.encode  <- list.files(pattern = 'ENC.*\\.txt$')
miRNA.count.matrix.encode <- NULL
column.encode             <- NULL
for( filename in miRNA.count.files.encode) {
    column.encode <- read.delim(filename, header = F, row.names = 1)
    miRNA.count.matrix.encode  <- cbind(miRNA.count.matrix.encode , column.encode$V2 )
}

rownames(miRNA.count.matrix.encode) <- rownames(column.encode)
#colnames(miRNA.count.matrix) <- sub('\\.txt','',miRNA.count.files)
miRNA.count.files.encode
colnames(miRNA.count.matrix.encode) <- c('parietal_lobe','temporal_lobe','spinal_cord',
                                         'tongue','occipital_lobe',
                                         'spinal_cord','diencephalon','heart',
                                         'skin_of_body','skeletal_muscle','temporal_lobe','uterus',
                                         'occipital_lobe','metanephros','parietal_lobe',
                                         'lung','skeletal_muscle','skin_of_body',
                                         'heart','cerebellum', 'frontal_cortex',
                                         'uterus','cerebellum','frontal_cortex',
                                         'lung','thyroid_gland','diencephalon',
                                         'thyroid_gland','metanephros',
                                         'tongue','urinary_bladder','liver')

data.colName                        <- colnames(miRNA.count.matrix.encode)[c(1:5,7:10,12,16,21,26,31,32)]
colnames(miRNA.count.matrix.encode) <- NULL
dge.tmm.miRNA.encode    <- (miRNA.count.matrix.encode[,c(1:5,7:10,12,16,21,26,31,32)]) %>% 
                               DESeqDataSetFromMatrix(colData = data.colName, design = ~ 1 ) %>%
                               DESeq() %>% counts(normalized = TRUE) %>% add(1) %>% log

rownames(miRNA.exprs.log.encode) <- rownames(column.encode)
miRNA.exprs.sd.encode            <- apply(miRNA.exprs.log.encode, 1, sd)
hist(miRNA.exprs.sd.encode)
estimator.sd.encode                <- shorth(miRNA.exprs.sd.encode)
miRNA.exprs.filtered.encode        <- miRNA.exprs.log.encode[miRNA.exprs.sd.encode > 0.7,]

heatmap.result <- heatmap.2( miRNA.exprs.filtered.encode  , col = greenred(75), scale  = 'row', 
						         Rowv = TRUE,Colv = FALSE, density.info = 'none',
                                 key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
                                 cexCol = 0.8,srtCol = 30,
                                 distfun = function(d) as.dist(1 - cor(t(d), method = 'spearman')),
						         hclustfun  = function(d) hclust(d, method = 'complete'),
						         dendrogram = 'no',labRow = NA);


#---
# index of heart is 7
#---
tissue.max <- apply(miRNA.exprs.filtered.encode, 1, which.max)
tissue.max[tissue.max != 7] <- 0
tissue.max[tissue.max == 7] <- 1

#------------------------------------------------------------------------
# PMID: 15388519 
# Bioinformatics. 2005 Mar 1;21(5):650-9. Epub 2004 Sep 23.
# Genome-wide midrange transcription profiles reveal expression 
# level relationships in human tissue specification.
# see also:
# PMID: 26891983
# see also:
# https://www.biostars.org/p/209984/
#------------------------------------------------------------------------
tao.func <- function(x) {
    max.val  <- max(x)
    size     <- length(x)
    norm.val <- x/max.val
    norm.val <- 1 - norm.val
    norm.val <- sum(norm.val)
    norm.val <- norm.val/(size - 1)
    return(norm.val)
}

tao.index    <- apply(miRNA.exprs.filtered.encode, 1, tao.func)
tissue.score <- tao.index * tissue.max


miRNA.heart.names <- rownames(miRNA.exprs.filtered.encode)[tissue.score > 0.8]
miRNA.heart.exprs <- miRNA.exprs.filtered.encode[tissue.score > 0.8,]

summary(tissue.score)

heatmap.result <- heatmap.2( miRNA.heart.exprs , col = greenred(75), scale  = 'row', 
						         Rowv = TRUE,Colv = FALSE, density.info = 'none',
                                 key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
                                 cexCol = 0.8,srtCol = 30,
                                 distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						         hclustfun  = function(d) hclust(d, method = 'complete'),
						         dendrogram = 'no',labRow = miRNA.heart.names );


color.bar <- colorRampPalette(c('green4', 'yellow', 'red'))(10)
#color.bar <- colorRampPalette(c('darkblue', 'darkolivegreen1', 'darkorange'))(10)
set.seed(123)
row.length <- dim(miRNA.heart.exprs)[1]
colors.ind <- rep('grey87',row.length)
names(colors.ind) <- rownames(miRNA.heart.exprs)
colors.ind[c(2,3,15,16)] <- 'darkslateblue'
annotation.df     <- data.frame( GeneClass = factor(rownames(miRNA.heart.exprs)) )
ann_colors        <- list(GeneClass = colors.ind)
#pheatmap( miRNA.heart.exprs, cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = 'correlation', color = color.bar, 
#          scale = 'none', fontsize_col = 15, annotation_legend = FALSE,
#          annotation_row = annotation.df, annotation_colors = ann_colors[1] )


#---
# this script is from 
# http://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps
# rotation of pheatmap x-labels
#---
draw_colnames_45 <- function (coln, gaps, ...) {
    coord <- pheatmap:::find_coordinates( length(coln), gaps)
    x     <- coord$coord - 0.5 * coord$size
    res   <- textGrob(coln, x = x, y = unit(1, 'npc') - unit(3,'bigpts'), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
    return(res)
}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x = 'draw_colnames', value = 'draw_colnames_45',
                  ns = asNamespace('pheatmap'))
figureS1 <- pheatmap( miRNA.heart.exprs, cluster_rows = TRUE, cluster_cols = FALSE,
                      clustering_distance_rows = 'correlation', color = color.bar, 
                      scale = 'none', fontsize_col = 12 ) 

"
setwd(figures.path)

jpeg('FigureS1.jpeg',unit = 'mm', res = 600)
grid.draw(figureS1$gtable)
dev.off()
"

#--- end for tissue QC for pipeline


setwd('/home/zhenyisong/biodata/wanglab/wangdata/Rdata')
save.image('yaoyan.encode.Rdata')
setwd('D:\\wangli_data\\Rdata')
save.image('yaoyan.miRNA.Rdata')
q('no')