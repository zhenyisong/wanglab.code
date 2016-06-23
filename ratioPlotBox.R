# @author Yisong Zhen
# @since  2016-06-23
# @update 2016-06-23
library(xlsx)
library(gdata)

# please see the post here to beatify the boxplot
# output 
# http://stackoverflow.com/questions/28889977/r-boxplot-how-to-customize-the-appearance-of-the-box-and-whisker-plots-e-g-r

"
set the working directory; the raw data in their Excel format
was deposited in this dir
"

setwd("/home/zhenyisong/data/wanglilab/projects/2016-06-21")

"
the gene differentially expressed in CM1 AND CM6 two status
"
file.deg = 'Heart_RNAseq_Analysis_20140624.xlsx'
"
the gene with different binding affinity
"
file.rip = '15_09_RIP-Seq_log2_fold_change.xlsx'


# deprecated 
# the reading speed is too slow!
# further information should be noticed
# http://www.milanor.net/blog/read-excel-files-from-r/
# deg.data = read.xlsx(file.deg, 1) # read the second sheet

"
get the data from Excel files
"

deg.data = read.xls(file.deg, sheet = 1, header = TRUE)
rip.data = read.xls(file.rip, sheet = 1, header = TRUE)

rip.ratio           = rip.data$log2_fold_change
names(rip.ratio)    = rip.data$symbol
rip.ratio           = na.omit(rip.ratio)

"
the original fold change was calculated by CM6/CM1
this is the reverse to our selection criteria
I thus reverse the order
"
fold.change.pos     = 1/1.5
fold.change.neg     = 1.5
prolif.pos.gene     = deg.data$gene[deg.data$Fold_change < fold.change.pos]
prolif.pos.gene     = na.omit(prolif.pos.gene)
prolif.neg.gene     = deg.data$gene[deg.data$Fold_change > fold.change.neg]
prolif.neg.gene     = na.omit(prolif.neg.gene)

rip.ratio.pos       = rip.ratio[prolif.pos.gene]
rip.ratio.pos       = as.numeric(as.character(rip.ratio.pos))
rip.ratio.pos       = na.omit(rip.ratio.pos)
rip.ratio.neg       = rip.ratio[prolif.neg.gene]
rip.ratio.neg       = as.numeric(as.character(rip.ratio.neg))
rip.ratio.neg       = na.omit(rip.ratio.neg)


par(mfrow = c(2,2))

rip.ratio.non.pos       = rip.ratio[!names(rip.ratio) %in% prolif.pos.gene ]
rip.ratio.non.neg       = rip.ratio[!names(rip.ratio) %in% prolif.neg.gene ]

rip.ratio.pos   = as.numeric(as.character(rip.ratio.pos))
rip.ratio.neg   = as.numeric(as.character(rip.ratio.neg))
boxplot( list(A = rip.ratio.pos, B = rip.ratio.neg),main = 'A vs. B', 
         boxlwd = 4, notch = TRUE, col = c("blue","red"),
         staplelwd = 4, outline = FALSE)
p.result       = wilcox.test(rip.ratio.pos,rip.ratio.neg)
text(1,2,p.result$p.value)


rip.ratio.pos     = as.numeric(as.character(rip.ratio.pos))
rip.ratio.non.pos = as.numeric(as.character(rip.ratio.non.pos))
boxplot(list(A = rip.ratio.pos, nonA = rip.ratio.non.pos),main = 'A vs. non A', outline = FALSE)
p.result  = wilcox.test(rip.ratio.pos,rip.ratio.non.pos)
text(1,2,p.result$p.value)

rip.ratio.neg     = as.numeric(as.character(rip.ratio.neg))
rip.ratio.non.neg = as.numeric(as.character(rip.ratio.non.neg))
boxplot(list(B = rip.ratio.neg, nonB = rip.ratio.non.neg),main = 'B vs. nonB', outline = FALSE)
p.result = wilcox.test(rip.ratio.neg,rip.ratio.non.neg)
text(1,2,p.result$p.value)

#qqnorm(rip.ratio.pos);qqline(rip.ratio.pos, col = 2)


"
the following step is to calssify the gene expression level
into the 4 according to their expression in CM_1
"
gene.exprs        = deg.data$CM1_FPKM
names(gene.exprs) = deg.data$gene
gene.exprs        = na.omit(gene.exprs)
gene.exprs        = sort(gene.exprs)

size              = 4

gene.group        = split(gene.exprs, 
                          ceiling(seq_along(gene.exprs)/(length(gene.exprs)/size)))

Q1                = names(gene.group$'1') 
Q2                = names(gene.group$'2')
Q3                = names(gene.group$'3')
Q4                = names(gene.group$'4')


par(mfrow = c(2,2))
# A1 vs. C1

A1                = intersect(Q1,prolif.pos.gene)
C1                = setdiff(Q1,A1)
A1.ratio          = as.numeric(as.character(rip.ratio[A1]))
C1.ratio          = as.numeric(as.character(rip.ratio[C1]))
boxplot(list(A1 = A1.ratio, C1 = C1.ratio),main = 'A1 vs. C1', outline = FALSE)
p.result = wilcox.test(A1.ratio,C1.ratio)
text(1,2,p.result$p.value)

# A2 vs. C2
A2                = intersect(Q2,prolif.pos.gene)
C2                = setdiff(Q2,A2)
A2.ratio          = as.numeric(as.character(rip.ratio[A2]))
C2.ratio          = as.numeric(as.character(rip.ratio[C2]))
boxplot(list(A2 = A2.ratio, C2 = C2.ratio),main = 'A2 vs. C2', outline = FALSE)
p.result = wilcox.test(A2.ratio,C2.ratio)
text(1,2,p.result$p.value)

# A3 vs. C3
A3                = intersect(Q3,prolif.pos.gene)
C3                = setdiff(Q3,A3)
A3.ratio          = as.numeric(as.character(rip.ratio[A3]))
C3.ratio          = as.numeric(as.character(rip.ratio[C3]))
boxplot(list(A3 = A3.ratio, C3 = C3.ratio),main = 'A3 vs. C3', outline = FALSE)
p.result = wilcox.test(A3.ratio,C3.ratio)
text(1,2,p.result$p.value)


# A4 vs. C4
A4                = intersect(Q4,prolif.pos.gene)
C4                = setdiff(Q4,A4)
A4.ratio          = as.numeric(as.character(rip.ratio[A4]))
C4.ratio          = as.numeric(as.character(rip.ratio[C4]))
boxplot(list(A4 = A4.ratio, C4 = C4.ratio),main = 'A4 vs. C4', outline = FALSE)
p.result = wilcox.test(A4.ratio,C4.ratio)
text(1,2,p.result$p.value)



par(mfrow = c(2,2))
# B1 vs. D1
B1                = intersect(Q1,prolif.neg.gene)
D1                = setdiff(Q1,B1)
B1.ratio          = as.numeric(as.character(rip.ratio[B1]))
D1.ratio          = as.numeric(as.character(rip.ratio[D1]))
boxplot(list(B1 = B1.ratio, D1 = D1.ratio),main = 'B1 vs. D1', outline = FALSE)
p.result = wilcox.test(B1.ratio,D1.ratio)
text(1,2,p.result$p.value)


# B2 vs. D2
B2                = intersect(Q2,prolif.neg.gene)
D2                = setdiff(Q2,B2)
B2.ratio          = as.numeric(as.character(rip.ratio[B2]))
D2.ratio          = as.numeric(as.character(rip.ratio[D2]))
boxplot(list(B2 = B2.ratio, D2 = D2.ratio),main = 'B2 vs. D2', outline = FALSE)
p.result = wilcox.test(B2.ratio,D2.ratio)
text(1,2,p.result$p.value)

# B3vs. D3
B3                = intersect(Q3,prolif.neg.gene)
D3                = setdiff(Q3,B3)
B3.ratio          = as.numeric(as.character(rip.ratio[B3]))
D3.ratio          = as.numeric(as.character(rip.ratio[D3]))
boxplot(list(B3 = B3.ratio, D3 = D3.ratio),main = 'B3 vs. D3', outline = FALSE)
p.result = wilcox.test(B3.ratio,D3.ratio)
text(1,2,p.result$p.value)


# B4 vs. D4
B4                = intersect(Q4,prolif.neg.gene)
D4                = setdiff(Q4,B4)
B4.ratio          = as.numeric(as.character(rip.ratio[B4]))
D4.ratio          = as.numeric(as.character(rip.ratio[D4]))
boxplot(list(B4 = B4.ratio, D4 = D4.ratio),main = 'B4 vs. D4', outline = FALSE)
p.result = wilcox.test(B4.ratio,D4.ratio)
text(1,2,p.result$p.value)

par(mfrow = c(2,2))
# A1 vs. B1
boxplot(list(A1 = A1.ratio, B1 = B1.ratio),main = 'A1 vs. B1', outline = FALSE)
p.result = wilcox.test(A1.ratio,B1.ratio)
text(1,2,p.result$p.value)

# A2 vs. B2
boxplot(list(A2 = A2.ratio, B2 = B2.ratio),main = 'A2 vs. B2', outline = FALSE)
p.result = wilcox.test(A2.ratio,B2.ratio)
text(1,2,p.result$p.value)

# A3 vs. B3
boxplot(list(A3 = A3.ratio, B3 = B3.ratio),main = 'A3 vs. B3', outline = FALSE)
p.result = wilcox.test(A3.ratio,B3.ratio)
text(1,2,p.result$p.value)

# A4 vs. B4
boxplot(list(A4 = A4.ratio, B4 = B4.ratio),main = 'A4 vs. B4', outline = FALSE)
p.result = wilcox.test(A4.ratio,B4.ratio)
text(1,2,p.result$p.value)



