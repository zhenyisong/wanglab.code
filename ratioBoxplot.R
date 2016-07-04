# @author Yisong Zhen
# @since  2016-06-23
# @update 2016-06-30
library(xlsx)
library(gdata)
library(qvalue)

# please see the post here to beatify the boxplot
# output 
# http://stackoverflow.com/questions/28889977/r-boxplot-how-to-customize-the-appearance-of-the-box-and-whisker-plots-e-g-r

# Genes were sorted using their expression values (FPKM) from cardiomyocyte 1 Month old RNA-seq data in decreasing order.
# Those genes were then catetirgoried into four groups, Q1, Q2, Q3 and Q4. Q1 genes have lowest gene expression values and Q4
# have hihgest gene epxression values.
# Cardiac proliferation genes were grouped based on thier expression ratio on two stages: cardiomyocyte 1 Month old vs.
# cardiomyocyte 6 Month old. If the ratio of FPKM value at two stages is above 1.5, then the genes is classified as 
# a proliferation gene in group A. If the ratio is below the 1/1.5, then the gene is an anti-proliferation gene in group B.
# the overlap between A and Q1 is labeled as group A1. The same is kept for B1, which includes the intersect genes from B and Q1.
# Non-parameter tests (wilcox test) were performed to assess the significance of the difference between A1 vs. B1, A2 vs. B2, A3 vs. B3,
# A4 vs. B4 and other comparsions between different groups were also performed. The p value were then adjusted using  
# Benjamini Hochberg method.


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

q.list              = c()
k.list              = c()

# A vs. B
rip.ratio.non.pos       = rip.ratio[!names(rip.ratio) %in% prolif.pos.gene ]
rip.ratio.non.neg       = rip.ratio[!names(rip.ratio) %in% prolif.neg.gene ]

rip.ratio.pos   = as.numeric(as.character(rip.ratio.pos))
rip.ratio.pos   = na.omit(rip.ratio.pos)
rip.ratio.neg   = as.numeric(as.character(rip.ratio.neg))
rip.ratio.neg   = na.omit(rip.ratio.neg)
p.result        = wilcox.test(rip.ratio.pos,rip.ratio.neg)
k.result        = kruskal.test(list(rip.ratio.pos,rip.ratio.neg))
q.list          = c(p.result$p.value)
k.list          = c(k.result$p.value)


# A vs. nonA
rip.ratio.pos     = as.numeric(as.character(rip.ratio.pos))
rip.ratio.pos     = na.omit(rip.ratio.pos)
rip.ratio.non.pos = as.numeric(as.character(rip.ratio.non.pos))
rip.ratio.non.pos = na.omit(rip.ratio.non.pos)
p.result          = wilcox.test(rip.ratio.pos,rip.ratio.non.pos)
k.result          = kruskal.test(list(rip.ratio.pos,rip.ratio.non.pos))
q.list            = c(q.list,p.result$p.value)
k.list            = c(k.list,k.result$p.value)

# B vs. nonB
rip.ratio.neg     = as.numeric(as.character(rip.ratio.neg))
rip.ratio.neg     = na.omit(rip.ratio.neg)
rip.ratio.non.neg = as.numeric(as.character(rip.ratio.non.neg))
rip.ratio.non.neg = na.omit(rip.ratio.non.neg)
p.result          = wilcox.test(rip.ratio.neg,rip.ratio.non.neg)
k.result          = kruskal.test(list(rip.ratio.neg,rip.ratio.non.neg))
q.list            = c(q.list,p.result$p.value)
k.list            = c(k.list,k.result$p.value)


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



# A1 vs. C1

A1                = intersect(Q1,prolif.pos.gene)
C1                = setdiff(Q1,A1)
A1.ratio          = as.numeric(as.character(rip.ratio[A1]))
A1.ratio          = na.omit(A1.ratio)
C1.ratio          = as.numeric(as.character(rip.ratio[C1]))
C1.ratio          = na.omit(A1.ratio)
p.result          = wilcox.test(A1.ratio,C1.ratio)
k.result          = kruskal.test(list(A1.ratio,C1.ratio))
q.list            = c(q.list,p.result$p.value)
k.list            = c(k.list,k.result$p.value)


# A2 vs. C2
A2                = intersect(Q2,prolif.pos.gene)
C2                = setdiff(Q2,A2)
A2.ratio          = as.numeric(as.character(rip.ratio[A2]))
A2.ratio          = na.omit(A2.ratio)
C2.ratio          = as.numeric(as.character(rip.ratio[C2]))
C2.ratio          = na.omit(C2.ratio)
p.result          = wilcox.test(A2.ratio,C2.ratio)
k.result          = kruskal.test(list(A2.ratio,C2.ratio))
q.list            = c(q.list,p.result$p.value)
k.list            = c(k.list,k.result$p.value)


# A3 vs. C3
A3                = intersect(Q3,prolif.pos.gene)
C3                = setdiff(Q3,A3)
A3.ratio          = as.numeric(as.character(rip.ratio[A3]))
A3.ratio          = na.omit(A3.ratio)
C3.ratio          = as.numeric(as.character(rip.ratio[C3]))
C3.ratio          = na.omit(C3.ratio)
p.result          = wilcox.test(A3.ratio,C3.ratio)
k.result          = kruskal.test(list(A3.ratio,C3.ratio))
q.list            = c(q.list,p.result$p.value)
k.list            = c(k.list,k.result$p.value)



# A4 vs. C4
A4                = intersect(Q4,prolif.pos.gene)
C4                = setdiff(Q4,A4)
A4.ratio          = as.numeric(as.character(rip.ratio[A4]))
A4.ratio          = na.omit(A4.ratio)
C4.ratio          = as.numeric(as.character(rip.ratio[C4]))
C4.ratio          = na.omit(C4.ratio)
p.result          = wilcox.test(A4.ratio,C4.ratio)
k.result          = kruskal.test(list(A4.ratio,C4.ratio))
q.list            = c(q.list,p.result$p.value)
k.list            = c(k.list,k.result$p.value)



# B1 vs. D1
B1                = intersect(Q1,prolif.neg.gene)
D1                = setdiff(Q1,B1)
B1.ratio          = as.numeric(as.character(rip.ratio[B1]))
B1.ratio          = na.omit(B1.ratio)
D1.ratio          = as.numeric(as.character(rip.ratio[D1]))
D1.ratio          = na.omit(D1.ratio)
p.result          = wilcox.test(B1.ratio,D1.ratio)
k.result          = kruskal.test(list(B1.ratio,D1.ratio))
q.list            = c(q.list,p.result$p.value)
k.list            = c(k.list,k.result$p.value)


# B2 vs. D2
B2                = intersect(Q2,prolif.neg.gene)
D2                = setdiff(Q2,B2)
B2.ratio          = as.numeric(as.character(rip.ratio[B2]))
B2.ratio          = na.omit(B2.ratio)
D2.ratio          = as.numeric(as.character(rip.ratio[D2]))
D2.ratio          = na.omit(D2.ratio)
p.result          = wilcox.test(B2.ratio,D2.ratio)
k.result          = kruskal.test(list(B2.ratio,D2.ratio))
q.list            = c(q.list,p.result$p.value)
k.list            = c(k.list,k.result$p.value)


# B3vs. D3
B3                = intersect(Q3,prolif.neg.gene)
D3                = setdiff(Q3,B3)
B3.ratio          = as.numeric(as.character(rip.ratio[B3]))
B3.ratio          = na.omit(B3.ratio)
D3.ratio          = as.numeric(as.character(rip.ratio[D3]))
D3.ratio          = na.omit(D3.ratio)
p.result          = wilcox.test(B3.ratio,D3.ratio)
k.result          = kruskal.test(list(B3.ratio,D3.ratio))
q.list            = c(q.list,p.result$p.value)
k.list            = c(k.list,k.result$p.value)



# B4 vs. D4
B4                = intersect(Q4,prolif.neg.gene)
D4                = setdiff(Q4,B4)
B4.ratio          = as.numeric(as.character(rip.ratio[B4]))
B4.ratio          = na.omit(B4.ratio)
D4.ratio          = as.numeric(as.character(rip.ratio[D4]))
D4.ratio          = na.omit(D4.ratio)
p.result          = wilcox.test(B4.ratio,D4.ratio)
k.result          = kruskal.test(list(B4.ratio,D4.ratio))
q.list            = c(q.list,p.result$p.value)
k.list            = c(k.list,k.result$p.value)


# A1 vs. B1
x = A1.ratio
y = B1.ratio
f <- as.factor(c(rep("a",length(x)), rep("b",length(y))) )
d <- c(x,y)
kruskal.test(d ~ f)  
p.result = wilcox.test(A1.ratio,B1.ratio)
k.result = kruskal.test(list(A1.ratio,B1.ratio))
q.list   = c(q.list,p.result$p.value)
k.list   = c(k.list,k.result$p.value)

# A2 vs. B2
x = A2.ratio
y = B2.ratio
f <- as.factor(c(rep("a",length(x)), rep("b",length(y))) )
d <- c(x,y)
kruskal.test(d ~ f) 
p.result = wilcox.test(A2.ratio,B2.ratio)
k.result = kruskal.test(list(A2.ratio,B2.ratio))
q.list   = c(q.list,p.result$p.value)
k.list   = c(k.list,k.result$p.value)


# A3 vs. B3
p.result = wilcox.test(A3.ratio,B3.ratio)
k.result = kruskal.test(list(A3.ratio,B3.ratio))
q.list   = c(q.list,p.result$p.value)
k.list   = c(k.list,k.result$p.value)


# A4 vs. B4

p.result = wilcox.test(A4.ratio,B4.ratio)
k.result = kruskal.test(list(A3.ratio,B3.ratio))
q.list   = c(q.list,p.result$p.value)
k.list   = c(k.list,k.result$p.value)


#--------------------------------------------------------------------
# use BH method to generate adj.p
#--------------------------------------------------------------------


adj.p             = p.adjust(q.list, method = "BH")
adj.k             = p.adjust(k.list, method = "BH")


#--------------------------------------------------------------------
# output all the figures 
#--------------------------------------------------------------------

#par(mfrow = c(2,2))
#boxplot( list(A = rip.ratio.pos, B = rip.ratio.neg),main = 'A vs. B', 
#         boxlwd = 2, boxwex = 0.5, notch = TRUE, col = c("blue","red"),
#         staplelwd = 4, outline = FALSE)
#text(1.5,2,"adj.p < ", adj.p[1])
#
#boxplot( list(A = rip.ratio.pos, nonA = rip.ratio.non.pos),
#         main = 'A vs. non A', outline = FALSE)
#text(1,2,adj.p[2])
#
#boxplot(list(B = rip.ratio.neg, nonB = rip.ratio.non.neg),
#        main = 'B vs. nonB', outline = FALSE)
#text(1,2,adj.p[3])
#
#
#par(mfrow = c(2,2))
#
#boxplot(list(A1 = A1.ratio, C1 = C1.ratio),main = 'A1 vs. C1', outline = FALSE)
#text(1,2,adj.p[4])
#boxplot(list(A2 = A2.ratio, C2 = C2.ratio),main = 'A2 vs. C2', outline = FALSE)
#text(1,2,adj.p[5])
#boxplot(list(A3 = A3.ratio, C3 = C3.ratio),main = 'A3 vs. C3', outline = FALSE)
#text(1,2,adj.p[6])
#boxplot(list(A4 = A4.ratio, C4 = C4.ratio),main = 'A4 vs. C4', outline = FALSE)
#text(1,2,adj.p[7])
#
#par(mfrow = c(2,2))
#boxplot(list(B1 = B1.ratio, D1 = D1.ratio),main = 'B1 vs. D1', outline = FALSE)
#text(1,2,adj.p[8])
#boxplot(list(B2 = B2.ratio, D2 = D2.ratio),main = 'B2 vs. D2', outline = FALSE)
#text(1,2,adj.p[9])
#boxplot(list(B3 = B3.ratio, D3 = D3.ratio),main = 'B3 vs. D3', outline = FALSE)
#text(1,2,adj.p[10])
#boxplot(list(B4 = B4.ratio, D4 = D4.ratio),main = 'B4 vs. D4', outline = FALSE)
#text(1,2,adj.p[11])

#
#
# this is the output the final figure
# I dicard the multiple-test correction
#
#

# multiple statistical tests/Multiple comparisons problem

par(mfrow = c(2,2))
boxplot( list(A1 = A1.ratio, B1 = B1.ratio),
         boxlwd = 2, boxwex = 0.5, 
         notch = FALSE, col = c("blue","red"),
         staplelwd = 3, outline = FALSE,
         main = 'A1 vs. B1',
         ylab = 'Ratio( RIP/input )')
#p.string = format(adj.p[12],width = 9,digits = 3, justify = c("centre"))
p.string = format(q.list[12],width = 6,digits = 3)
pvalue.string = paste0("p=", p.string)
text(1.3,3,pvalue.string)

boxplot( list(A2 = A2.ratio, B2 = B2.ratio),
         boxlwd = 2, boxwex = 0.5, 
         notch = FALSE, col = c("blue","red"),
         staplelwd = 3, outline = FALSE,
         main = 'A2 vs. B2')
#p.string = format(adj.p[13],width = 9,digits = 3, justify = c("centre"))
p.string = format(q.list[13],width = 6,digits = 3)
pvalue.string = paste0("p=", p.string)
text(1.3,4,pvalue.string)

boxplot( list(A3 = A3.ratio, B3 = B3.ratio),
         boxlwd = 2, boxwex = 0.5, 
         notch = FALSE, col = c("blue","red"),
         staplelwd = 3, outline = FALSE,
         main = 'A3 vs. B3',
         ylab = 'Ratio( RIP/input )')
#p.string = format(adj.p[14],width = 9,digits = 3, justify = c("centre"))
p.string = format(q.list[14],width = 9,digits = 3, justify = c("centre"))
pvalue.string = paste0("p=", p.string)
text(1.48,3,pvalue.string)

boxplot( list(A4 = A4.ratio, B4 = B4.ratio),
         boxlwd = 2, boxwex = 0.5, 
         notch = FALSE, col = c("blue","red"),
         staplelwd = 3, outline = FALSE,
         main = 'A4 vs. B4')
#p.string = format(adj.p[15],width = 9,digits = 3, justify = c("centre"))
p.string = format(q.list[15],width = 9,digits = 3, justify = c("centre"))
pvalue.string = paste0("p=", p.string)
text(1.48,2,pvalue.string)

adj.f             = p.adjust(q.list[12:15], method = "BH")
print(adj.f)