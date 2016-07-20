library(multtest)
library(outliers)
library(nortest)
data(golub)

layout( matrix(c(1,0,2)), widths = c(3,0.5,3),
        heights = c(3,0.5,3), respect = TRUE)

gene.index <- unique(row(golub.gnames)[grep("CCND3", golub.gnames)])
golub.fac  <- factor(golub.cl, levels = 0:1, labels = c("ALL","AML"))
boxplot(golub[gene.index,] ~ gol.fac)

qqnorm(golub[gene.index, golub.fac == "ALL"])
qqline(golub[gene.index, golub.fac == "ALL"])

#grubbs.test(golub[gene.index,gol.fac == "ALL"])

shapiro.test(golub[gene.index,golub.fac == "ALL"])

#ad.test(golub[gene.index,golub.fac == "ALL"])

var.test(golub[gene.index,] ~ golub.fac)

t.test(golub[gene.index,] ~ golub.fac, var.equal = TRUE)

#wilcox.test(golub[gene.index,] ~ golub.fac)