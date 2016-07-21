# dummy guide for R programming
# since 2016-07-19
# Author Yisong Zhen
#

"
load the library neccessary for demo
you should install the package
"
"
how to install the R package
"
# https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages
# install.packages("ggplot2")
# install.packages("vioplot")
# install.packages("multtest")
# source("http://bioconductor.org/biocLite.R")
# biocLite("multtest")
library(ggplot2)
library(vioplot)
library(multtest)
data(golub)



## Simulate some data

## 3 Factor Variables
FacVar1 = as.factor(rep(c("level1", "level2"), 25))
FacVar2 = as.factor(rep(c("levelA", "levelB", "levelC"), 17)[-51])
FacVar3 = as.factor(rep(c("levelI", "levelII", "levelIII", "levelIV"), 13)[-c(51:52)])

## 4 Numeric Vars
set.seed(123)
NumVar1 = round(rnorm(n = 50, mean = 1000, sd = 50), digits = 2)  ## Normal distribution

set.seed(123)
NumVar2 = round(runif(n = 50, min = 500, max = 1500), digits = 2)  ## Uniform distribution

set.seed(123)
NumVar3 = round(rexp(n = 50, rate = 0.001))  ## Exponential distribution

NumVar4 = 2001:2050

simData = data.frame(FacVar1, FacVar2, FacVar3, NumVar1, NumVar2, NumVar3, NumVar4)

boxplot(NumVar1 ~ FacVar2, data = simData)
plot(simData$FacVar3)  ## bar plot
bartable = table(simData$FacVar2, simData$FacVar3)  ## get the cross tab
barplot(bartable, beside = TRUE, legend = levels(unique(simData$FacVar2)))  ## plot


#ggplot(  data = simData, aes(x = FacVar2, y = NumVar1, fill = FacVar1)) +
#         geom_bar(stat="identity", position=position_dodge())

x <- rnorm(100)
y <- rnorm(100)
plot(x, y, xlim=c(-5,5), ylim=c(-5,5))
vioplot( x, col = "tomato", horizontal = TRUE, 
         at = -4, add = TRUE,lty = 2, rectCol = "gray")
vioplot( y, col = "cyan", horizontal = FALSE, 
         at = -4, add = TRUE,lty = 2)


boxplot( NumVar1 ~ FacVar2, horizontal = F, 
         xlab = "God bless us",outpch = NA)
stripchart( NumVar1 ~ FacVar2, vertical = T, 
            method = "jitter", pch = 21, col = "blue", 
            bg = "yellow", add = TRUE)


## NumVar4 is 2001 through 2050... possibly, 
## a time variable - use that as
## the x-axis
plot( simData$NumVar4, simData$NumVar1, type = "o", 
      ylim = c(0, max(simData$NumVar1, simData$NumVar2)))  ## join dots with lines

lines(simData$NumVar4, simData$NumVar2, type = "o", lty = 2, col = "red")  ## add another line


hist(simData$NumVar3, prob = TRUE, col = 'grey')
lines(density(simData$NumVar3), col = "blue", lwd = 2) 