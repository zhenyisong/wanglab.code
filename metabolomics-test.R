#install.packages("readMzXmlData")
#install.packages("MALDIquant")
# wget https://regis-web.systemsbiology.net/rawfiles/lcq/7MIX_STD_110802_1.mzXML
# yahoo tiny.pwiz.1.1.mzML
# and download the original data
library(MALDIquant)
library(readMzXmlData)
library(mzR)
setwd('E:\\FuWai\\wangli.lab\\metabolomics')
metabolic.data <- readMzXmlFile('7MIX_STD_110802_1.mzXML')
class(metabolic.data)
length(metabolic.data)
head(summary(metabolic.data))
dim(summary(metabolic.data))

#metabolic.data <- openMSfile('tiny.pwiz.1.1.mzML')

mass.spectrum <- createMassSpectrum( mass = metabolic.data[[1]]$spectrum$mass, 
                                     intensity = metabolic.data[[1]]$spectrum$intensity)

plot(mass.spectrum,col = 'black')

