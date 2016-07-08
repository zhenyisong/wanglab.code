library(affy)
library(annotate)
library(oligo)

library(rae230a.db)
library(rat2302.db)
library(hgu133a.db)
library(hgu133a2.db)
library(hgu133plus2.db)
library(mouse4302.db)
library(pd.mogene.1.0.st.v1)
library(pd.mogene.2.0.st)
library(pd.ragene.1.0.st.v1)
library(pd.huex.1.0.st.v2)
library(pd.hugene.1.0.st.v1)

# https://bioinformaticsagatha.wordpress.com/2011/03/10/1-expression-analysis-of-affymetrix-arrays-using-bioconductor/

setwd("D:\\wangli_data\\GSE15713")
raw.data   = ReadAffy();
rma.data   = rma(raw.data);
exprs1.data = exprs(rma.data)
probes     = rownames(exprs1.data)
symbols    = unlist(mget(probes, rae230aSYMBOL, ifnotfound = NA))


setwd("D:\\wangli_data\\GSE66624")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs2.data = exprs(rma.data)

# combine
norm.data  = cbind(symbols,exprs1.data,exprs2.data)
setwd("D:\\wangli_data\\read_results")
write.table( norm.data, file = "rae230a.db.txt", quote = FALSE, 
             sep = "\t", row.names = FALSE, col.names = TRUE)

#export rat2302.db data
setwd("D:\\wangli_data\\GSE13594")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs1.data = exprs(rma.data)
probes      = rownames(exprs1.data)
symbols     = unlist(mget(probes, rat2302SYMBOL, ifnotfound = NA))

setwd("D:\\wangli_data\\GSE19106")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs2.data = exprs(rma.data)

setwd("D:\\wangli_data\\GSE21573")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs3.data = exprs(rma.data)


setwd("D:\\wangli_data\\GSE31080")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs4.data = exprs(rma.data)

setwd("D:\\wangli_data\\GSE15841")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs5.data = exprs(rma.data)

setwd("D:\\wangli_data\\GSE19909")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs6.data = exprs(rma.data)


setwd("D:\\wangli_data\\GSE21403")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs7.data = exprs(rma.data)


norm.data   = cbind( symbols,exprs1.data,exprs2.data,
                     exprs3.data, exprs4.data,
                     exprs5.data, exprs6.data,
                     exprs7.data   )



# export hgu133a.db

setwd("D:\\wangli_data\\GSE36487")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs1.data = exprs(rma.data)
probes      = rownames(exprs1.data)
symbols     = unlist(mget(probes, hgu133aSYMBOL, ifnotfound = NA))

setwd("D:\\wangli_data\\GSE4725")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs2.data = exprs(rma.data)



norm.data   = cbind( symbols,exprs1.data,exprs2.data  )

setwd("D:\\wangli_data\\read_results")
write.table( norm.data, file = "hgu133a.db.txt", quote = FALSE, 
             sep = "\t", row.names = FALSE, col.names = TRUE)

# export hgu133a2.db

setwd("D:\\wangli_data\\GSE29955")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs1.data = exprs(rma.data)
probes      = rownames(exprs1.data)
symbols     = unlist(mget(probes, hgu133a2SYMBOL, ifnotfound = NA))

norm.data   = cbind( symbols,exprs1.data)

setwd("D:\\wangli_data\\read_results")
write.table( norm.data, file = "hgu133a2.db.txt", quote = FALSE, 
             sep = "\t", row.names = FALSE, col.names = TRUE)

# export hgu133plus2.db

setwd("D:\\wangli_data\\GSE12261")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs1.data = exprs(rma.data)
probes      = rownames(exprs1.data)
symbols     = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound = NA))

setwd("D:\\wangli_data\\GSE17543")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs2.data = exprs(rma.data)

setwd("D:\\wangli_data\\GSE13791")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs3.data = exprs(rma.data)

setwd("D:\\wangli_data\\GSE11367")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs4.data = exprs(rma.data)

norm.data   = cbind( symbols,exprs1.data,
                     exprs2.data,exprs3.data,exprs4.data)

setwd("D:\\wangli_data\\read_results")
write.table( norm.data, file = "hgu133plus2.db.txt", quote = FALSE, 
             sep = "\t", row.names = FALSE, col.names = TRUE)

# export mouse4302.db

setwd("D:\\wangli_data\\GSE47744")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs1.data = exprs(rma.data)
probes      = rownames(exprs1.data)
symbols     = unlist(mget(probes, mouse4302SYMBOL, ifnotfound = NA))

setwd("D:\\wangli_data\\GSE42813")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs2.data = exprs(rma.data)

setwd("D:\\wangli_data\\GSE60447")
raw.data    = ReadAffy();
rma.data    = rma(raw.data);
exprs3.data = exprs(rma.data)

norm.data   = cbind( symbols,exprs1.data,
                     exprs2.data,exprs3.data)

setwd("D:\\wangli_data\\read_results")
write.table( norm.data, file = "mouse4302.db.txt", quote = FALSE, 
             sep = "\t", row.names = FALSE, col.names = TRUE)

# export pd.mogene.1.0.st.v1
setwd("D:\\wangli_data\\GSE50251")
celPath  = list.files( getwd(), pattern = "\\.CEL|\\.CEL\\.gz", 
                      full.names = TRUE, ignore.case = TRUE)
raw.data = read.celfiles(celPath, pkgname = "pd.mogene.1.0.st.v1")
rma.data = rma(raw.data,target = 'core')
featureData(rma.data) =  getNetAffx(rma.data,'transcript')
gene.symbol           = pData(featureData(rma.data))[, 'geneassignment']
exprs1.data           = exprs(rma.data)

setwd("D:\\wangli_data\\GSE66538")
celPath  = list.files( getwd(), pattern = "\\.CEL|\\.CEL\\.gz", 
                      full.names = TRUE, ignore.case = TRUE)
raw.data = read.celfiles(celPath, pkgname = "pd.mogene.1.0.st.v1")
rma.data = rma(raw.data,target = 'core')
featureData(rma.data) =  getNetAffx(rma.data,'transcript')
exprs2.data           = exprs(rma.data)


norm.data             = cbind( gene.symbol,exprs1.data,exprs2.data)

setwd("D:\\wangli_data\\read_results")
write.table( norm.data, file = "pd.mogene.1.0.st.v1.txt", quote = FALSE, 
             sep = "\t", row.names = FALSE, col.names = TRUE)

# export pd.mogene.2.0.st

setwd("D:\\wangli_data\\GSE66280")
celPath  = list.files( getwd(), pattern = "\\.CEL|\\.CEL\\.gz", 
                      full.names = TRUE, ignore.case = TRUE)
raw.data = read.celfiles(celPath, pkgname = "pd.mogene.2.0.st")
rma.data = rma(raw.data,target = 'core')
featureData(rma.data) =  getNetAffx(rma.data,'transcript')
gene.symbol           = pData(featureData(rma.data))[, 'geneassignment']
exprs1.data           = exprs(rma.data)

norm.data             = cbind( gene.symbol,exprs1.data)

setwd("D:\\wangli_data\\read_results")
write.table( norm.data, file = "pd.mogene.2.0.st.txt", quote = FALSE, 
             sep = "\t", row.names = FALSE, col.names = TRUE)

# export pd.ragene.1.0.st.v1

setwd("D:\\wangli_data\\GSE35627")
celPath  = list.files( getwd(), pattern = "\\.CEL|\\.CEL\\.gz", 
                      full.names = TRUE, ignore.case = TRUE)
raw.data = read.celfiles(celPath, pkgname = "pd.ragene.1.0.st.v1")
rma.data = rma(raw.data,target = 'core')
featureData(rma.data) =  getNetAffx(rma.data,'transcript')
gene.symbol           = pData(featureData(rma.data))[, 'geneassignment']
exprs1.data           = exprs(rma.data)

setwd("D:\\wangli_data\\E-MTAB-1888")
celPath  = list.files( getwd(), pattern = "\\.CEL|\\.CEL\\.gz", 
                      full.names = TRUE, ignore.case = TRUE)
raw.data = read.celfiles(celPath, pkgname = "pd.mogene.1.0.st.v1")
rma.data = rma(raw.data,target = 'core')
featureData(rma.data) =  getNetAffx(rma.data,'transcript')
exprs2.data           = exprs(rma.data)

norm.data             = cbind( gene.symbol,exprs1.data,exprs2.data)

setwd("D:\\wangli_data\\read_results")
write.table( norm.data, file = "pd.ragene.1.0.st.v1.txt", quote = FALSE, 
             sep = "\t", row.names = FALSE, col.names = TRUE)




