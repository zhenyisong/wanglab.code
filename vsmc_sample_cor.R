library(gplots)
library(xlsx)
library(Rsubread)
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(DESeq2)
library(gplots)
library(genefilter)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)

# @parent program
#    mergeAffy.pl
#    mergeRNAseq.pl

#setwd('/home/zhenyisong/data/wanglilab/vsmc_db');
#load(file = 'vsmc.Rdata')

setwd('/home/zhenyisong/data/wanglilab/vsmc_db');

"
these are processed high-through-put data, those 
data needed the Perl script to extract the 
information. the data were inputed by read.table method

RNA-seq data is imported from the Perl script mergeRNAseq.pl
Affy data is imported from Perl script mergeAffy.pl
"
rna.seq.filename = 'final_rna_seq.cos';
rna.raw.filename = 'final_rna_seq.norm'
non.db.filename  = 'final_nons.cos';
affy.db.filename = 'final_affy.cos';
rna.seq.db = read.table( rna.seq.filename, header = TRUE, sep = "\t",
                         row.names = 1)
rna.raw.db = read.table( rna.raw.filename, header = TRUE, sep = "\t",
                         row.names = 1)
non.db     = read.table( non.db.filename , header = TRUE, sep = "\t",
                         row.names = 1)
affy.db    = read.table( affy.db.filename , header = TRUE, sep = "\t",
                         row.names = 1)
"
transform the data to the data.frame
"
rna.matrix     = as.matrix(rna.seq.db)
rna.raw.matrix = as.matrix(rna.raw.db)

lfc.cutoff     = 0.58

#----------------------------------------------------------------------------
# the follwing code is to extract the common dge gene names
#----------------------------------------------------------------------------
rna.log.matrix = log(rna.raw.matrix + 1)
wangli.data.name    = apply(rna.log.matrix[,c('SRR01','SRR02')], 1, median) - 
                      apply(rna.log.matrix[,c('SRR03','SRR04')], 1, median)
wangli.data.name    = names(wangli.data.name)[wangli.data.name > lfc.cutoff | wangli.data.name < -lfc.cutoff]

GSE35664_1.name     = rna.log.matrix[,'SRR407420'] - rna.log.matrix[,'SRR407419']
GSE35664_1.name     = names(GSE35664_1.name)[GSE35664_1.name > lfc.cutoff | GSE35664_1.name < -lfc.cutoff]

GSE35664_2.name     = rna.log.matrix[,'SRR407421'] - rna.log.matrix[,'SRR407419']
GSE35664_2.name     = names(GSE35664_2.name)[GSE35664_2.name > lfc.cutoff | GSE35664_2.name < -lfc.cutoff]

GSE35664_3.name     = rna.log.matrix[,'SRR407422'] - rna.log.matrix[,'SRR407419']
GSE35664_3.name     = names(GSE35664_3.name)[GSE35664_3.name > lfc.cutoff | GSE35664_3.name < -lfc.cutoff]

GSE38056.name       = rna.log.matrix[,'SRR498452'] - rna.log.matrix[,'SRR498451']
GSE38056.name       = names(GSE38056.name )[GSE38056.name  > lfc.cutoff | GSE38056.name < -lfc.cutoff]

GSE44461.name       = apply(rna.log.matrix[,c('SRR748305','SRR748306','SRR748307')], 1, median) - 
                      apply(rna.log.matrix[,c('SRR748302','SRR748303','SRR748304')], 1, median)
GSE44461.name       = names(GSE44461.name )[GSE44461.name  > lfc.cutoff | GSE44461.name < -lfc.cutoff]

GSE51878_1.name     = apply(rna.log.matrix[,c('SRR1020592','SRR1020593','SRR1020594')], 1, median) - 
                      apply(rna.log.matrix[,c('SRR1020595','SRR1020596','SRR1020597')], 1, median)
GSE51878_1.name     = names(GSE51878_1.name )[GSE51878_1.name  > lfc.cutoff | GSE51878_1.name < -lfc.cutoff]

GSE51878_2.name     = apply(rna.log.matrix[,c('SRR1020598','SRR1020598','SRR1020600')], 1, median) - 
                      apply(rna.log.matrix[,c('SRR1020595','SRR1020596','SRR1020597')], 1, median)
GSE51878_2.name     = names(GSE51878_2.name )[GSE51878_2.name  > lfc.cutoff | GSE51878_2.name < -lfc.cutoff]


GSE60642_1.name     = apply(rna.log.matrix[,c('SRR1555818','SRR1555819')], 1, median) - 
                      apply(rna.log.matrix[,c('SRR1555820','SRR1555821')], 1, median)
GSE60642_1.name     = names(GSE60642_1.name )[GSE60642_1.name  > lfc.cutoff | GSE60642_1.name < -lfc.cutoff]


GSE60642_2.name     = apply(rna.log.matrix[,c('SRR1555816','SRR1555817')], 1, median) - 
                      apply(rna.log.matrix[,c('SRR1555820','SRR1555821')], 1, median)
GSE60642_2.name     = names(GSE60642_2.name )[GSE60642_2.name  > lfc.cutoff | GSE60642_2.name < -lfc.cutoff]

GSE60641.name       = apply(rna.log.matrix[,c('SRR1556006','SRR1556007','SRR1556008')], 1, median) - 
                      apply(rna.log.matrix[,c('SRR1556003','SRR1556004','SRR1556005')], 1, median)
GSE60641.name       = names(GSE60641.name)[GSE60641.name  > lfc.cutoff | GSE60641.name < -lfc.cutoff]

GSE65354_1.name     = rna.log.matrix[,'SRR1776527'] - rna.log.matrix[,'SRR1776529']
GSE65354_1.name     = names(GSE65354_1.name)[GSE65354_1.name  > lfc.cutoff | GSE65354_1.name < -lfc.cutoff]

GSE65354_2.name     = rna.log.matrix[,'SRR1776526'] - rna.log.matrix[,'SRR1776528']
GSE65354_2.name     = names(GSE65354_2.name)[GSE65354_2.name  > lfc.cutoff | GSE65354_2.name < -lfc.cutoff]


# union gene names
RNA.seq.common      = union(wangli.data.name, GSE35664_1.name)
RNA.seq.common      = union(RNA.seq.common, GSE35664_2.name)
RNA.seq.common      = union(RNA.seq.common, GSE35664_3.name)
RNA.seq.common      = union(RNA.seq.common, GSE38056.name)
RNA.seq.common      = union(RNA.seq.common, GSE44461.name)
RNA.seq.common      = union(RNA.seq.common, GSE51878_1.name)
RNA.seq.common      = union(RNA.seq.common, GSE51878_2.name)
RNA.seq.common      = union(RNA.seq.common, GSE60642_1.name)
RNA.seq.common      = union(RNA.seq.common, GSE60642_2.name)
RNA.seq.common      = union(RNA.seq.common, GSE60641.name)
RNA.seq.common      = union(RNA.seq.common, GSE60642_1.name)
RNA.seq.common      = union(RNA.seq.common, GSE60642_2.name)
RNA.seq.common      = union(RNA.seq.common, GSE65354_1.name)
RNA.seq.common      = union(RNA.seq.common, GSE65354_2.name)
#--- the end of RNA-seq DGE common names



wangli.data    = apply(rna.matrix[,c('SRR01','SRR02')], 1, median) - apply(rna.matrix[,c('SRR03','SRR04')], 1, median)
	
GSE35664_1     = rna.matrix[,'SRR407420'] - rna.matrix[,'SRR407419']
GSE35664_2     = rna.matrix[,'SRR407421'] - rna.matrix[,'SRR407419']
GSE35664_3     = rna.matrix[,'SRR407422'] - rna.matrix[,'SRR407419']

GSE38056       = rna.matrix[,'SRR498452'] - rna.matrix[,'SRR498451']
GSE44461       = apply(rna.matrix[,c('SRR748305','SRR748306','SRR748307')], 1, median) - 
                 apply(rna.matrix[,c('SRR748302','SRR748303','SRR748304')], 1, median)
GSE51878_1     = apply(rna.matrix[,c('SRR1020592','SRR1020593','SRR1020594')], 1, median) - 
                 apply(rna.matrix[,c('SRR1020595','SRR1020596','SRR1020597')], 1, median)
GSE51878_2     = apply(rna.matrix[,c('SRR1020598','SRR1020598','SRR1020600')], 1, median) - 
                 apply(rna.matrix[,c('SRR1020595','SRR1020596','SRR1020597')], 1, median)
GSE60642_1     = apply(rna.matrix[,c('SRR1555818','SRR1555819')], 1, median) - 
                 apply(rna.matrix[,c('SRR1555820','SRR1555821')], 1, median)
GSE60642_2     = apply(rna.matrix[,c('SRR1555816','SRR1555817')], 1, median) - 
                 apply(rna.matrix[,c('SRR1555820','SRR1555821')], 1, median)

GSE60641       = apply(rna.matrix[,c('SRR1556006','SRR1556007','SRR1556008')], 1, median) - 
                 apply(rna.matrix[,c('SRR1556003','SRR1556004','SRR1556005')], 1, median)

GSE65354_1     = rna.matrix[,'SRR1776527'] - rna.matrix[,'SRR1776529']
GSE65354_2     = rna.matrix[,'SRR1776526'] - rna.matrix[,'SRR1776528']

pseudo_row     = seq(1:length(wangli.data))
rna.kd.matrix  = c()
rna.kd.matrix  = rbind(pseudo_row,wangli.data,GSE35664_1,GSE35664_2,GSE35664_3)
rna.kd.matrix  = rbind(rna.kd.matrix, GSE38056,GSE44461,GSE51878_1,GSE51878_2 )
rna.kd.matrix  = rbind(rna.kd.matrix, GSE60642_1,GSE60642_2, GSE60641, GSE65354_1,GSE65354_2)
rna.kd.matrix  = rna.kd.matrix[-1,]
rna.kd.matrix  = t(rna.kd.matrix)

colnames(rna.kd.matrix) = c('H2AZ KD','AngII-1h','AngII-3h','AngII-24h','Ang 2','TCF21',
                             'SENCR_kd','SENCR_mk','Jag_kn','Jag_hr','JAG','Hsp60','TNFa')

"
transform the data.frame into the log format
this step is deprecated
the raw RNA-seq data have been already processed
by rlog fucntion, stablize the variance
"
#rna.log.matrix = log(rna.matrix + 1)
affy.matrix    = as.matrix(affy.db)
#colnames(affy.matrix)

GSE29955_1     = apply(affy.matrix[,c(4:6)], 1, median) - 
                 apply(affy.matrix[,c(1:3)], 1, median)
GSE29955_2     = apply(affy.matrix[,c(10:12)], 1, median) - 
                 apply(affy.matrix[,c(7:9)], 1, median)
GSE29955_3     = apply(affy.matrix[,c(16:18)], 1, median) - 
                 apply(affy.matrix[,c(13:15)], 1, median)
GSE36487_1     = affy.matrix[,20] - affy.matrix[,19]
GSE36487_2     = affy.matrix[,21] - affy.matrix[,19]
GSE12261_1     = apply(affy.matrix[,c(25:27)], 1, median) - 
                 apply(affy.matrix[,c(22:24)], 1, median)
GSE12261_2     = apply(affy.matrix[,c(31:33)], 1, median) - 
                 apply(affy.matrix[,c(28:30)], 1, median)
GSE17543       = affy.matrix[,35] - affy.matrix[,34]
GSE13791_1     = apply(affy.matrix[,c(54,55)], 1, median) - 
                 apply(affy.matrix[,c(52,53)], 1, median)
GSE13791_2     = apply(affy.matrix[,c(59:61)], 1, median) - 
                 apply(affy.matrix[,c(56,58)], 1, median)
GSE11367       = apply(affy.matrix[,c(63,65,67)], 1, median) - 
                 apply(affy.matrix[,c(62,64,66)], 1, median)
GSE47744_1     = apply(affy.matrix[,c(70,71)], 1, median) - 
                 apply(affy.matrix[,c(68,69)], 1, median)
GSE47744_2     = apply(affy.matrix[,c(74,75)], 1, median) - 
                 apply(affy.matrix[,c(72,73)], 1, median)
GSE47744_3     = apply(affy.matrix[,c(76,77)], 1, median) - 
                 apply(affy.matrix[,c(72,73)], 1, median)
GSE42813_1     = apply(affy.matrix[,c(78,79)], 1, median) - 
                 apply(affy.matrix[,c(80,81)], 1, median)
GSE42813_2     = apply(affy.matrix[,c(80,81)], 1, median) - 
                 apply(affy.matrix[,c(82,83)], 1, median)
GSE60447_1     = apply(affy.matrix[,c(84,85)], 1, median) - 
                 apply(affy.matrix[,c(86:88)], 1, median)
GSE60447_2     = apply(affy.matrix[,c(82,94)], 1, median) - 
                 apply(affy.matrix[,c(89:91)], 1, median)

GSE19441       = apply(affy.matrix[,c(96,98,100)], 1, median) - 
                 apply(affy.matrix[,c(95,97,99)], 1, median)
GSE56819_1     = apply(affy.matrix[,c(104:106)], 1, median) - 
                 apply(affy.matrix[,c(101:103)], 1, median)
GSE56819_2     = apply(affy.matrix[,c(107:109)], 1, median) - 
                 apply(affy.matrix[,c(101:103)], 1, median)

GSE30004       = apply(affy.matrix[,c(110:121)], 1, median) - 
                 apply(affy.matrix[,c(122:133)], 1, median)

GSE17556_1     = apply(affy.matrix[,c(137:139)], 1, median) - 
                 apply(affy.matrix[,c(134:136)], 1, median)
GSE17556_2     = apply(affy.matrix[,c(140:142)], 1, median) - 
                 apply(affy.matrix[,c(134:136)], 1, median)
GSE17556_3     = apply(affy.matrix[,c(146:148)], 1, median) - 
                 apply(affy.matrix[,c(143:145)], 1, median)
GSE17556_4     = apply(affy.matrix[,c(149:151)], 1, median) - 
                 apply(affy.matrix[,c(146:148)], 1, median)
GSE63425_1     = apply(affy.matrix[,c(152,153,159,163)], 1, median) - 
                 apply(affy.matrix[,c(154,155,158,162)], 1, median)
GSE63425_2     = apply(affy.matrix[,c(156,157,160,161)], 1, median) - 
                 apply(affy.matrix[,c(154,155,158,162)], 1, median)

GSE68021_1     = apply(affy.matrix[,c(167:169)], 1, median) - 
                 apply(affy.matrix[,c(164:166)], 1, median)
GSE68021_2     = apply(affy.matrix[,c(170:172)], 1, median) - 
                 apply(affy.matrix[,c(164:166)], 1, median)
GSE68021_3     = apply(affy.matrix[,c(173:175)], 1, median) - 
                 apply(affy.matrix[,c(164:166)], 1, median)
GSE68021_4     = apply(affy.matrix[,c(176:178)], 1, median) - 
                 apply(affy.matrix[,c(164:166)], 1, median)
GSE68021_5     = apply(affy.matrix[,c(179:181)], 1, median) - 
                 apply(affy.matrix[,c(164:166)], 1, median)
GSE68021_6     = apply(affy.matrix[,c(182:184)], 1, median) - 
                 apply(affy.matrix[,c(164:166)], 1, median)
GSE50251       = apply(affy.matrix[,c(185:187)], 1, median) - 
                 apply(affy.matrix[,c(188:190)], 1, median)
GSE66538       = apply(affy.matrix[,c(195:198)], 1, median) - 
                 apply(affy.matrix[,c(191:194)], 1, median)
GSE66280       = apply(affy.matrix[,c(205:210)], 1, median) - 
                 apply(affy.matrix[,c(199:204)], 1, median)
GSE15713       = affy.matrix[,212] - affy.matrix[,211]
GSE66624       = apply(affy.matrix[,c(219:224)], 1, median) - 
                 apply(affy.matrix[,c(213:218)], 1, median)
GSE13594       = apply(affy.matrix[,c(227,228)], 1, median) - 
                 apply(affy.matrix[,c(225,226)], 1, median)
GSE19106       = apply(affy.matrix[,c(231,232)], 1, median) - 
                 apply(affy.matrix[,c(229,230)], 1, median)
GSE21573_1     = apply(affy.matrix[,c(235,236)], 1, median) - 
                 apply(affy.matrix[,c(233,234)], 1, median)
GSE21573_2     = apply(affy.matrix[,c(237,238)], 1, median) - 
                 apply(affy.matrix[,c(233,234)], 1, median)
GSE21573_3     = apply(affy.matrix[,c(239,240)], 1, median) - 
                 apply(affy.matrix[,c(233,234)], 1, median)
GSE21573_4     = apply(affy.matrix[,c(241,242)], 1, median) - 
                 apply(affy.matrix[,c(233,234)], 1, median)
GSE31080_1     = affy.matrix[,244] - affy.matrix[,243]
GSE31080_2     = affy.matrix[,246] - affy.matrix[,245]
GSE15841_1     = apply(affy.matrix[,c(253,254)], 1, median) - 
                 apply(affy.matrix[,c(247,248)], 1, median)
GSE15841_2     = apply(affy.matrix[,c(255,256)], 1, median) - 
                 apply(affy.matrix[,c(249,250)], 1, median)
GSE15841_3     = apply(affy.matrix[,c(257,258)], 1, median) - 
                 apply(affy.matrix[,c(251,252)], 1, median)
GSE19909       = apply(affy.matrix[,c(259,261)], 1, median) - 
                 apply(affy.matrix[,c(262,264)], 1, median)
GSE21403       = apply(affy.matrix[,c(267,268)], 1, median) - 
                 apply(affy.matrix[,c(265,266)], 1, median)  


pseudo_row     = seq(1:length(GSE21403))
affy.kd.matrix = c()
affy.kd.matrix = rbind(pseudo_row, GSE29955_1,GSE29955_2,GSE29955_3)
affy.kd.matrix = rbind(affy.kd.matrix, GSE36487_1,GSE36487_2,GSE12261_1,GSE12261_2)
affy.kd.matrix = rbind(affy.kd.matrix, GSE17543,GSE13791_1,GSE13791_2,GSE11367)
affy.kd.matrix = rbind(affy.kd.matrix, GSE47744_1,GSE47744_2,GSE47744_3,GSE42813_1,GSE42813_2)
affy.kd.matrix = rbind(affy.kd.matrix, GSE60447_1,GSE60447_2,GSE19441,GSE56819_1,GSE56819_2)
affy.kd.matrix = rbind(affy.kd.matrix, GSE30004,GSE17556_1,GSE17556_2,GSE17556_3,GSE17556_4)
affy.kd.matrix = rbind(affy.kd.matrix, GSE63425_1,GSE63425_2,GSE68021_1,GSE68021_2,GSE68021_3)
affy.kd.matrix = rbind(affy.kd.matrix, GSE68021_4,GSE68021_5,GSE68021_6,GSE50251,GSE66538,GSE66280)
affy.kd.matrix = rbind(affy.kd.matrix, GSE15713,GSE66624,GSE13594,GSE19106)
affy.kd.matrix = rbind(affy.kd.matrix, GSE21573_1,GSE21573_2,GSE21573_3,GSE21573_4)
affy.kd.matrix = rbind(affy.kd.matrix, GSE31080_1,GSE31080_2,GSE15841_1,GSE15841_2,GSE15841_3,GSE19909,GSE21403)


affy.kd.matrix  = affy.kd.matrix[-1,]
affy.kd.matrix  = t(affy.kd.matrix)


# here I changed according to the request by wang-li ox-LDL_24 to ox_LDL
colnames(affy.kd.matrix) = c('OPG',   'RANKL',        'TRAIL',      'moxLDL_3h','moxLDL_21h','ME-treated_4h','ME-treated_30h',
                            'FoxM1', 'T.cruzi_24h', 'T.cruzi_48','IL-17','Cholesterol.1','Cholesterol.1_II','HDL',
                            'APOE+VE', 'APOE','Zyxin','Zyxin_Stretch','CASMC', 'rock1', 'ZIPT','TGF-B_Rx','HG_TSP',
                            'HG1-LG1','Man+TSP1','Man','TAB+GCA','TAB-GCA','LDL_1h','LDL_5h','LDL_24h','ox-LDL_1h','ox-LDL_5h',
                            'ox-LDL', 'embryonic origin-specific','jasplakinolide','glucose','Glut1','Versican',
                             'CD9','PDGF-BB','BMPR2','cdBMPR2','knBMPR2','edBMPR2','IL-1b','PDGF-DD','thrombin','TRAP',
                             'TRAP+PTX','fluid stress','IL-1b')

# extract common gene names from 

lcf.cutoff = lfc.cutoff

GSE29955_1.name     = names(GSE29955_1)[GSE29955_1 > lcf.cutoff | GSE29955_1 < -lcf.cutoff]
GSE29955_2.name     = names(GSE29955_2)[GSE29955_2 > lcf.cutoff | GSE29955_2 < -lcf.cutoff]
GSE29955_3.name     = names(GSE29955_3)[GSE29955_3 > lcf.cutoff | GSE29955_3 < -lcf.cutoff]                 

GSE36487_1.name     = names(GSE36487_1)[GSE36487_1 > lcf.cutoff | GSE36487_1 < -lcf.cutoff] 
GSE36487_2.name     = names(GSE36487_2)[GSE36487_2 > lcf.cutoff | GSE36487_2 < -lcf.cutoff] 
GSE12261_1.name     = names(GSE12261_1)[GSE12261_1 > lcf.cutoff | GSE12261_1 < -lcf.cutoff] 
GSE12261_2.name     = names(GSE12261_2)[GSE12261_2 > lcf.cutoff | GSE12261_2 < -lcf.cutoff] 
GSE17543.name       = names(GSE17543)[GSE17543 > lcf.cutoff | GSE17543 < -lcf.cutoff] 
GSE13791_1.name     = names(GSE13791_1)[GSE13791_1 > lcf.cutoff | GSE13791_1 < -lcf.cutoff] 
GSE13791_2.name     = names(GSE13791_2)[GSE13791_2 > lcf.cutoff | GSE13791_2 < -lcf.cutoff] 

GSE11367.name       = names(GSE11367)[GSE11367 > lcf.cutoff | GSE11367 < -lcf.cutoff] 
GSE47744_1.name     = names(GSE47744_1)[GSE47744_1 > lcf.cutoff | GSE47744_1 < -lcf.cutoff] 
GSE47744_2.name     = names(GSE47744_2)[GSE47744_2 > lcf.cutoff | GSE47744_2 < -lcf.cutoff]
GSE47744_3.name     = names(GSE47744_3)[GSE47744_3 > lcf.cutoff | GSE47744_3 < -lcf.cutoff]  

GSE42813_1.name     = names(GSE42813_1)[GSE42813_1 > lcf.cutoff | GSE42813_1 < -lcf.cutoff]  
GSE42813_2.name     = names(GSE42813_2)[GSE42813_1 > lcf.cutoff | GSE42813_2 < -lcf.cutoff] 
GSE60447_1.name     = names(GSE60447_1)[GSE60447_1 > lcf.cutoff | GSE60447_1 < -lcf.cutoff] 
GSE60447_2.name     = names(GSE60447_2)[GSE60447_2 > lcf.cutoff | GSE60447_2 < -lcf.cutoff] 

GSE19441.name       = names(GSE19441)[GSE19441 > lcf.cutoff | GSE19441 < -lcf.cutoff] 
GSE56819_1.name     = names(GSE56819_1)[GSE56819_1 > lcf.cutoff | GSE56819_1 < -lcf.cutoff] 
GSE56819_2.name     = names(GSE56819_2)[GSE56819_2 > lcf.cutoff | GSE56819_2 < -lcf.cutoff] 

GSE30004.name       = names(GSE30004)[GSE30004 > lcf.cutoff | GSE30004 < -lcf.cutoff] 

GSE17556_1.name     = names(GSE17556_1)[GSE17556_1 > lcf.cutoff | GSE17556_1 < -lcf.cutoff] 
GSE17556_2.name     = names(GSE17556_2)[GSE17556_2 > lcf.cutoff | GSE17556_2 < -lcf.cutoff] 
GSE17556_3.name     = names(GSE17556_3)[GSE17556_3 > lcf.cutoff | GSE17556_3 < -lcf.cutoff] 
GSE17556_4.name     = names(GSE17556_4)[GSE17556_4 > lcf.cutoff | GSE17556_4 < -lcf.cutoff] 
GSE63425_1.name     = names(GSE63425_1)[GSE63425_1 > lcf.cutoff | GSE63425_1 < -lcf.cutoff] 
GSE63425_2.name     = names(GSE63425_1)[GSE63425_2 > lcf.cutoff | GSE63425_2 < -lcf.cutoff] 

GSE68021_1.name     = names(GSE68021_1)[GSE68021_1 > lcf.cutoff | GSE68021_1 < -lcf.cutoff] 
GSE68021_2.name     = names(GSE68021_2)[GSE68021_2 > lcf.cutoff | GSE68021_2 < -lcf.cutoff] 
GSE68021_3.name     = names(GSE68021_3)[GSE68021_3 > lcf.cutoff | GSE68021_3 < -lcf.cutoff] 
GSE68021_4.name     = names(GSE68021_4)[GSE68021_4 > lcf.cutoff | GSE68021_4 < -lcf.cutoff] 
GSE68021_5.name     = names(GSE68021_5)[GSE68021_5 > lcf.cutoff | GSE68021_5 < -lcf.cutoff] 
GSE68021_6.name     = names(GSE68021_6)[GSE68021_6 > lcf.cutoff | GSE68021_6 < -lcf.cutoff] 
GSE50251.name       = names(GSE50251)[GSE50251 > lcf.cutoff | GSE50251 < -lcf.cutoff] 
GSE66538.name       = names(GSE66538)[GSE66538 > lcf.cutoff | GSE66538 < -lcf.cutoff] 
GSE66280.name       = names(GSE66280)[GSE66538 > lcf.cutoff | GSE66280 < -lcf.cutoff] 
GSE15713.name       = names(GSE15713)[GSE15713 > lcf.cutoff | GSE15713 < -lcf.cutoff] 
GSE66624.name       = names(GSE66624)[GSE66624 > lcf.cutoff | GSE66624 < -lcf.cutoff]
GSE13594.name       = names(GSE13594)[GSE13594 > lcf.cutoff | GSE13594 < -lcf.cutoff]
GSE19106.name       = names(GSE19106)[GSE19106 > lcf.cutoff | GSE19106 < -lcf.cutoff]
GSE21573_1.name     = names(GSE21573_1)[GSE21573_1 > lcf.cutoff | GSE21573_1 < -lcf.cutoff]
GSE21573_2.name     = names(GSE21573_2)[GSE21573_2 > lcf.cutoff | GSE21573_2 < -lcf.cutoff]
GSE21573_3.name     = names(GSE21573_3)[GSE21573_3 > lcf.cutoff | GSE21573_3 < -lcf.cutoff]
GSE21573_4.name     = names(GSE21573_4)[GSE21573_4 > lcf.cutoff | GSE21573_4 < -lcf.cutoff]
GSE31080_1.name     = names(GSE31080_1)[GSE31080_1 > lcf.cutoff | GSE31080_1 < -lcf.cutoff]
GSE31080_2.name     = names(GSE31080_2)[GSE31080_2 > lcf.cutoff | GSE31080_2 < -lcf.cutoff]
GSE15841_1.name     = names(GSE15841_1)[GSE15841_1 > lcf.cutoff | GSE15841_1 < -lcf.cutoff]
GSE15841_2.name     = names(GSE15841_2)[GSE15841_2 > lcf.cutoff | GSE15841_2 < -lcf.cutoff]
GSE15841_3.name     = names(GSE15841_3)[GSE15841_3 > lcf.cutoff | GSE15841_3 < -lcf.cutoff]
GSE19909.name       = names(GSE19909)[GSE19909 > lcf.cutoff | GSE19909 < -lcf.cutoff]
GSE21403.name       = names(GSE21403)[GSE21403 > lcf.cutoff | GSE21403 < -lcf.cutoff]

#
affy.common         = union(GSE29955_1.name,GSE29955_2.name)
affy.common         = union(affy.common,GSE29955_3.name)
affy.common         = union(affy.common,GSE36487_1.name)
affy.common         = union(affy.common,GSE36487_2.name)
affy.common         = union(affy.common,GSE12261_1.name)
affy.common         = union(affy.common,GSE12261_2.name)
affy.common         = union(affy.common,GSE17543.name)
affy.common         = union(affy.common,GSE13791_1.name)
affy.common         = union(affy.common,GSE13791_2.name)
affy.common         = union(affy.common,GSE11367.name)
affy.common         = union(affy.common,GSE47744_1.name)
affy.common         = union(affy.common,GSE47744_2.name)
affy.common         = union(affy.common,GSE47744_3.name)
affy.common         = union(affy.common,GSE42813_1.name)
affy.common         = union(affy.common,GSE42813_2.name)
affy.common         = union(affy.common,GSE60447_1.name)
affy.common         = union(affy.common,GSE60447_2.name)
affy.common         = union(affy.common,GSE19441.name)
affy.common         = union(affy.common,GSE56819_1.name)
affy.common         = union(affy.common,GSE56819_2.name)
affy.common         = union(affy.common,GSE30004.name)
affy.common         = union(affy.common,GSE17556_1.name)
affy.common         = union(affy.common,GSE17556_2.name)
affy.common         = union(affy.common,GSE17556_3.name)
affy.common         = union(affy.common,GSE17556_4.name)
affy.common         = union(affy.common,GSE63425_1.name)
affy.common         = union(affy.common,GSE63425_2.name)
affy.common         = union(affy.common,GSE68021_1.name)
affy.common         = union(affy.common,GSE68021_2.name)
affy.common         = union(affy.common,GSE68021_3.name)
affy.common         = union(affy.common,GSE68021_4.name)
affy.common         = union(affy.common,GSE68021_5.name)
affy.common         = union(affy.common,GSE68021_6.name)
affy.common         = union(affy.common,GSE50251.name)
affy.common         = union(affy.common,GSE66538.name)
affy.common         = union(affy.common,GSE66280.name)
affy.common         = union(affy.common,GSE15713.name)
affy.common         = union(affy.common,GSE66624.name)
affy.common         = union(affy.common,GSE13594.name)
affy.common         = union(affy.common,GSE19106.name)
affy.common         = union(affy.common,GSE21573_1.name)
affy.common         = union(affy.common,GSE21573_2.name)
affy.common         = union(affy.common,GSE21573_3.name)
affy.common         = union(affy.common,GSE21573_4.name)
affy.common         = union(affy.common,GSE31080_1.name)
affy.common         = union(affy.common,GSE31080_2.name)
affy.common         = union(affy.common,GSE15841_1.name)
affy.common         = union(affy.common,GSE15841_2.name)
affy.common         = union(affy.common,GSE15841_3.name)
affy.common         = union(affy.common,GSE19909.name)
affy.common         = union(affy.common,GSE21403.name)
#-- the end for affy common name extraction




"    
the following module is to extract the differential expression
gene name from VSMC data analysis.
see the whole script rsubread_limma.R
get the gene names from DEG analysis
"


"
this is the targets file which indicate the fastq file path and other experiment
inforamtion regarding the sequence
"
targets.file      = '/home/zhenyisong/data/wanglilab/projects/2016-05-26/targets.txt'
reads.files       = read.table(targets.file,header = F)

"
output path, where the resuls are saved
"
reads.path        = '/home/zhenyisong/data/wanglilab/projects/2016-05-26/'
output.path       = '/home/zhenyisong/data/wanglilab/projects/2016-05-26/results/'


"
generate the path vectors
"
reads.paths       = paste0(reads.path,reads.files$V1)
outputs.files     = paste0(output.path,reads.files$V1,'.sam')


# get gene's counts
gene = featureCounts( outputs.files, useMetaFeatures = TRUE, 
                      annot.inbuilt = "mm10", allowMultiOverlap = TRUE)

gene.counts  = gene$counts
gene.ids     = gene$annotation$GeneID
colnames(gene.counts) = c( 'VSMC_H2AZ44_Day_4_1','VSMC_H2AZ44_Day_4_2',
                           'VSMC_NT_Day_4_1','VSMC_NT_Day_4_1');


keytypes(org.Mm.eg.db)

columns  = c("ENTREZID","SYMBOL", "MGI", "GENENAME");
GeneInfo = select( org.Mm.eg.db, keys= as.character(gene.ids), 
                   keytype="ENTREZID", columns = columns);
m        = match(gene$annotation$GeneID, GeneInfo$ENTREZID);
Ann      = cbind( gene$annotation[, c("GeneID", "Chr","Length")],
                          GeneInfo[m, c("SYMBOL", "MGI", "GENENAME")]);

rownames(gene.counts) = GeneInfo[m,'SYMBOL'];
#write.table( gene.counts, file = "vsmc.counts.txt", quote = FALSE, 
#             sep = "\t", row.names = TRUE, col.names = TRUE);

Ann$Chr  =  unlist( lapply(strsplit(Ann$Chr, ";"), 
                    function(x) paste(unique(x), collapse = "|")))
Ann$Chr  = gsub("chr", "", Ann$Chr)

gene.exprs          = DGEList(counts = gene.counts, genes = Ann)
gene.exprs.nofilter = calcNormFactors(gene.exprs)
A          = rowSums(gene.exprs$counts)
isexpr     = A > 50

hasannot   = rowSums(is.na(gene.exprs$genes)) == 0
gene.exprs = gene.exprs[isexpr & hasannot, , keep.lib.size = FALSE]
gene.exprs = calcNormFactors(gene.exprs)

d          = gene.exprs

group  = factor(c('CT','CT','TR','TR'));
design = model.matrix(~ 0 + group);
colnames(design) = c('Control','Treatment')
contrast.matrix  = makeContrasts(Treatment - Control, levels = design)
d.norm          = voom(d, design = design)
fit             = lmFit(d.norm, design)
fit2            = contrasts.fit(fit,contrast.matrix)
fit2            = eBayes(fit2)


#gene.result = topTable( fit2, coef = ncol(design), 
#                        number = Inf, adjust.method="BH", sort.by="p");
gene.result = topTable( fit2, number  = Inf, 
                        adjust.method = "BH", 
                        sort.by = "p",
                        lfc     = 0.58,
                        p.value = 0.05);
"
get the DEG gene names and transform the 
gene name to upper case.
"
deg.names   = rownames(gene.result)
deg.names   = toupper(deg.names)
#
# -- module end
#

"
get the common name
and gene expression matrix
"
common.name          = intersect(rownames(rna.kd.matrix), rownames(affy.kd.matrix))
common.name          = intersect(common.name, affy.common)
common.name          = intersect(common.name, RNA.seq.common)
colnames.vector      = c(colnames(rna.kd.matrix),colnames(affy.kd.matrix))
common.matrix        = seq(1:length(colnames.vector))

for( gene in common.name) {
    gene.exprs    = matrix( c( rna.kd.matrix[gene,],
                               affy.kd.matrix[gene,]),
                            byrow = F, nrow = 1)
    common.matrix = rbind(common.matrix,gene.exprs)
}
length(common.name)
exprs.matrix            = common.matrix[-1,]
colnames(exprs.matrix)  = colnames.vector

"
sample name and list is from Wang lis preference
"
sample.names  <- c( 'H2AZ KD', 'PDGF-BB', 'PDGF-DD','ox-LDL', 
                 'APOE+VE', 'JAG', 'fluid stress', 'Ang 2', 
                 'TNFa', 'rock1', 'OPG', 'CD9', 'glucose')
final.data    <- exprs.matrix[,sample.names]

results.final <- cor(final.data,method = 'spearman')

#rna.cor = cor(rna.log.matrix, method = 'spearman')


#----------------------------------------------
# random shuffle the matrix and generate
# the Spearman p.value distribution
# simulation begin
#----------------------------------------------
m_len      = dim(final.data)[2]
pseudo_num = seq(m_len *m_len)
random.matrix = final.data
shuffle.times = 1000
for( j in 1:shuffle.times) {
    for ( i in 1:dim(final.data)[1]) {
        random.matrix[1,] = final.data[1,sample(dim(final.data)[2])]
    }
    pseudo.cor = cor(random.matrix, method = 'spearman')
    pseudo_num = cbind(pseudo_num,as.vector(pseudo.cor))
}
pseudo_num = pseudo_num[,-1]

pdf("pseudo.pvalue.pdf")
hist(as.vector(pseudo_num))
dev.off()
cutoff   = mean(as.vector(pseudo_num) < -0.13)
fileConn = file("cutoff.txt")
writeLines(c("hello","world",cutoff), fileConn)
close(fileConn)



# ---simulation end

#-------------------------------------------------------------
# Figure
# heatmap of sample correlation
#-------------------------------------------------------------


my_palette     = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
#heatmap.result = heatmap.2(results, col = my_palette, scale  = 'none', 
#						   Rowv = T,Colv = T, density.info = 'none',
#                           key  = TRUE, trace='none', symm = T,symkey = F,symbreaks = T,
#						   margins = c(5.5,5.5),dendrogram = 'none',
#                           cexRow  = 0.4, cexCol = 0.4,
#                           labRow = rownames(results),
#                           labCol = colnames(results)
#						  );
#order.sample.names <- rownames(results)[heatmap.result$rowInd]
#partial.name       <- c('VSMC.wang', 'PDGF-BB','PDGF-DD','AngII','TRAP+PTX',
#                         'APOE+VE', 'JAG','oLDL_24h','Jag_kn', 'TNFA-a',
#                         'TGF-B_Rx', 'TAB-GCA', 'LDL_1h', 'LDL_5h',
#                          'fluid Stress')
#
#order.sample.names <- setdiff(order.sample.names, partial.name)
#blackCol           <- rep('black',length(order.sample.names))
#brownCol           <- rep('green4',length(partial.name))
#order.sample.names <- union(partial.name, order.sample.names)
#colRow.vector      <- c(brownCol,blackCol)
#colCol.vector      <- colRow.vector
#
#heatmap.2(results[order.sample.names, order.sample.names], col = my_palette, scale  = 'none', 
#						   Rowv = F,Colv = F, density.info = 'none',
#                           key  = TRUE, trace='none', symm = T,symkey = F,symbreaks = T,
#						   margins = c(5.5,5.5),dendrogram = 'none',
#                           cexRow  = 0.4, cexCol = 0.4,
#                           labRow = order.sample.names,
#                           labCol = order.sample.names,
#                           colRow = colRow.vector,
#                           colCol = colCol.vector
#						  );
#


blackCol           <- rep('black',length(results.final) - 1)
brownCol           <- rep('green4',1)
colRow.vector      <- c(brownCol,blackCol)
colCol.vector      <- colRow.vector

heatmap.2(results.final, col = my_palette, scale  = 'none', 
						   Rowv = F,Colv = F, density.info = 'none',
                           key  = TRUE, trace='none', symm = T,symkey = F,symbreaks = T,
						   margins = c(5.5,5.5),dendrogram = 'none',
                           cexRow  = 0.8, cexCol = 0.8, srtCol = 30,
                           labRow = sample.names,
                           labCol = sample.names,
                           colRow = colRow.vector,
                           colCol = colCol.vector
						  );
plot.new()
results.map = signif(results.final, digits = 2)
grid.table( results.map,
            theme = ttheme_default(base_size = 3.5))

##my_palette  = colorRampPalette(c("white","green","green4","violet","purple"))(255)

#whole.heatmap = heatmap( results,  margins = c(10, 10),
#                         cexCol = 0.4, cexRow = 0.4);
#partial.map = results[results['wangli.data',] > 0.9,results['wangli.data',] > 0.9]
#
#partial.matrix = results[ results['VSMC.wang',] > 0.1 | results['VSMC.wang',] < -0.1, 
#                       results[,'VSMC.wang'] > 0.1 | results[,'VSMC.wang'] < -0.1]
#partial.map    = signif(partial.matrix, digits = 3)
#write.xlsx(partial.map, file = "cutoff_0.8.xls")
#heatmap(  partial.map,  margins = c(10, 10),
#           symm = TRUE, 
#           cexCol = 1, cexRow = 1);

#partial.distance = 1 - partial.map
#heatmap.result = heatmap.2(partial.map, col = my_palette, scale  = 'none', 
#						   Rowv = T,Colv = T, density.info = 'none',
#                           key  = TRUE, trace='none', symm = T,symkey = F,symbreaks = T,
#						   cexRow  = 0.6, cexCol = 0.6,srtCol = 30, cellnote = partial.map,
#                           notecex = 0.5, notecol= "red",
#						   margins = c(10,5),dendrogram = 'none',
#                           labRow = rownames(partial.map),
#                           labCol = colnames(partial.map)
#						  );
#
#rownames(results)[results['SRR01',] > 0.9]
#colnames.vector[whole.heatmap$rowInd]
#summary(as.vector(rna.cor))

#-------------------------------------------------------------
# Spearman cor distribution
#-------------------------------------------------------------

spearman.d = as.vector(pseudo_num)
spearman.d = spearman.d[spearman.d != 1]
hist(spearman.d, prob = TRUE, n = 200, col = 'grey')
lines(density(spearman.d), col = "blue", lwd = 2) # add a density estimate with defaults
lines(density(spearman.d, adjust=2), lty = "dotted", col = "darkgreen", lwd = 2) 


#-------------------------------------------------------------
# Figure
# heatmap of smooth muscle differentiation markers
#-------------------------------------------------------------

vcms.markers      = 'SM-markers.xlsx' # this data is manually curated
vcms.dif.table    = read.xlsx(vcms.markers,header = TRUE, stringsAsFactors = FALSE, sheetIndex = 1)
vcms.sec.table    = read.xlsx(vcms.markers,header = TRUE, stringsAsFactors = FALSE, sheetIndex = 2)
vsmc.dif.genename = vcms.dif.table$GeneSymbol
vsmc.sec.genename = vcms.sec.table$GeneSymbol

dge.tmm                  = t(t(gene.exprs.nofilter$counts) * gene.exprs.nofilter$samples$norm.factors)
#dge.tmm.counts <- round(dge.tmm, digits = 0)
dge.tmm.counts           = apply(dge.tmm,2, as.integer)

# I have no feeling of why this happens
# error message: https://support.bioconductor.org/p/84562/
# rownames(dge.tmm.counts) = as.character(gene.exprs.nofilter$genes$SYMBOL)

sample.info              = data.frame( treat  = c('TR','TR','CT','CT') )
dds                      = DESeqDataSetFromMatrix( countData = dge.tmm.counts,
                                                   colData   = sample.info,
                                                   design    = ~ treat)
vsd                      = varianceStabilizingTransformation(dds, blind = FALSE);
vsd.expr                 = assay(vsd)
colnames(vsd.expr)       = c('VSMC_H2AZ44_1','VSMC_H2AZ44_2','VSMC_NT_1','VSMC_NT_2')
rownames(vsd.expr)       = as.character(gene.exprs.nofilter$genes$SYMBOL)
vsd.dif.markers              = vsd.expr[vsmc.dif.genename,]
vsd.sec.markers              = vsd.expr[vsmc.sec.genename,]

# this code is extracted from 
# http://stackoverflow.com/questions/17820143/how-to-change-heatmap-2-color-range-in-r
# http://seqanswers.com/forums/archive/index.php/t-12022.html

#colors      = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))
#my_palette  = colorRampPalette(c("red", "black", "green"))(n = 299)
#my_palette  = colorRampPalette(c("white","green","green4","violet","purple"))(255)
my_palette  = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
#my_palette  = colorRampPalette(c("red","white","blue"))(256)
heatmap.result = heatmap.2(vsd.dif.markers, col = my_palette, scale  = 'row', 
						   Rowv = TRUE,Colv = FALSE, density.info = 'none',
                           key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
						   cexRow  = 1, cexCol = 1,srtCol = 30,
                           distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						   hclustfun  = function(d) hclust(d, method = 'complete'),
						   dendrogram = 'row',margins = c(12,6),labRow = vsmc.dif.genename,
						  );

heatmap.result = heatmap.2(vsd.sec.markers, col = my_palette, scale  = 'row', 
						   Rowv = TRUE,Colv = FALSE, density.info = 'none',
                           key  = TRUE, trace='none', symm = F,symkey = F,symbreaks = T,
						   cexRow  = 1, cexCol = 1,srtCol = 30,
                           distfun = function(d) as.dist(1-cor(t(d),method = 'pearson')),
						   hclustfun  = function(d) hclust(d, method = 'complete'),
						   dendrogram = 'row',margins = c(12,6),labRow = vsmc.sec.genename,
						  );


#-------------------------------------------------------------
# Figure Gene expression analysis
# Figure Marker assay
#-------------------------------------------------------------

# colnames(dge.tmm)
NT.control.x     <- apply(dge.tmm[,c(3,4)], 1, median)
H2AS.y           <- apply(dge.tmm[,c(1,2)], 1, median)
NT.control.x.log <- log(NT.control.x + 1)
H2AS.y.log       <- log(H2AS.y + 1)

vsmc.exprs.cordinate             <- data.frame(control = NT.control.x.log, treat = H2AS.y.log)

# error message : https://www.biostars.org/p/62988/
# see HTML help page at http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/
#
row.names(vsmc.exprs.cordinate)  <- make.names(gene.exprs.nofilter$genes$SYMBOL, unique = TRUE)
vsmc.exprs.cordinate$type        <- rep(0,dim(vsmc.exprs.cordinate)[1])
vsmc.exprs.cordinate[row.names(vsmc.exprs.cordinate) %in% vsmc.dif.genename,3] <- 1
vsmc.exprs.cordinate[row.names(vsmc.exprs.cordinate) %in% vsmc.sec.genename,3] <- 2
vsmc.dif.df          <- vsmc.exprs.cordinate[vsmc.dif.genename,]
vsmc.sec.df          <- vsmc.exprs.cordinate[vsmc.sec.genename,]
vsmc <- ggplot(vsmc.exprs.cordinate)
vsmc + ylab('Knockdown of H2AZ') +
      geom_point(aes( x = control, y = treat, label = row.names(vsmc.exprs.cordinate),
                      color = as.factor(type), size = as.factor(type), alpha = as.factor(type) ) )+
       scale_size_manual(values = c(2, 2, 2), guide = FALSE) + 
       scale_alpha_manual(values = c(1/200, 1, 1), guide = FALSE) + 
       scale_colour_manual(name = 'gene groups',values = c("black", "blue", "red"), 
                           labels = c('non-related genes','differentiation marker','secretory marker')) + 

       geom_abline(intercept = 0.58, slope = 1, size = 1, alpha = 1/10) +
       geom_abline(intercept = -0.58, slope = 1, size = 1, alpha = 1/10) +
       geom_text( aes(x = control, y = treat), label = rownames(vsmc.dif.df), 
                  data = vsmc.dif.df,hjust = 1,vjust = 1, size = 2, col = 'blue') +
       geom_text( aes(x = control, y = treat), label = rownames(vsmc.sec.df), 
                  data = vsmc.sec.df,hjust = 1,vjust = 1, size = 2, col = 'red') +
       theme(legend.position = c(0.8, 0.2),legend.title.align = 0.5)


#-------------------------------------------------------------
# data is exported and can be replicated in Window system
#
#-------------------------------------------------------------

setwd('/home/zhenyisong/data/wanglilab/vsmc_db');
save.image(file = 'vsmc.Rdata')
quit("no")