## GSM1638233 SRR1922516
## GSM1638234 SRR1922517
setwd('/home/zhenyisong/data/cardiodata/GSE67074')
file.name <- 'GSM1638233_MISO_0006_1_unique_counts_miRBasev20_bwa.txt'
standard.miRNA.count.df <- read.table(file.name, header = TRUE)

## tail -n 2593 SRR1922516.trimed.txt > SRR1922516.trimed.filtered.txt
setwd('/home/zhenyisong/data/cardiodata/SRP056382/bwa')
file.name <- 'SRR1922516.txt'
bwa.miRNA.count.df    <- read.table(file.name, header = FALSE)

setwd('/home/zhenyisong/data/cardiodata/SRP056382/bowtie/adaptor')
file.name <- 'SRR1922516.txt'
bowtie.adaptor.df <- read.table(file.name, header = FALSE)

setwd('/home/zhenyisong/data/cardiodata/SRP056382/miRDeep2/adaptor')
file.name <- 'miRNAs_expressed_all_samples_05_01_2017_t_10_07_56.csv'
miRDeep2.adaptor.df <- read.table(file.name, header = FALSE)

setwd('/home/zhenyisong/data/cardiodata/SRP056382/bowtie/trim')
file.name <- 'SRR1922516.txt'
bowtie.trim.df <- read.table(file.name, header = FALSE)

common.names <- intersect(standard.miRNA.count.df$miRNA, bwa.miRNA.count.df$V1)
common.count <- merge(standard.miRNA.count.df, bwa.miRNA.count.df, by.x = 'miRNA', by.y = 'V1')
common.count <- merge(standard.miRNA.count.df, bwa.miRNA.count.df, by.x = 'miRNA', by.y = 'V1')
common.count <- merge(miDeep2.miRNA.count.df, bowtie.miRNA.count.df, by.x = 'V1', by.y = 'V1')
common.count <- merge(miRDeep2.adaptor.df, bowtie.adaptor.df, by.x = 'V1', by.y = 'V1')
common.count <- merge(miRDeep2.adaptor.df, standard.miRNA.count.df, by.x = 'V1', by.y = 'miRNA')
common.count <- merge(standard.miRNA.count.df, bowtie.trim.df, by.x = 'miRNA', by.y = 'V1')
common.count <- merge(bowtie.trim.df, bowtie.adaptor.df, by.x = 'V1', by.y = 'V1')
common.count <- merge(standard.miRNA.count.df, bowtie.adaptor.df, by.x = 'miRNA', by.y = 'V1')
## head(common.count)
cor(common.count$count,common.count$V2, method = 'spearman')
cor(common.count$count,common.count$V2, method = 'pearson')
cor(log(common.count$count + 0.01),log(common.count$V2 + 0.01), method = 'spearman')
cor(log(common.count$count + 0.01),log(common.count$V2 + 0.01), method = 'spearman')
cor(log(common.count$V2.x + 0.01),log(common.count$V2.y + 0.01), method = 'spearman')
cor(log(common.count$V7 + 1), log(common.count$count + 1), method = 'spearman')
cor(common.count$V7, common.count$count, method = 'spearman')
cor(log(common.count$V7 + 0.01), log(common.count$V2.y + 0.01), method = 'spearman')
cor(common.count$V7, common.count$V2.y, method = 'pearson')
cor(common.count$V2.x, common.count$V2.y, method = 'pearson')
cor(log(common.count$count + 1), log(common.count$V2 + 1), method = 'kendall')