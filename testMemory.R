library(gplots)
library(xlsx)
library(Rsubread)
library(edgeR)

# cd /home/zhenyisong/wanglab/wangcode
# nohup R CMD BATCH testMemory.R &

genome_ref.path     <- "/bioware/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"

"
GSE49906/SRP029464
PMID:24752171
"

setwd('/home/zhenyisong/cardiodata/SRP029464')
reads.files.names    <- list.files(pattern = '*.fastq')
read.path.1          <- reads.files.names[grep("_1",reads.files.names)]
read.path.2          <- reads.files.names[grep("_2",reads.files.names)]
reads.paths.1        <- paste0(getwd(),'/',read.path.1)
reads.paths.2        <- paste0(getwd(),'/',read.path.2)
unlink('rsubread')
dir.create('rsubread')
output.path          <- '/home/zhenyisong/cardiodata/SRP029464/rsubread'
setwd('/home/zhenyisong/cardiodata/SRP029464/rsubread')
outputs.files        <- paste0(output.path,'/', read.path.1,'.bam')


"
the base index name
"
base.string         <- 'mm10_index'

"
use the Rsubread command to generate index file
this index file will be generated and saved at getwd()
you do not need to generate the script
"

buildindex( basename = base.string, reference = genome_ref.path )

"
this is the function which is called to align the genome
sequence
"

align( index          = base.string, 
       readfile1      = reads.paths.1, 
       readfile2      = reads.paths.2, 
       input_format   = "FASTQ", 
       type           = 'rna',
       output_file    = outputs.files, 
       output_format  = "BAM",
       PE_orientation = 'fr', 
       nthreads       = 12, 
       indels         = 1,
       maxMismatches  = 3,
       phredOffset    = 33,
       unique         = T )

cooper.gene    <-   featureCounts( outputs.files, useMetaFeatures = TRUE,
                                   countMultiMappingReads = FALSE,
                                   strandSpecific         = 1, 
                                   isPairedEnd            = TRUE,
                                   requireBothEndsMapped  = TRUE,
                                   autosort               = TRUE,
                                   nthreads               = 12,
                                   annot.inbuilt = "mm10", allowMultiOverlap = TRUE)
setwd("/home/zhenyisong/cardiodata/SRP029464")
save.image(file = 'testmemory.Rdata')
quit("no")