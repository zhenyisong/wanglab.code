#!/bin/bash

# the script is to demultiplexing the result from next-500
# @rawdate path
# @samplesheet edited by 
#---
#baseCalls='/mnt/date/Sequencing/RawData/0314/170314_NB501912_0005_AH3WNYBGX2/Data/Intensities/BaseCalls'
#sampleSheet='/mnt/date/Sequencing/FastQ/sampleSheet/2017_03_14_totalRNA.csv'
#dataOutput='/mnt/date/Sequencing/FastQ/totalRNA_2017_03_14'
#runFolder='/mnt/date/Sequencing/RawData/0314/170314_NB501912_0005_AH3WNYBGX2'
#bcl2fastq --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions \
#          --no-lane-splitting -R $runFolder -i $baseCalls -r 10 -d 10 -p 10 \
#          --sample-sheet $sampleSheet -o $dataOutput 2>> /dev/null &
#
# hisat2 to assemble the mRNA 
# generate the GTF file for hg19
# data was downloaded from NONCODE database
# shasum *.gtf
# bfbd03d5d09880628e673bae6d128bafe6ba4b67  NONCODE2016_human_hg38_lncRNA.gtf
#  shasum *hg19.gtf
#  1e5600869854fb1b7a408df841d2c32f4049a669  NONCODE2016.hg19.gtf
# gffread NONCODE2016.hg19.gtf -g $hg19_genome -w hg19_lncRNA.fa


# strandness prediction
#hg19_RefSeq='/bioware/RSeQClib/hg19_RefSeq.bed'           
#infer_experiment.py -i A1-70OM_S2.bam -r $hg19_RefSeq
#Reading reference gene model /bioware/RSeQClib/hg19_RefSeq.bed ... Done
#Loading SAM/BAM file ...  Total 200000 usable reads were sampled
#
#
#This is PairEnd Data
#Fraction of reads failed to determine: 0.0159
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.0169
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.9672

#hg19_mRNA_GTF='/mnt/date/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf'
#hg19_lncRNA_GTF='/mnt/date/genomelib/annotation/NONCODE2016.hg19.gtf'
#hg19_genome='/mnt/date/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
#rawdata='/mnt/date/Sequencing/FastQ/totalRNA_2017_03_14'
#bamdata='/home/zhenyisong/biodata/wanglab/wangdata/lncRNAClinic/hisat2'
#hg19_lncGeneGTF='/mnt/date/genomelib/annotation/human_hg19_GenePusLncRNA.gtf'
#merge_gtf='/home/zhenyisong/biodata/wanglab/wangdata/lncRNAClinic/hisat2/stringtie_merged.gtf'
#mergelist='mergelist.txt'
#cd /mnt/date/genomelib/annotation
#cat $hg19_mRNA_GTF > human_hg19_GenePusLncRNA.gtf
#cat $hg19_lncRNA_GTF >> human_hg19_GenePusLncRNA.gtf
#shasum human_hg19_GenePusLncRNA.gtf
#1f46090b603fbc5a4a4fc4574f5c3aa68de6e3f3  human_hg19_GenePusLncRNA.gtf

#cd hisat2
#hisat2-build -q -f -p 15 $hg19_genome genome &
cd $rawdata
files1=(*R1_001.fastq)
files2=(*R2_001.fastq)
len=${#files1[@]}

cd $bamdata

##for (( i=0; i<$((len - 1)); i++ ));
##do
##    forward=${files1[$i]}
##    backward=${files2[$i]}
##    base=${forward%_R*_001.fastq}
##    hisat2 -p 15 --dta --fr --rna-strandness FR -x genome -1 $rawdata/$forward -2 $rawdata/$backward -S $base.sam
##    samtools view -Sb $base.sam -@ 15 -o $base.bam
##    samtools sort -@ 15 $base.bam $base
##    stringtie -p 15 -e -B -G $hg19_lncGeneGTF -o ballgown/$base/$base.gtf $base.bam  
##done

for (( i=0; i<$((len - 1)); i++ ));
do
    forward=${files1[$i]}
    base=${forward%_R*_001.fastq}
    hisat2 -p 15 --dta --rna-strandness FR -x genome -U $rawdata/$forward -S $base.sam
    samtools view -Sb $base.sam -@ 15 -o $base.bam
    samtools sort -@ 15 $base.bam $base
    stringtie -p 15 -e -B -G $hg19_lncGeneGTF -o ballgown/$base/$base.gtf $base.bam  
done

prepDE.py -i ballgown