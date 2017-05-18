#!/bin/bash
set -e
set -u 
set -o pipefail



# the path to this script
# /home/zhenyisong/biodata/wanglab/wangcode
# nohup xargs -n 1 curl -O -C - -L < files.txt &

# get the NCBI assembly from ftp://mirbase.org/pub/mirbase/CURRENT/genomes/mmu.gff3
# download date, 2017-05-09
# please see all public annotation files at /mnt/date/genomelib
# shasum mmu.gff3
# bc41f38751d3d11f25d9643c744d2c7da28c98ac  mmu.gff3

#miRNA_mmu="/mnt/date/genomelib/annotation/mmu.gff3"
#mmu_genome='/mnt/date/igenomes/Mus_musculus/NCBI/GRCm38/Sequence/WholeGenomeFasta/genome.fa'


#
# I used the ENCODE human sample, miRNA to collect
# download the corresponding data from ENCODE
# 
#---

miRNA_hsa="/mnt/date/genomelib/annotation/hsa.gff3"
hsa_genome='/mnt/date/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa'


#
# public data
# from ENCODE
# /home/zhenyisong/biodata/cardiodata/encodeHeart/humanHeart/miRNA
# 

##cd /home/zhenyisong/biodata/cardiodata/encodeHeart/humanHeart/miRNA
##rm -rf miRNA_trimmed_data
##mkdir miRNA_trimmed_data
##
##miRNA_lane_QC=($(find ./ -name "*.fastq" -type f))
##for filename in ${miRNA_lane_QC[@]}; do
##    trim_galore --small_rna --fastqc --length 18 --phred33 --no_report_file --output_dir miRNA_trimmed_data $filename
##done
##
##rm -rf miRNA_bwa  
##mkdir miRNA_bwa 
##cd miRNA_bwa 
##rm -rf * 
##ln -s $hsa_genome ./
##nohup bwa index genome.fa > bwa.log 2>> bwa.log &
##
##
cd /home/zhenyisong/biodata/cardiodata/encodeHeart/humanHeart/miRNA/miRNA_bwa


miRNA_lane=($(find ../miRNA_trimmed_data -name "*trimmed.fq" -type f))
file_count=${#miRNA_lane[@]}

for (( i=0; i<${file_count}; i++ )); do  
    lane=${miRNA_lane[$i]}
    base_lane=`basename "${lane}"`
    base=${base_lane%_trimmed.fq}
    bwa aln -t 15 -n 1 -o 0 -e 0 -l 8 -k 0 genome.fa ${lane} > ${base_lane}.sai
    bwa samse genome.fa ${base_lane}.sai ${lane} > ${base_lane}.sam
    samtools view -Sb -h -@ 15 ${base_lane}.sam > ${base_lane}.bam
    samtools sort -n -@ 10 -m 1G  ${base_lane}.bam ${base}
    python -m HTSeq.scripts.count -s no --stranded=no -r name -f bam -t miRNA -m union -i Name -q ${base}.bam ${miRNA_hsa} > ${base}.txt
done
##
##
##
### public data set
### GSE67885
###---
##
####cd /home/zhenyisong/biodata/cardiodata/SRP057170
####
####rm -rf miRNA_trimmed_data
####mkdir miRNA_trimmed_data
#####
####miRNA_lane_QC=($(find ./ -name "*.fastq" -type f))
####for filename in ${miRNA_lane_QC[@]}; do
####    trim_galore --small_rna --fastqc --length 16 --phred33 --no_report_file --output_dir miRNA_trimmed_data $filename
####done
####
####multiqc miRNA_trimmed_data
####
####rm -rf miRNA_bwa
####mkdir miRNA_bwa
####cd miRNA_bwa
####rm -rf *
####ln -s $mmu_genome ./
####nohup bwa index genome.fa > bwa.log 2>> bwa.log &
##
##cd /home/zhenyisong/biodata/cardiodata/SRP057170/miRNA_bwa
##
##miRNA_lane=($(find ../miRNA_trimmed_data -name "*trimmed.fq" -type f))
##file_count=${#miRNA_lane[@]}
##
##for (( i=0; i<${file_count}; i++ )); do  
##    lane=${miRNA_lane[$i]}
##    base_lane=`basename "${lane}"`
##    base=${base_lane%_trimmed.fq}
##    bwa aln -t 15 -n 1 -o 0 -e 0 -l 8 -k 0 genome.fa ${lane} > ${base_lane}.sai
##    bwa samse genome.fa ${base_lane}.sai ${lane} > ${base_lane}.sam
##    samtools view -Sb -h -@ 15 ${base_lane}.sam > ${base_lane}.bam
##    samtools sort -n -@ 15 -m 1G -f ${base_lane}.bam ${base}.bam
##    python -m HTSeq.scripts.count -s no --stranded=no -r name -f bam -t miRNA -m union -i Name -q ${base}.bam ${miRNA_mmu} > ${base}.txt
##done
