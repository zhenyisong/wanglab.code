#!/bin/bash
set -e
set -u 
set -o pipefail



# get the NCBI assembly from ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3
# download date, 2017-03-24
# please see all public annotation files at /mnt/date/genomelib
# shasum hsa.gff3
# 444bca9047c25e8151e7845c2d50d020bfe550ec  hsa.gff3
miRNA_hsa="/mnt/date/genomelib/annotation/hsa.gff3"
hsa_genome='/mnt/date/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa'
cd /home/zhenyisong/biodata/wanglab/wangdata/yaoyan/rawdata
#rm -rf miRNA_bwa
#mkdir miRNA_bwa

#rm -rf miRNA_trimmed_data
#mkdir miRNA_trimmed_data
#
#miRNA_lane_QC=($(find ./2015-1-28-miRNA -name "*.clean.fq" -type f))
#for filename in ${miRNA_lane_QC[@]:36:36}; do
#    trim_galore --illumina --no_report_file --output_dir miRNA_trimmed_data $filename
#done
#
#miRNA_lane_QC=($(find ./2015-1-28-miRNA -name "*.clean.fq" -type f))
#mkdir miRNA_fastQC;
#for filename in ${miRNA_lane_QC[@]}; do
#    echo $filename
#    fastqc -t 15 -f fastq -o miRNA_fastQC $filename
#done
#
#multiqc miRNA_fastQC

cd miRNA_bwa
rm -rf *
ln -s $hsa_genome ./
bwa index genome.fa > bwa.log 2>> bwa.log

#find ../2015-1-28-miRNA -name "*clean.fq.gz" -type f | xargs gunzip &
# ( $(find /path/to/toplevel/dir -type f) )
miRNA_lane_1=($(find ../miRNA_trimmed_data -name "*_lane1.clean_trimmed.fq" -type f))
miRNA_lane_2=($(find ../miRNA_trimmed_data -name "*_lane2.clean_trimmed.fq" -type f))
file_count=${#miRNA_lane_1[@]}

for (( i=0; i<${file_count}; i++ )); do  
    lane1=${miRNA_lane_1[$i]}
    base_lane1=`basename "${lane1}"`
    lane2=${miRNA_lane_2[$i]}
    base_lane2=`basename "${lane2}"`
    base=${base_lane1%_lane1.clean_trimmed.fq}
    bwa aln -t 15 -n 1 -o 0 -e 0 -l 8 -k 0 genome.fa ${lane1} > ${base_lane1}.sai
    bwa aln -t 15 -n 1 -o 0 -e 0 -l 8 -k 0 genome.fa ${lane2} > ${base_lane2}.sai
    bwa samse genome.fa ${base_lane1}.sai ${lane1} > ${base_lane1}.sam
    bwa samse genome.fa ${base_lane2}.sai ${lane2} > ${base_lane2}.sam
    samtools view -Sb -h -@ 15 ${base_lane1}.sam > ${base_lane1}.bam
    samtools sort -n -@ 15 -m 1G -f ${base_lane1}.bam ${base_lane1}.sorted.bam
    samtools view -Sb -h -@ 15 ${base_lane2}.sam > ${base_lane2}.bam
    samtools sort -n -@ 15 -m 1G -f ${base_lane2}.bam ${base_lane2}.sorted.bam
    samtools merge -n -f $base.bam ${base_lane1}.sorted.bam ${base_lane2}.sorted.bam
    python -m HTSeq.scripts.count -s no --stranded=no -r name -f bam -t miRNA -m union -i Name -q ${base}.bam ${miRNA_hsa} > ${base}.txt
done
