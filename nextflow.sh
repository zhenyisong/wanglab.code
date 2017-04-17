#!/bin/bash

#set -e
#set -o pipefail
# https://github.com/ewels/NGI-RNAseq
#cd /temp


# Define the constant in the whole QC pipeline setting
# we could add or upgrade the content here later
# outside data and annotation file, path
# test command:
# nextflow.sh -p 'test123' -g mm10 -r "data/*.fastq"
#----

# How to get execution time of a script effectively?
mm10_UCSC_genome='/mnt/date/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa'
mm9_UCSC_genome='/mnt/date/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa'
hg19_UCSC_genome='/mnt/date/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
hg38_UCSC_genome=''
mm10_UCSC_GTF='/mnt/date/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
mm9_UCSC_GTF='/mnt/date/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf'
hg19_UCSC_GTF='/mnt/date/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf'
hg38_UCSC_GTF=''

#---

# define the config
# base + 
#---
config_file='/bioware/lib/wanglab.config'

genome=
project=
strand='none'
GTFfile=

while [[ -n $1 ]];do
    case $1 in 
        -g | --genome)    shift
                          genome=$1
                          ;;
                  
        -p | --project)   shift
                          project=$1
                          ;;
        -r | --reads)     shift
                          reads=$1
                          ;;
        -s | --strand)    shift
                          strand=$1
                          ;;
        esac
        shift
done


case $genome in 
    hg19)     echo "you are now using the UCSC hg19 as the genome assemble."
              genome="$hg19_UCSC_genome"
              GTFfile="$hg19_UCSC_GTF"
              ;;
    hg38)     echo "you are now using the UCSC hg38 as the genome assemble."
              genome="$hg38_UCSC_genome"
              GTFfile="$hg38_UCSC_GTF"
              ;;
    mm9)      echo "you are now using the UCSC mm9 as the genome assemble."
              genome="$mm9_UCSC_genome"
              GTFfile="$mm9_UCSC_GTF"
              ;;
    mm10)     echo "you are now using the UCSC mm10 as the genome assemble."
              genome="$mm10_UCSC_genome"
              GTFfile="$mm10_UCSC_GTF"
              ;;
    *)        echo " Invalid genome assemble version: hg19 hg38 mm9 mm10" >&2
              echo $genome
              exit 1
              ;;
esac


case $strand in 
    FR)     echo 'you use the second strand to construct library, --fr-secondstrand'
            strand='--forward_stranded'
            ;;
    RF)     echo 'you use the first strand to construct library, antisense'
            strand='--reverse_stranded'
            ;;
    none)   echo ' your library is unstranded'
            strand='--unstranded'
            ;;
    *)      echo 'Invalid strandness parameter, please check it, FR/RF/none'
            exit 1
            ;;
esac


random_projectID=`mktemp -u  NextflowQC_XXXXXX`
project=${project:="$random_projectID"}
echo "your project ID is ${project} created for this case, please write down in your worksheet\n"
echo $genome
echo $GTFfile

#nextflow run SciLifeLab/NGI-RNAseq -profile base -c $config_file \
#         --project $project --reads $reads --fasta $genome --gtf $GTFfile

#  work            # Directory containing the nextflow working files
#  results         # Finished results (configurable, see below)
#  .nextflow_log   # Log file from Nextflow
#  # Other nextflow hidden files, eg. history of pipeline runs and old logs.

nextflow run SciLifeLab/NGI-RNAseq -profile base -c /bioware/lib/wanglab.config --saveReference \
         --project $project --reads "$reads" --fasta $genome --gtf $GTFfile
#nextflow run SciLifeLab/NGI-RNAseq -profile docker --project 'test123' --saveReference \
#         --project $project --reads "$reads" --fasta $genome --gtf $GTFfile
# cd /home/zhenyisong/biodata/wanglab/wangdata/yaoyan/clean/nextflow/mRNAQC
# nextflow.sh -p 'yaoyanQC' -g hg19 -r '../../mRNA/*.{1,2}.clean.fq.gz' 2>nextflow.error &