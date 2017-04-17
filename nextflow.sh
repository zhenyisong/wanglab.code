#!/bin/bash
# https://github.com/ewels/NGI-RNAseq
#cd /temp
random_projectID=`mktemp -u  'NextflowQC_XXXXXX`
echo "your project ID is ${random_projectID} created for this case, please write down in your worksheet\n"
mm10_UCSC_genome='/mnt/date/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa'
mm9_UCSC_genome='/mnt/date/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa'
hg19_UCSC_genome='/mnt/date/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
hg38_UCSC_genome=''
mm10_UCSC_GTF='/mnt/date/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/gene.gtf'
mm9_UCSC_GTF='/mnt/date/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/gene.gtf'
hg19_UCSC_GTF='/mnt/date/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf'
hg38_UCSC_GTF=''

#nextflow run SciLifeLab/NGI-RNAseq -profile base --project 'test123' --reads 'reads/*.fastq' --fasta '/mnt/date/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa' --gtf '/mnt/date/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'