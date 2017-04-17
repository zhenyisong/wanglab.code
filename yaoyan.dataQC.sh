#!/bin/bash
cd /home/zhenyisong/biodata/wanglab/wangdata/yaoyan
mkdir clean
cd clean
mkdir mRNA
cd mRNA
mRNAcleanPath='/home/zhenyisong/biodata/wanglab/wangdata/yaoyan/rawdata/2015-1-13-mRNA/'

find $mRNAcleanPath -name "*.clean.fq.gz" | xargs -P 4 -i cp {} ./

mkdir -p nextflow/mRNAQC
cd /home/zhenyisong/biodata/wanglab/wangdata/yaoyan/clean/nextflow/mRNAQC

nohup nextflow.sh -p 'yaoyanQC' -g hg19 -r '../../mRNA/*.{1,2}.clean.fq.gz' 2>nextflow.error &
# tail -f nohup.out