#!/bin/bash

# this is the protocol for mRNA
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


# this is the protocol for miRNA
cd /home/zhenyisong/biodata/wanglab/wangdata/yaoyan/clean
mkdir miRNA
cd miRNA
miRNAcleanPath='/home/zhenyisong/biodata/wanglab/wangdata/yaoyan/rawdata/2015-1-28-miRNA/141219_SN7001347_0224_AHB61HADXX/'
#  find $miRNAcleanPath -name "*.clean.fq"|wc -l
find $miRNAcleanPath -name "*.clean.fq" | xargs -P 4 -i cp {} ./ &


