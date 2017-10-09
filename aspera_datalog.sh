#!/bin/bash
# nohup /home/zhenyisong/biodata/wanglab/wangcode/aspera_datalog.sh &
#---

#
# data query
#
# 
# data annotation
#


set -euo pipefail
ASPC='/bioware/bin/.aspera/connect/bin/ascp'
KEY='/bioware/bin/.aspera/connect/etc/asperaweb_id_dsa.openssh'
HOST='ftp-private.ncbi.nlm.nih.gov'

# this is for zongna
DATA22='/sra/sra-instant/reads/ByStudy/sra/SRP/SRP106/SRP106502'
DATA23='/sra/sra-instant/reads/ByStudy/sra/SRP/SRP076/SRP076227'

# this is for wangyin
# source GEO
# GSE35583
#---
DATA1='/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189972'
DATA2='/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189966'
DATA3='/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189967'
DATA4='/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189968'
DATA5='/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189969'
DATA6='/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189973'
DATA7='/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189974'
DATA8='/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189975'
DATA9='/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189976'
DATA10='/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189977'
DATA11='/sra/sra-instant/reads/ByExp/sra/SRX/SRX190/SRX190098'
DATA12='/sra/sra-instant/reads/ByExp/sra/SRX/SRX190/SRX190099'
DATA13='/sra/sra-instant/reads/ByExp/sra/SRX/SRX190/SRX190102'
DATA14='/sra/sra-instant/reads/ByExp/sra/SRX/SRX190/SRX190103'
DATA15='/sra/sra-instant/reads/ByExp/sra/SRX/SRX190/SRX190104'
DATA16='/sra/sra-instant/reads/ByExp/sra/SRX/SRX190/SRX190105'
DATA17='/sra/sra-instant/reads/ByExp/sra/SRX/SRX190/SRX190106'
DATA18='/sra/sra-instant/reads/ByExp/sra/SRX/SRX190/SRX190107'
DATA19='/sra/sra-instant/reads/ByExp/sra/SRX/SRX190/SRX190108'
DATA20='/sra/sra-instant/reads/ByExp/sra/SRX/SRX190/SRX190109'

# GSE19090 
# microarray data
#---


# GSE58363
# 
DATA21='/sra/sra-instant/reads/ByStudy/sra/SRP/SRP043/SRP043076'

$ASPC -T -l 200M -k 3 -i $KEY --host=$HOST --user=anonftp --mode=recv --retry-timeout=10 -d ${DATA21} ./

#find . -name "*.sra" -type f | xargs basename -s ".sra"| xargs -I{} -n 1 -P 6 fastq-dump {}.sra | xargs -I{} mv {}.fastq ./
#find . -name "*.sra" -type f | xargs  -n 1 -P 6 fastq-dump
# nohup find . -name '*.sra' -type f | xargs  -n 1 -P 2 fastq-dump --split-3 &
# nohup find . -name '*.sra' -type f | xargs  -n 1 -P 4 fastq-dump &