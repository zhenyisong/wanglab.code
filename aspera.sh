#!/bin/bash
set -e
set -u
set -o pipefail
# need change the mod
# chmod 755
# create cro
# crontab -e
# 0 0-23/2 * * * /home/zhenyisong/biodata/wanglab/wangcode/aspera.sh
# crontab -l
# @version
#   1.001
# @since 
#   2016-05-18
# @update
#   2016-05-18
# @ reference site
#   https://gist.github.com/cwvhogue/735e78b88a1b19243186
#   http://www.chenlianfu.com/?p=2319
#   http://blog.sciencenet.cn/home.php?mod=space&uid=1224852&do=blog&id=874249
#   http://www.plob.org/2012/07/31/3013.html
#   http://bencane.com/2015/09/22/preventing-duplicate-cron-job-executions/
# @ how to install
#   wget -c -t 0 http://download.asperasoft.com/download/sw/connect/3.5/aspera-connect-3.5.1.92523-linux-64.sh
#   sh aspera-connect-3.5.1.92523-linux-64.sh
# @specify the enviroment
#   deprecated:
#   echo 'export PATH=/home/zhenyisong/.aspera/connect/bin/:$PATH' >> .bash_profile
#   echo 'export aspera_license=/home/zhenyisong/.aspera/connect/etc/asperaweb_id_dsa.openssh' >> .bash_profile
# @parameters
#    -l max_rate	
#       Set the target transfer rate in Kbps (default: 10000 Kbps). 
#       If the ascp client does not specify a target rate, it will 
#       be acquired from aspera.conf (server-side, as the local aspera.
#       conf target rate setting doesn't apply). If local or server 
#       aspera.conf rate caps are specified, the "starting" (default) 
#       rates will be not higher than the cap.
#    -T	
#       Disable encryption for maximum throughput.
#    -k{0|1|2|3}	
#       Enable resuming partially transferred files at the 
#       specified resume level (default: 0). Note that this 
#       must be specified for your first transfer; otherwise, 
#       it will not work for subsequent transfers. Resume levels:
#       0 - Always retransfer the entire file.
#       1 - Check file attributes and resume if the current and original 
#            attributes match.
#       2 - Check file attributes and do a sparse file checksum; 
#            resume if the current and original attributes/checksums match.
#       3 - Check file attributes and do a full file checksum; 
#           resume if the current and original attributes/checksums match.
#       Note that when a complete file exists at the destination (no .aspx), 
#       the source file size is compared with the destination file size. 
#       When a partial file and a valid .aspx file exist at the destination, 
#       the source file size is compared with the file size recorded 
#       inside the .aspx file.
#     -i private_key_file	
#        Use public key authentication and specify the private key file. 
#        Typically, the private key file is in the directory
#     -d	
#        Create target directory if it doesn't already exist.
#     -Q 
#        Enables the fair transfer policy, which ensures that the 
#        available bandwidth is shared amongst other traffic and 
#        transfers at a fair rate

#     defalut usage is from NCBI SRA database
#     alternatively,we should use EBI mirror
#     --host=fasp.sra.ebi.ac.uk --user=era-fasp
#     

# full command line of ascp usage
# you need to specify the file directory and host/user
# the output is set at the current working directory


# this is an ugly commandline, too loog parameters.
# I prune the length using bash enviroment and 
# define the private variables here to short the length
#
# /home/zhenyisong/.aspera/connect/bin/ascp -T -l 200M -k 3 -i \
# /home/zhenyisong/.aspera/connect/etc/asperaweb_id_dsa.openssh \
#  --host=ftp-private.ncbi.nlm.nih.gov --user=anonftp --mode=recv \
#  -d /sra/sra-instant/reads/ByStudy/sra/SRP/SRP006/SRP006947/ ./

# aspera parameters
# you need to update these two parameters

NCBI_DIR='/sra/sra-instant/reads/ByStudy/sra/SRP/SRP033/SRP033009/'
OUTPUT_DIR='/home/zhenyisong/biodata/cardiodata'


# to shorten the command line code.
ASCP='/bioware/bin/.aspera/connect/bin/ascp'
PASS_WORD='/bioware/bin/.aspera/connect/etc/asperaweb_id_dsa.openssh'

PIDFILE=/home/zhenyisong/aspera.pid
if [ -f $PIDFILE ]
then
  PID=$(cat $PIDFILE)
  ps -p $PID > /dev/null 2>&1
  # using the $? special variable to check the 
  # exit code of the last command executed
  if [ $? -eq 0 ]
  then
    echo "Job is already running"
    exit 1
  else
    ## Process not found assume not running
    # The $$ variable, is a special variable within 
    # BASH that returns the current process ID.
    echo $$ > $PIDFILE
    if [ $? -ne 0 ]
    then
      echo "Could not create PID file"
      exit 1
    fi
  fi
else
  echo $$ > $PIDFILE
  if [ $? -ne 0 ]
  then
    echo "Could not create PID file"
    exit 1
  fi
fi



$ASCP -T -l 200M -k 3 -i $PASS_WORD --host=ftp-private.ncbi.nlm.nih.gov --user=anonftp --mode=recv -d $NCBI_DIR $OUTPUT_DIR

rm $PIDFILE

