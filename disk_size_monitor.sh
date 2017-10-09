#!/bin/bash
set -u
set -e
set -o pipefailed
# this is for checking the disk size and send error message if 
# the disk size is overloaded in

# df -ah
# 192.168.1.11:/date        51T   29T   23T  57% /mnt/date
# /dev/mapper/centos-root  1.5T  323G  1.1T  23% /
LIMIT=25
cutoff=`df -ah /mnt/date | grep -Po '(\d+)\%' | sed -r 's/\%//'`
if [[ $cutoff gt $LIMIT ]]; then
    mail -s 'dangerous disk overflow, please erase them!' zhenyisong@hotmail.com
    exit 0
fi