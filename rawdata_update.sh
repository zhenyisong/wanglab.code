#!/bin/bash
#set -e
set -u
set -o pipefail
# @author Yisong Zhen
# @since 2017-03-13
# @update 2017-03-13
# this is for cronb and update all new seqeuncing data in 
# source machine: 172.31.120.44
# running worker is user1
# destination machine: 172.31.205.11
# destiniation :username is -- rawdata -- runing robotic
# @result
#   capture all the error output and send to admin to fix it
#   sychronize the raw data in soource machine and destiniation machine
#---

# option -a: enables wrsync's archive mode
#        -v makes rsync's progress more verbose so we can see what's been transferred
#        -z: enalbes file transfer compression
#        -e: connecting the remote host throught SSH
#        -r:recurse into directories
#        -q:suppress non-error messages
#        -t:preserve modification times
# rsync -avz -e ssh /date/user1/Sequencing/*.temp rawdata@172.31.205.11:/mnt/date/Sequencing/RawData/

# crontab -e
# I edited the cron configuration file
# which at /etc/crontab
# and add the PATH: /bioware/bin
# MAILTO="" not root here
# https://help.ubuntu.com/community/CronHowto
# http://stackoverflow.com/questions/5200551/how-to-set-a-cron-job-to-run-at-a-exact-time
# http://stackoverflow.com/questions/4114716/rsync-permissions-question-destination-perms-not-properly-applying
# http://serverfault.com/questions/44400/run-a-shell-script-as-a-different-user


# @bug1
#   I have to create a unique file for any special user to create his/her own output log.
#   Otherwise, i have to chmod the log file, this is dangerous using sudo fucntion
# @bug2
#   when synchrinze the data, the permission need to be reset to rxw
# @bug3
#   I cancelled the usage of set -e setting, this will result in no-error throw
#   out catched by following command
#   if we debug the progam, we can open this settting
# @bug4 (?)
#  I have found that on large files (14G-ish or more) the z option (compress) messes up the transfer. 
#  I've had the exact symptoms as the poster happen on a few computers and
#  removing the -z option solved the problem. 
#  Please see this post:
#  http://www.linuxquestions.org/questions/linux-networking-3/rsync-error-23-a-922673/
# @bug5
#    sudo
#    chown -R rawdata
#    chgrp -R rawdata
#    I am not sure if using the -O parameter will result
#    http://stackoverflow.com/questions/667992/rsync-error-failed-to-set-times-on-foo-bar-operation-not-permitted
#-----
rsync -avzr -e ssh --chmod=a=wrx --perms user1@172.31.120.44:/date/user1/Sequencing/* /mnt/date/Sequencing/RawData/ > /tmp/rsync_${USER}.log 2>&1
if [ $? -ne 0 ]; then
    echo "Notice! rawdata transfer failed, please check the server side at 172.31.205.11 (destination) and 172.31.120.44(source)" | mail -s "Alert! rsync falied!" < /tmp/rsync_${USER}.log zhenyisong@hotmail.com
    #echo "Notice! rawdata transfer failed, please check the server side at 172.31.205.11 (destination) and 172.31.120.44(source)" | mail -s "Alert! rsync falied!" yup@pku.edu.cn
else
    echo "Robotic have finished mirroring raw data from Next-500." | mail -s "Good, file transfer successfully completed!" zhenyisong@hotmail.com
    #echo "Robotic have finished mirroring raw data from Next-500." | mail -s "Good, file transfer successfully completed!" yup@pku.edu.cn
fi    