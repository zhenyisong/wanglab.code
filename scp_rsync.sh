# report the disk size of the specified dir in human readable format
du -sh cardiodata

# if transfer large files , scp is unsafe, as it cannot resume the 
# copy if corrupted, instead, using rsync
# pwd
# cd /home/zhenyisong/data
# https://bugzilla.samba.org/show_bug.cgi?id=5478
nohup scp -r cardiodata/* 172.31.205.11:/home/zhenyisong/biodata/cardiodata 2>scp.log&
nohup rsync -avr --no-compress --inplace --partial -e ssh cardiodata/* 172.31.205.11:/home/zhenyisong/biodata/cardiodata &
