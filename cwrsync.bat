@ECHO ON
REM *****************************************************************
REM
REM CWRSYNC.CMD - Batch file template to start your rsync command (s).
REM
REM *****************************************************************

REM Make environment variable changes local to this batch file
SETLOCAL

REM Specify where to find rsync and related files
REM Default value is the directory of this batch file
REM This need to be user specific, I changed it according
REM my computer setting, this defined the working dir
REM of the program downloaded from
REM https://www.itefix.net/cwrsync
REM SHA256: 37e8ef21ac975d4ee86c9d3be40c8935e8b9d0ba84e9302fc106b9452296cb85
REM version   5.5.0_x86
REM and unzipped the FREE edition
REM no quoted, otherwise will throw out error!
SET CWRSYNCHOME=C:\Users\Yisong\Desktop\cwRsync

REM in the bin dir and run the code to generate 
REM fingerprint file id_rsa & id_rsa.pub
REM
REM I must use this command to generate these files
REM ssh-keygen.exe -N "" -f id_rsa -C "wanglab_window"

REM run the code below and create the corresponding dir
REM please uncommmented the phrase below
REM Create a home directory for .ssh 
REM IF NOT EXIST %CWRSYNCHOME%\home\%USERNAME%\.ssh MKDIR %CWRSYNCHOME%\home\%USERNAME%\.ssh
REM enter into the home/yisong/.ssh dir
REM and cut & paste the generated id_rsa & id_rsa.pub



REM Make cwRsync home as a part of system PATH to find required DLLs
SET CWOLDPATH=%PATH%
SET PATH=%CWRSYNCHOME%\bin;%PATH%

REM Windows paths may contain a colon (:) as a part of drive designation and 
REM backslashes (example c:\, g:\). However, in rsync syntax, a colon in a 
REM path means searching for a remote host. Solution: use absolute path 'a la unix', 
REM replace backslashes (\) with slashes (/) and put -/cygdrive/- in front of the 
REM drive letter:
REM 
REM Example : C:\WORK\* --> /cygdrive/c/work/*
REM 
REM Example 1 - rsync recursively to a unix server with an openssh server :
REM
REM       rsync -r /cygdrive/c/work/ remotehost:/home/user/work/
REM
REM Example 2 - Local rsync recursively 
REM
REM       rsync -r /cygdrive/c/work/ /cygdrive/d/work/doc/
REM
REM Example 3 - rsync to an rsync server recursively :
REM    (Double colons?? YES!!)
REM
REM       rsync -r /cygdrive/c/doc/ remotehost::module/doc
REM
REM Rsync is a very powerful tool. Please look at documentation for other options. 
REM

REM ** CUSTOMIZE ** Enter your rsync command(s) here


REM remember to change the permission of files
REM only this user is cared even the root must be excluded from the
REM the read/write right

REM --------------------------------------
REM added by Yisong
REM you need specify the source machine address
REM dstiniation machine already is specified
REM

rsync -vzr -I /cygdrive/d/tempdata/ rawdata@172.31.205.11:/mnt/date/Sequencing/RawData/

REM --------------------------------------

REM you must install 
REM https://www.activestate.com/activeperl/downloads
REM baidu_yun??
REM
REM install the sendEmail package
REM software version
REM sendEmail-v156
REM download from 
REM http://caspian.dotconf.net/menu/Software/SendEmail/
REM http://caspian.dotconf.net/menu/Software/SendEmail/sendEmail-v156.zip
REM ppm install/upgrade Net::SSLeay
REM ppm  IO::Socket::SSL

if errorlevel 1 (
   echo Failure Reason Given is %errorlevel%
   sendEmail -f zhenyisong@fuwaihospital.org -t zhenyisong@hotmail.com ^
   -o tls=auto -s mail.fuwaihospital.org:25 -xu zhenyisong@fuwaihospital.org ^
   -xp 789uio2008 -u "Window raw data tranfer failed! " -m "you should check the next-500"
   exit /b %errorlevel%
)
else (
   sendEmail -f zhenyisong@fuwaihospital.org -t zhenyisong@hotmail.com ^
   -o tls=auto -s mail.fuwaihospital.org:25 -xu zhenyisong@fuwaihospital.org ^
   -xp 789uio2008 -u "Window raw data tranfer succeed! " -m "you can sleep next-500"
   exit /b %errorlevel%
)
