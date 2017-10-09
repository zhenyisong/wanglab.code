#!/bin/bash
# How to install Torque on ubuntu 8.04 
# on a single multicore machine
# https://ubuntuforums.org/showthread.php?t=1372508
# the above title
# 

# https://superuser.com/questions/541500/can-we-use-torque-4-2-0-for-single-machine-acting-as-both-server-and-node
# Caveats: This in and of itself is not practical

# http://stackoverflow.com/questions/28513379/unknown-queue-msg-requested-queue-not-found
# 


# system default setting
##  [zhenyisong wanglab:] $ qmgr -c 'p s'
##  #
##  # Create queues and set their attributes.
##  #
##  #
##  # Create and define queue high
##  #
##  create queue high
##  set queue high queue_type = Execution
##  set queue high max_running = 100
##  set queue high acl_user_enable = False
##  set queue high acl_users = root
##  set queue high resources_default.walltime = 24:00:00
##  set queue high enabled = True
##  set queue high started = True
##  #
##  # Create and define queue middle
##  #
##  create queue middle
##  set queue middle queue_type = Execution
##  set queue middle max_running = 100
##  set queue middle acl_user_enable = False
##  set queue middle acl_users = root
##  set queue middle resources_default.walltime = 24:00:00
##  set queue middle enabled = True
##  set queue middle started = True
##  #
##  # Create and define queue batch
##  #
##  create queue batch
##  set queue batch queue_type = Execution
##  set queue batch max_running = 200
##  set queue batch acl_user_enable = False
##  set queue batch acl_users = root
##  set queue batch resources_default.walltime = 24:00:00
##  set queue batch enabled = True
##  set queue batch started = True
##  #
##  # Create and define queue low
##  #
##  create queue low
##  set queue low queue_type = Execution
##  set queue low max_running = 200
##  set queue low acl_user_enable = False
##  set queue low acl_users = root
##  set queue low resources_default.walltime = 24:00:00
##  set queue low enabled = True
##  set queue low started = True
##  #
##  # Create and define queue defaultApp
##  #
##  create queue defaultApp
##  set queue defaultApp queue_type = Execution
##  set queue defaultApp max_running = 200
##  set queue defaultApp acl_user_enable = False
##  set queue defaultApp acl_users = root
##  set queue defaultApp resources_default.walltime = 24:00:00
##  set queue defaultApp enabled = True
##  set queue defaultApp started = True
##  #
##  # Set server attributes.
##  #
##  set server scheduling = True
##  set server acl_hosts = WangLab
##  set server acl_users = root@*
##  set server acl_roots = root@*
##  set server managers = root@*
##  set server operators = root@*
##  set server default_queue = low
##  set server log_events = 255
##  set server mail_from = adm
##  set server query_other_jobs = True
##  set server scheduler_iteration = 600
##  set server node_check_rate = 150
##  set server tcp_timeout = 90
##  set server job_stat_rate = 120
##  set server poll_jobs = True
##  set server log_level = 7
##  set server owner_purge = False
##  set server mom_job_sync = True
##  set server keep_completed = 300
##  set server submit_hosts = WangLab
##  set server log_file_max_size = 512000
##  set server log_file_roll_depth = 3
##  set server log_keep_days = 30
##  set server next_job_number = 7
##  set server job_start_timeout = 120
##  set server moab_array_compatible = True
##  set server nppcu = 1
##
sudo service trqauthd restart
sudo service pbs_server restart
sudo service pbs_mom restart
##sudo service qterm restart
sudo service maui.d restart

sudo service trqauthd stop
sudo service pbs_server stop
sudo service pbs_mom stop
sudo service qterm restart
sudo service pbs_sched stop


ps -aux | grep pbs #check all is running
qstat -q #check the presence of the queue
qmgr -c 'p s' #check server & queue settings
pbsnodes -a  #check if the nodes are listed and up


sudo qmgr -c "set queue batch max_running = 35"
sudo qmgr -c "set queue batch resources_max.ncpus = 35"
cd /var/lib/torque


sudo qmgr -c "set queue batch max_running = 8"
sudo qmgr -c "set queue batch resources_max.ncpus = 8"
sudo qmgr -c "set queue batch resources_min.ncpus = 1"
sudo qmgr -c "set queue batch resources_max.nodes = 1"
sudo qmgr -c "set queue batch resources_default.ncpus = 1"
sudo qmgr -c "set queue batch resources_default.neednodes = 1:ppn=1"
sudo qmgr -c "set queue batch resources_default.nodect = 1"
sudo qmgr -c "set queue batch resources_default.nodes = 1"