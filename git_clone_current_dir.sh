# http://stackoverflow.com/questions/651038/how-do-you-clone-a-git-repository-into-a-specific-folder
# https://gist.github.com/rakuishi/eeeeb75ab1cc3bc062c02b8ffd8381b0
git init
git remote add -t \* -f origin https://github.com/zhenyisong/wanglab.code.git
git checkout master

# enter the github store without password
# you should copy id_rsa.pub to the github
# https://github.com/settings/keys
# the IP adddress should be updated to the outside address
# 123.124.148.35
# the above IP address is the machine wanglab IP address
git remote set-url origin ssh://git@github.com/zhenyisong/wanglab.code.git