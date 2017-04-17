# run this will result in creation of id_rsa.pub
# and id_rsa
# in $HOME/.ssh/, after generation this file
# please use vi to edit it using IP address to
# 

# this is created at wanglab side
# we finaly want to log into the fuwai side
# we upload the public key at wanglab side
# to the authorized_file at fuwai side
#--

ssh-keygen -b 2048
#
# cd ~/.ssh/
# vi id_rsa.pub and change the @ to IP address.
# for example, wanglab, 172.31.205.11

# permission requirement
# chmod 700 /home/user
# chmod 700 ~/.ssh
# chmod 600 ~/.ssh/authorized_keys
# chmod 600 ~/.ssh/config
# chmod 600 ~/.ssh/privatekey
# chmod 644 ~/.ssh/publickey.pub

ssh-copy-id -i ~/.ssh/id_rsa.pub  zhenyisong@172.31.202.11