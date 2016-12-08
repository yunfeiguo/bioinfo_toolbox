#!/bin/bash
#ssh -f -N guoy28@10.27.4.2 -L 1986:guoy28@tehran-4:22 
sudo umount $HOME/Desktop/bina_t-rex_home
#sudo sshfs -p 1986 -o follow_symlinks,allow_other,defer_permissions,IdentityFile=$HOME/.ssh/id_rsa guoy28@127.0.0.1:/net/hippo/volumes/wadi/shared/prj/guoy28 $HOME/Desktop/bina_t-rex_home
sudo sshfs -o follow_symlinks,allow_other,defer_permissions,IdentityFile=$HOME/.ssh/id_rsa guoy28@10.27.4.2:/net/hippo/volumes/wadi/shared/prj/guoy28 $HOME/Desktop/bina_t-rex_home
