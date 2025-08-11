#! /bin/bash
#
source $HOME/.bashrc
pushd /home/user/UEDGE
git clean -d -f
git pull
if [ -n "$1" ] 
then
     echo "Executing the extra command: git $@"
     git $@
fi
cd builder
./dsys config linux_docker
./dsys clean
./dsys build 
./dsys load
./dsys test
exit $?


