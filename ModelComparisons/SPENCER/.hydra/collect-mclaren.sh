#!/bin/bash

scp mc:~/hydra-tmp/$4 hydra-tmp/$4.tmp

if [ $? != 0 ]
  then
    exit
fi

echo $4 collected successfully

mv -f hydra-tmp/$4.tmp hydra-tmp/$4

ssh mc "rm -f ~/hydra-tmp/$1"
ssh mc "rm -f ~/hydra-tmp/$2"
ssh mc "rm -f ~/hydra-tmp/$4"
ssh mc "rm -f ~/hydra-tmp/$1.sh"
ssh mc "rm -f ~/hydra-tmp/$1.o"
ssh mc "rm -f ~/hydra-tmp/$1.e"

