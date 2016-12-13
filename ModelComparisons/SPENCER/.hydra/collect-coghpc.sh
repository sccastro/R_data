#!/bin/bash

# 2011-10-29

scp ares:~/hpc/hydra-tmp/$4 hydra-tmp/$4.tmp

SCP_STATUS=$?

if [ $SCP_STATUS == 0 ]; then

  echo $4 collected successfully

  mv hydra-tmp/$4.tmp hydra-tmp/$4

  ssh ares "rm -f ~/hpc/hydra-tmp/$1     \
                  ~/hpc/hydra-tmp/$2     \
                  ~/hpc/hydra-tmp/$4     \
                  ~/hpc/hydra-tmp/$1.sh  \
                  ~/hpc/hydra-tmp/$1*.o   \
                  ~/hpc/hydra-tmp/$1*.e"
  exit 0
fi

QSTAT_OUTPUT=`ssh ares "source /grid/sge6/default/common/settings.sh ; /grid/sge6/bin/lx24-amd64/qstat -u \\\$USER | grep -v Eqw | tail +3"`
QSTAT_STATUS=$?

if [ $QSTAT_STATUS != 0 ]; then
  exit 2
fi

IFS=$'\n'

for line in $QSTAT_OUTPUT
do
  if [[ ${line:16:10} == $1 ]]
  then
    exit 2
    break
  fi
done

exit 1
