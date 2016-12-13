#!/bin/bash

# 2012-06-13

scp ares:~/hydra-tmp/$4 hydra-tmp/$4.tmp

SCP_STATUS=$?

if [ $SCP_STATUS == 0 ]; then

  echo $4 collected successfully

  mv hydra-tmp/$4.tmp hydra-tmp/$4

  ssh ares "rm -f ~/hydra-tmp/$1     \
                  ~/hydra-tmp/$2     \
                  ~/hydra-tmp/$4     \
                  ~/hydra-tmp/$1.sh  \
                  ~/hydra-tmp/$1*.o   \
                  ~/hydra-tmp/$1*.e"
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
