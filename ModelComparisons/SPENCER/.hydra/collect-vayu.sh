#!/bin/bash

PROJECT=y52

scp vayu:/short/$PROJECT/hydra-tmp/$4 hydra-tmp/$4.tmp

SCP_STATUS=$?

if [ $SCP_STATUS == 0 ]; then

  echo $4 collected successfully

  mv -f hydra-tmp/$4.tmp hydra-tmp/$4

  ssh vayu "rm -f /short/$PROJECT/hydra-tmp/$1 ;    \
            rm -f /short/$PROJECT/hydra-tmp/$2 ;    \
            rm -f /short/$PROJECT/hydra-tmp/$4 ;    \
            rm -f /short/$PROJECT/hydra-tmp/$1.sh ; \
            rm -f /short/$PROJECT/hydra-tmp/$1*.o ; \
            rm -f /short/$PROJECT/hydra-tmp/$1*.e   "
  exit 0
fi

QSTAT_OUTPUT=`ssh vayu "qstat -u \\\$USER | grep -v Eqw | tail +3"`
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
