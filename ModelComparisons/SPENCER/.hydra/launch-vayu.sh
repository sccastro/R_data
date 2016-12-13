#!/bin/bash

echo launch-vayu.sh version 1.3

PROJECT=y52

status=0

ssh vayu "mkdir /short/$PROJECT/hydra-tmp 2>/dev/null"
status=$?
if [ status = 0 ]
then
  echo "hydra-tmp created"
fi

for i in 1 2 3 4 5 6 7 8 9
do
  #echo copying $1
  scp hydra-tmp/$1 vayu:/short/$PROJECT/hydra-tmp/$1
  status=$?
  if [ $status -eq 0 ]
  then
#    echo $1 copied successfully
    break
  else
    echo $1 could not be copied, error code $status
    sleep 5
  fi
done

if [ $status -ne 0 ]
then
  echo after $i attempts, giving up
  exit $status
fi


for i in 1 2 3 4 5 6 7 8 9
do
#  echo copying $2
  scp hydra-tmp/$2 vayu:/short/$PROJECT/hydra-tmp/$2
  status=$?
  if [ $status -eq 0 ]
  then
#    echo $2 copied successfully
    break
  else
    echo $2 could not be copied, error code $status
    sleep 5
  fi
done

if [ $status -ne 0 ]
then
  echo after $i attempts, giving up
  exit $status
fi

for i in 1 2 3 4 5 6 7 8 9
do
#  echo writing $1.sh

  ssh vayu "echo -e \#!/bin/bash \\\n\
    set -e \\\n\
    module load R/2.11.1 \\\n\
    cd /short/$PROJECT/hydra-tmp \\\n\
    Rscript $1 \;\\\n\
    mv /short/$PROJECT/hydra-tmp/$3 /short/$PROJECT/hydra-tmp/$4 \\\n\
    exit \\\n\
  > /short/$PROJECT/hydra-tmp/$1.sh"

  status=$?

  if [ $status -eq 0 ]
  then
#    echo $1.sh written successfully
    break
  else
    echo $1.sh could not be written, error code $status
    sleep 5
  fi
done

if [ $status -ne 0 ]
then
  echo after $i attempts, giving up
  exit $status
fi


for i in 1 2 3 4 5 6 7 8 9
do
#  echo submitting $1.sh

  est_time=$5
  if [[ $est_time == "" ]]
  then
    est_time=3600
  fi

  options=$6
  if [[ $options == "" ]]
  then
    options="vmem=500MB"
  fi

  ssh vayu "/opt/pbs/bin/qsub -wd -q normal -l walltime=$est_time,$options -N $1 -o /short/$PROJECT/hydra-tmp/$1-$$.o -e /short/$PROJECT/hydra-tmp/$1-$$.e /short/$PROJECT/hydra-tmp/$1.sh"

  status=$?

  if [ $status -eq 0 ]
  then
#    echo $1.sh submitted successfully
    break
  else
    echo could submit job $1.sh, error code $status
    sleep 5
  fi
done

if [ $status -ne 0 ]
then
  echo after $i attempts, giving up
  exit $status
fi

echo Job started successfully

exit 0

