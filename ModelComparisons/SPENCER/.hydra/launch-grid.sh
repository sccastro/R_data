#!/bin/bash

echo launch-grid.sh version 1.4

status=0

ssh ares "mkdir ~/hpc/hydra-tmp 2>/dev/null"
status=$?
if [ status = 0 ]
then
  echo "hydra-tmp created"
fi

for i in 1 2 3 4 5 6 7 8 9
do
  echo copying $1
  scp hydra-tmp/$1 ares:~/hpc/hydra-tmp/$1
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
  echo copying $2
  scp hydra-tmp/$2 ares:~/hpc/hydra-tmp/$2
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



#echo writing run.sh
ssh ares "echo bash \\\$1 > ~/hpc/hydra-tmp/run.sh"
#echo chmodding run.sh
ssh ares "chmod 700 ~/hpc/hydra-tmp/run.sh"

for i in 1 2 3 4 5 6 7 8 9
do
#  echo writing $1.sh

  ssh ares "echo -e \#!/bin/bash \\\n\
    set -e \\\n\
    cd ~/hpc/hydra-tmp \\\n\
    if [ -f /usr/bin/Rscript ]\; then \\\n\
      /usr/bin/Rscript $1 \;\\\n\
    elif [ \\\"\\\`uname\\\`\\\" = \\\"SunOS\\\" ]\; then \\\n\
      /opt/sfw/R-2.7.1/bin/Rscript $1 \;\\\n\
    else \\\n\
      /grid/local/linux-64bit/R/R-2.10.0/bin/Rscript $1 \;\\\n\
    fi \\\n\
    rm -f ~/hpc/hydra-tmp/$1  \\\n\
    rm -f ~/hpc/hydra-tmp/$2  \\\n\
    mv ~/hpc/hydra-tmp/$3 ~/hpc/hydra-tmp/$4 \\\n\
    exit \\\n\
  > ~/hpc/hydra-tmp/$1.sh"

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

  ssh ares "set -e ; source /grid/sge6/default/common/settings.sh ; /grid/sge6/bin/lx24-amd64/qsub -R yes -clear -q all.q -q es409.q -q fist.q -N $1 -o ~/hpc/hydra-tmp/$1-$$.o -e ~/hpc/hydra-tmp/$1-$$.e ~/hpc/hydra-tmp/run.sh ~/hpc/hydra-tmp/$1.sh"

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

