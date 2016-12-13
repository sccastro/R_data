#!/bin/bash

echo launch-local.sh version 1.1

cd hydra-tmp
Rscript $1
mv $3 $4


