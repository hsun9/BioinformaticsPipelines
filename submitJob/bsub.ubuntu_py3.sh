#!/bin/sh

# Hua Sun

## USAGE
## sh bsub.sh <memory int> <threads int> <any name> <any command>


MEM=$1
THREADS=$2
NAME=$3

DIR=`pwd`
LOG=$DIR/logs
mkdir -p $LOG

bsub -G mackgrp -q general -M ${MEM}000000 -n ${THREADS} -R "select[mem>${MEM}000] span[hosts=1] rusage[mem=${MEM}000]" -oo ${LOG}/${NAME}.log -eo ${LOG}/${NAME}.err -a "docker(hsun9/ubuntu_xenial_py39)" "$@"


