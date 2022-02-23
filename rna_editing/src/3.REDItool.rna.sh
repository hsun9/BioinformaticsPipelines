#!/bin/bash

# Hua Sun
# 2021-12-20

# getOptions
while getopts "C:S:B:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    S)
      SAMPLE=$OPTARG
      ;;
    B)
      RNA_BAM=$OPTARG
      ;;
    O)
      OUTDIR=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


source ${CONFIG}

OUT=$OUTDIR/$SAMPLE
mkdir -p $OUT


# 1.call RNA-seq rna-editing 
${REDItoolDnaRna} -i ${RNA_BAM} -o $OUT -f ${GENOME} -t 8 -g 2 -c 8,8 -v 3 -m 20,20 -q 25,25 -a 6-0 -n 0.001 -z -e -u -l -R

