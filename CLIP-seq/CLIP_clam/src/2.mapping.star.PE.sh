#!/bin/bash

# Hua Sun
# 2021-05-21 updated
# memory 64 Gb in RIS-server

TAG="star"
NAME=''
FQ1=''
FQ2=''
while getopts "C:T:N:1:2:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    T)
      TAG=$OPTARG
      ;;
    N)
      NAME=$OPTARG
      ;;
    1)
      FQ1=$OPTARG
      ;;
    2)
      FQ2=$OPTARG
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

if [ ! -e $CONFIG ]; then
  echo "[ERROR] The $CONFIG not exists!" >&2
  exit 1
fi

if [ ! -e $FQ1 ]; then
  echo "[ERROR] The $FQ1 not exists!" >&2
  exit 1
fi

if [ ! -e $FQ2 ]; then
  echo "[ERROR] The $FQ2 not exists!" >&2
  exit 1
fi

if [ ! -d $OUTDIR ]; then
  echo "[ERROR] The $OUTDIR not exists!" >&2
  exit 1
fi

if [ -z "$NAME" ]; then
  echo "[ERROR] The Name is empty!" >&2
  exit 1
fi

source $CONFIG

OUT=$OUTDIR/$NAME
mkdir -p $OUT



# Reference
# https://pureclip.readthedocs.io/en/latest/GettingStarted/preprocessing.html#
#STAR stalled at the very beginning (cannot create fifo to read fastq.gz ?) #687
zcat $FQ1 > $OUT/$NAME.R1.fastq
zcat $FQ2 > $OUT/$NAME.R2.fastq


$STAR --runThreadN 10 \
--genomeDir $STAR_GENOME_DIR \
--readFilesIn $OUT/$NAME.R1.fastq $OUT/$NAME.R2.fastq \
--outSAMtype BAM Unsorted \
--alignEndsType EndToEnd \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterType BySJout \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \
--outFilterMultimapNmax 100 \
--scoreDelOpen -1 \
--genomeLoad NoSharedMemory \
--outFileNamePrefix $OUT/$TAG.


rm -f $OUT/$NAME.R1.fastq $OUT/$NAME.R2.fastq


# Remove rRNAs
$BEDTOOLS intersect -f 0.90 -abam $OUT/star.Aligned.out.bam -b ${rRNA_annotation} -v > $OUT/star.aligned.rm_rRNA.bam

rm -f $OUT/star.Aligned.out.bam


# remove dup.
$SAMBAMBA markdup -t 8 -r -l 5 --io-buffer-size 512 --overflow-list-size 1000000 --tmpdir $TMPDIR $OUT/star.aligned.rm_rRNA.bam $OUT/star.aligned.rm_rRNA.remDup.bam
# sort
$SAMTOOLS sort -m 8G -@8 -o $OUT/star.aligned.rm_rRNA.remDup.sorted.bam $OUT/star.aligned.rm_rRNA.remDup.bam
# Create indices for all the bam files for visualization and QC
$SAMTOOLS index $OUT/star.aligned.rm_rRNA.remDup.sorted.bam


rm -f $OUT/star.aligned.rm_rRNA.bam $OUT/star.aligned.rm_rRNA.remDup.bam

