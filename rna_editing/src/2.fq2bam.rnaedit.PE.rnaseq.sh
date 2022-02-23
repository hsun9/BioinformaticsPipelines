#!/bin/bash

# Hua Sun
# 2022-02-21 v0.2

# getOptions
while getopts "C:S:1:2:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    S)
      SAMPLE=$OPTARG
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
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done


source ${CONFIG}


OUT=$OUTDIR/$SAMPLE
mkdir -p $OUT


# align sequencing reads to the genome
# to solve ??????????? ? ?       ?               ?            ? tmp.fifo.read1
# DO NOT use  '--readFilesCommand zcat' 
zcat $FQ1 > $OUT/${SAMPLE}.R1.fastq
zcat $FQ2 > $OUT/${SAMPLE}.R2.fastq


# ref-1:Investigating RNA editing in deep transcriptome datasets with REDItools and REDIportal
# ref-2: https://hal.archives-ouvertes.fr/hal-03080231/document
${STAR} --runThreadN 8 \
--genomeDir ${STAR_GENOME_DIR} \
--sjdbGTFfile ${GTF} \
--readFilesIn $OUT/${SAMPLE}.R1.fastq $OUT/${SAMPLE}.R2.fastq \
--genomeLoad NoSharedMemory \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--twopassMode Basic \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outSAMattributes All \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \
--outFilterMultimapNmax 1 \
--limitBAMsortRAM 32000000000 \
--outFileNamePrefix $OUT/star.

# remove temp fastq
rm -f $OUT/${SAMPLE}.R1.fastq $OUT/${SAMPLE}.R2.fastq

# sort
$JAVA -Xmx16G -jar $PICARD SortSam \
    -CREATE_INDEX true \
    -I $OUT/star.Aligned.out.bam \
    -O $OUT/$SAMPLE.sorted.bam \
    -SORT_ORDER coordinate \
    -VALIDATION_STRINGENCY STRICT

# remove bam for save space
rm -f $OUT/star.Aligned.out.bam

# remove-duplication
$JAVA -Xmx16G -jar $PICARD MarkDuplicates \
    -I $OUT/$SAMPLE.sorted.bam \
    -O $OUT/$SAMPLE.remDup.bam \
    --REMOVE_DUPLICATES true \
    -M $OUT/$SAMPLE.remdup.metrics.txt


# index bam
$SAMTOOLS index $OUT/$SAMPLE.remDup.bam

# remove sorted bam
rm -f $OUT/$SAMPLE.sorted.bam
