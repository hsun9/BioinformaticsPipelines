#!/bin/bash

# Hua Sun
# 2021-10-18



CONFIG='config.htseq.ini'
## getOptions
while getopts "C:N:O:1:2:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    N)
      NAME=$OPTARG
      ;;
    O)
      OUTDIR=$OPTARG
      ;;
    1)
      FQ1=$OPTARG
      ;;
    2)
      FQ2=$OPTARG
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


OUT=$OUTDIR/$NAME
mkdir -p $OUT


# align sequencing reads to the genome
# to solve ??????????? ? ?       ?               ?            ? tmp.fifo.read1
# DO NOT use  '--readFilesCommand zcat' 
zcat $FQ1 > $OUT/${NAME}.R1.fastq
zcat $FQ2 > $OUT/${NAME}.R2.fastq


# GDC pipeline - Dr15plus
# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/


# https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
# --outSAMtype BAM Unsorted
# This ”unsorted” file can be directly used with downstream software such as HTseq, without the need of name sorting.
# I removed below commands from original 'GDC pipeline'
# --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
# --quantMode TranscriptomeSAM GeneCounts \
${STAR} --runThreadN 8 \
--genomeDir ${STAR_GENOME_DIR} \
--sjdbGTFfile ${GTF} \
--readFilesIn $OUT/${NAME}.R1.fastq $OUT/${NAME}.R2.fastq \
--outSAMattrRGline "ID:${NAME}\tSM:${NAME}\tLB:${NAME}.RNAseq\tPU:${NAME}\tPL:illumina" \
--alignIntronMax 1000000 \
--alignIntronMin 20 \
--alignMatesGapMax 1000000 \
--alignSJDBoverhangMin 1 \
--alignSJoverhangMin 8 \
--alignSoftClipAtReferenceEnds Yes \
--chimJunctionOverhangMin 15 \
--chimMainSegmentMultNmax 1 \
--chimSegmentMin 15 \
--genomeLoad NoSharedMemory \
--limitSjdbInsertNsj 1200000 \
--outFilterIntronMotifs None \
--outFilterMatchNminOverLread 0.33 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--outFilterMultimapNmax 20 \
--outFilterScoreMinOverLread 0.33 \
--outFilterType BySJout \
--outSAMattributes NH HI AS nM NM ch \
--outSAMstrandField intronMotif \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--twopassMode Basic \
--outFileNamePrefix $OUT/star.


# remove temp fastq
rm -f $OUT/${NAME}.R1.fastq $OUT/${NAME}.R2.fastq


# Run htseq-count and calculate gene-level counts
# https://htseq.readthedocs.io/en/release_0.9.1/count.html
${HTSEQ_COUNT} \
-f bam \
-r name \
-s no \
-a 10 \
-t exon \
-i gene_id \
-m intersection-nonempty \
$OUT/star.Aligned.out.bam \
${GTF} > $OUT/htseq_count.out
