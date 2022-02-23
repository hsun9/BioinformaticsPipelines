
# Hua Sun
# 2021-12-07


# CLAM  https://github.com/Xinglab/CLAM

TYPE='CLIP-Seq'

MULTI_BAM_STR=''
CONTROL_MULTI_BAM_STR=''
while getopts "C:T:N:1:2:3:4:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    T)
      TYPE=$OPTARG
      ;;
    N)
      NAME=$OPTARG
      ;;
    1)
      UNIQ_BAM_STR=$OPTARG
      ;;
    2)
      MULTI_BAM_STR=$OPTARG
      ;;
    3)
      CONTROL_UNIQ_BAM_STR=$OPTARG
      ;;
    4)
      CONTROL_MULTI_BAM_STR=$OPTARG
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

if [ ! -d $OUTDIR ]; then
  echo "[ERROR] The $OUTDIR not exists!" >&2
  exit 1
fi

if [ -z "$NAME" ]; then
  echo "[ERROR] The Name is empty!" >&2
  exit 1
fi


METHOD='start'
if [[ $TYPE == 'RIP-Seq' ]];then
  METHOD='median'
fi


source $CONFIG

OUT=$OUTDIR/$NAME
mkdir -p $OUT

# Negative-binomial model with control

if [[ $MULTI_BAM_STR != '' ]] || [[ $CONTROL_MULTI_BAM_STR != '' ]];then

# Peak calling
${CLAM} peakcaller -i ${UNIQ_BAM_STR} ${MULTI_BAM_STR} \
-c ${CONTROL_UNIQ_BAM_STR} ${CONTROL_MULTI_BAM_STR} \
-o ${OUT} --binsize 100 \
--gtf ${GTF}

else

# Peak calling
${CLAM} peakcaller -i ${UNIQ_BAM_STR} \
-c ${CONTROL_UNIQ_BAM_STR} \
-o ${OUT} --binsize 100 \
--gtf ${GTF}

fi

