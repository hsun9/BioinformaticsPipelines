
# Hua Sun
# 2021-11-10


# CLAM  https://github.com/Xinglab/CLAM



TYPE='CLIPSEQ'  # CLIPSEQ / RIPSEQ
Multimap=100
while getopts "C:T:N:B:M:O:" opt; do
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
    B)
      BAM=$OPTARG
      ;;
    M)
      Multimap=$OPTARG
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

if [ ! -e $BAM ]; then
  echo "[ERROR] The $BAM not exists!" >&2
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
if [[ $TYPE == 'RIPSEQ' ]];then
  METHOD='median'
fi


source $CONFIG

OUT=$OUTDIR/$NAME
mkdir -p $OUT


# CLAM preprocessor
${CLAM} preprocessor -i ${BAM} -o ${OUT} --read-tagger-method ${METHOD} --max-multihits ${Multimap}

