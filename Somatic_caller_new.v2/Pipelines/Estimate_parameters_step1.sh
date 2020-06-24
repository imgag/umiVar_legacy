#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

SAMTOOLS=samtools

####https://github.com/nhoffman/argparse-bash
##wget https://raw.githubusercontent.com/nhoffman/argparse-bash/master/argparse.bash
##chmod +x argparse.bash
ARGPARSE=$DIR/argparse.bash
source $ARGPARSE || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-bam','--bam', required=True, help='Tumor bam file')
parser.add_argument('-b','--bed', default='', help='Bed file of the targeted regions')
parser.add_argument('-r','--ref', required=True, help='Reference genome - fasta')
parser.add_argument('-o', '--out_file', default='$PWD/Parameters.txt', help='Beta-binomial parameters table out file.')
parser.add_argument('-mq','--mq', default=30, help='Minimum mapping quality')
parser.add_argument('-bq','--bq', default=20, help='Minimum base quality')
EOF

echo Parameters:
echo Tumor bam: "$BAM"
echo Out file: "$OUT_FILE"
echo Bed file: "$BED"
echo Reference file: $REF
echo Beta-binomial parameters table: "$PARAM"
echo Minimum mapping quality: "$MQ"
echo Minimum base quality: "$BQ"
#REF=/mnt/share/data/genomes/GRCh37.fa
#BED=/mnt/share/data/enrichment/HNSCC_cfDNA_hotspot_2018_12_20.bed

SCRIPTS=$DIR/../Tools
OUT=$(dirname $OUT_FILE)

mkdir -p $OUT
TEMP=$OUT/Temp
mkdir -p $TEMP

## 1. Split bam
# Split bam file by barcode
echo python $SCRIPTS/split_bam_step1.py -i $BAM -o $TEMP/DEDUPLICATED 
python $SCRIPTS/split_bam_step1.py -i $BAM -o $TEMP/DEDUPLICATED 


# Index new split bam files
for i in $TEMP/*DP*.bam; do
    $SAMTOOLS index $i
done

## 2. Mpileup files
# Check if bed file is provided
if [[ ! -z "$BED" ]];then
    for i in $TEMP/*DP*.bam; do
        NAME=$(echo $i | sed 's/bam$/pileup/')
        $SAMTOOLS mpileup -d 100000000 -f $REF \
            -l $BED \
            -q $MQ \
            -Q 1 \
            -x \
            $i > $NAME &   
    done
    wait
else
    for i in $TEMP/*DP*.bam; do
        NAME=$(echo $i | sed 's/bam$/pileup/')
        $SAMTOOLS mpileup -d 100000000 -f $REF \
            -q $MQ \
            -Q 1 \
            -x \
            $i > $NAME &   
    done
    wait
fi

## 3. Mpileup to tsv
for i in $TEMP/*DP*.pileup; do
    NAME2=$(echo $i | sed 's/pileup$/tsv/')
    python $SCRIPTS//pileup2tsv_08_01_19.py --pileup $i \
        -o $NAME2 \
        --minBQ $BQ \
        -m 10 &
done
wait

## 4. Parameters
Rscript $SCRIPTS/Parameters_step1.r -t1 $TEMP/DEDUPLICATED_DP1.tsv \
    -t2 $TEMP/DEDUPLICATED_DP2.tsv \
    -t3 $TEMP/DEDUPLICATED_DP3.tsv \
    -t4 $TEMP/DEDUPLICATED_DP4.tsv \
    -o $OUT_FILE
