#!/bin/bash

SAMTOOLS=/mnt/users/ahmuyaf1/programs/samtools-1.3//samtools


####https://github.com/nhoffman/argparse-bash
##wget https://raw.githubusercontent.com/nhoffman/argparse-bash/master/argparse.bash
##chmod +x argparse.bash
ARGPARSE=/mnt/users/ahmuyaf1/programs/argparse/argparse.bash
source $ARGPARSE || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-bam','--bam', required=True, help='Tumor bam file')
parser.add_argument('-nbam','--nbam', default='', required=False, help='Normal bam file')
parser.add_argument('-b','--bed', default='', help='Bed file of the targeted regions')
parser.add_argument('-r','--ref', required=True, help='Reference genome - fasta')
parser.add_argument('-o', '--out_file', default='$PWD/Sample.vcf', help='Out vcf file')
parser.add_argument('-p', '--param', default='', help='Beta-binomial parameters table')
parser.add_argument('-mq','--mq', default=30, help='Minimum mapping quality')
parser.add_argument('-bq','--bq', default=20, help='Minimum base quality')
parser.add_argument('-d','--dist', default=5, help='Minimum distance allowed between variants')
parser.add_argument('-ac','--ac', default=3, help='Minimum number of reads supporting a variant')
parser.add_argument('-ns','--num_sites', default=-1, help='Number of sites to be analysed')
parser.add_argument('-str','--strand', default=0, choices=['0','1'], help='Strand filter actiavtion. 0 for desactivating the filter. Default [1]')

EOF

echo Parameters:
echo Tumor bam: "$BAM"
echo Normal bam: "$NBAM"
echo Out vcf file: "$OUT_FILE"
echo Bed file: "$BED"
echo Reference file: "$REF"
echo Beta-binomial parameters table: "$PARAM"
echo Minimum mapping quality: "$MQ"
echo Minimum base quality: "$BQ"
echo Minimum reads supporting variant: "$AC"
echo Minimum distance allowed between variants: "$DIST"
echo Strand filter: "$STRAND"
echo Number of sites analysed: "$NUM_SITES"
echo

#REF=/mnt/share/data/genomes/GRCh37.fa
#BED=/mnt/share/data/enrichment/HNSCC_cfDNA_hotspot_2018_12_20.bed

SCRIPTS=/mnt/users/ahmuyaf1/scripts/Somatic_caller_new.v2/Tools

OUT=$(dirname $OUT_FILE)

mkdir -p $OUT
TEMP=$OUT/Temp
mkdir -p $TEMP

## 1. Split bam
## Split bam file by barcode
#echo python $SCRIPTS/split_bam_step1.py -i $BAM -o $TEMP/DEDUPLICATED 
#python $SCRIPTS/split_bam_step1.py -i $BAM -o $TEMP/DEDUPLICATED 

# If normal bam file is provided, we copy this into 
if [[ ! -z "$NBAM" ]];then
    cp $NBAM $TEMP/NORMAL_DP.bam
    IDn=$(basename $NBAM | sed 's/\.bam$//')
fi

IDt=$(basename $BAM | sed 's/\.bam$//')

## Index new split bam files
#for i in $TEMP/*DP*.bam; do
#    $SAMTOOLS index $i
#done
#
#
### 2. Mpileup files
## Check if bed file is provided
#if [[ ! -z "$BED" ]];then
#    for i in $TEMP/*DP*.bam; do
#        NAME=$(echo $i | sed 's/bam$/pileup/')
#        $SAMTOOLS mpileup -d 9999999 -f $REF \
#            -l $BED \
#            -q $MQ \
#            -Q 1 \
#            -x \
#            $i > $NAME &   
#    done
#    wait
#else
#    for i in $TEMP/*DP*.bam; do
#        NAME=$(echo $i | sed 's/bam$/pileup/')
#        $SAMTOOLS mpileup -d 9999999 -f $REF \
#            -q $MQ \
#            -Q 1 \
#            -x \
#            $i > $NAME &   
#    done
#    wait
#fi
#
### 3. Mpileup to tsv
#for i in $TEMP/*DP*.pileup; do
#    NAME2=$(echo $i | sed 's/pileup$/tsv/')
#    python $SCRIPTS/pileup2tsv_08_01_19.py --pileup $i \
#        -o $NAME2 \
#        --minBQ $BQ \
#        -m 10 &   
#done
#wait

## 4. Paramenters
if [[ -z "$PARAM" ]];then
    Rscript $SCRIPTS/Parameters_step1.r -t1 $TEMP/DEDUPLICATED_DP1.tsv \
        -t2 $TEMP/DEDUPLICATED_DP2.tsv \
        -t3 $TEMP/DEDUPLICATED_DP3.tsv \
        -t4 $TEMP/DEDUPLICATED_DP4.tsv \
        -o $OUT/Parameters.txt
    cp $OUT/Parameters.txt $TEMP/Parameters.txt
else
    cp $PARAM $TEMP/Parameters.txt
fi

## 5. Tsv stats
# If normal bam file is provided
if [[ ! -z "$NBAM" ]];then       
    Rscript $SCRIPTS/Combine_functions_step1.r -t1 $TEMP/DEDUPLICATED_DP1.tsv \
        -t2 $TEMP/DEDUPLICATED_DP2.tsv \
        -t3 $TEMP/DEDUPLICATED_DP3.tsv \
        -t4 $TEMP/DEDUPLICATED_DP4.tsv \
        -n $TEMP/NORMAL_DP.tsv \
        -params $TEMP/Parameters.txt \
        -num $NUM_SITES \
        -o $TEMP/Stats.tsv
else       
    Rscript $SCRIPTS/Combine_functions_step1.r -t1 $TEMP/DEDUPLICATED_DP1.tsv \
        -t2 $TEMP/DEDUPLICATED_DP2.tsv \
        -t3 $TEMP/DEDUPLICATED_DP3.tsv \
        -t4 $TEMP/DEDUPLICATED_DP4.tsv \
        -params $TEMP/Parameters.txt \
        -num $NUM_SITES \
        -o $TEMP/Stats.tsv
fi

## 6. TSV2VCF
# If normal bam file is provided
python $SCRIPTS/TSV2VCF.py -i $TEMP/Stats.tsv \
    -tID $IDt \
    -ref $REF \
    -o $OUT_FILE \
    -cov 10 \
    -ac $AC \
    --strand $STRAND \
    -variant_dist $DIST
#
## TSV
#cp $TEMP/Stats.tsv.tsv $OUT_FILE.tsv
#
## MRD
#cp $TEMP/Stats.tsv.mrd $OUT_FILE.mrd

#rm -r -f $TEMP
