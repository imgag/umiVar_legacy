#!/bin/bash

#REF=/mnt/share/data/genomes/GRCh37.fa
#BED=/mnt/share/data/enrichment/hpHSHNSCC_v1_2016_04_07.bed 
#BAM1=/mnt/users/ahmuyaf1/projects/ctDNA/Haloplex/HNSCC/GS160614_01/GS160614_01_before_dedup_filtered.bam
#BAM2=/mnt/users/ahmuyaf1/projects/ctDNA/Haloplex/HNSCC/GS160614_01/GS160614_01.bam
#OUT=/mnt/users/ahmuyaf1/projects/ctDNA/Haloplex/HNSCC/GS160614_01/STATS
#STEP=3


# Software
SAMTOOLS=/mnt/share/opt/samtools-1.9/samtools
BAMREADCOUNT=/mnt/share/opt/bam-readcount/bin/bam-readcount
ARGPARSE=/mnt/users/ahmuyaf1/programs/argparse/argparse.bash
BEDTOOLS=/mnt/share/opt/bedtools2/bin//bedtools
PAIRTOBED=/mnt/share/opt/bedtools2/bin/pairToBed

# In-house scripts
SCRIPTS=/mnt/users/ahmuyaf1/scripts/Somatic_caller_new.v2/Tools

# Parameters needed for running pipeline
source $ARGPARSE || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-bam','--bam', required=True, help='Deduplicated bam')
parser.add_argument('-r','--ref', required=True, help='Reference genome - fasta')
parser.add_argument('-b','--bed', default='', help='Bed file of the targeted regions')
parser.add_argument('-o', '--out_folder', default='$PWD', help='Out folder (Full path)')
parser.add_argument('-bq', '--bq', default=20, help='Base quality threshold [20]')
parser.add_argument('-mq', '--mq', default=30, help='Mapping quality threshold [30]')
parser.add_argument('-sam', '--sampling', default=1, choices=['0', '1'], help='Sampling: 1 for YES; 0 for NO. Default = 1')
EOF

#echo $READ1, $READ2, $BARCODE1, $BARCODE2, $REF, $BED, $OUT_FOLDER, $SAMPLE_ID
echo Parameters:
echo Required bam: "$BAM"
echo Reference genome: "$REF"
echo Targed bed file: "$BED"
echo Out folder: "$OUT_FOLDER"
echo Base quality: "$BQ"
echo Mapping quality: "$MQ"
#echo Step protocol: ""
echo Sampling: "$SAMPLING"
echo

echo "START"
date

# Sample name and main directory to save data
MAIN_FOLDER=$OUT_FOLDER
TEMP=$MAIN_FOLDER/temp
mkdir -p $MAIN_FOLDER

mkdir -p $MAIN_FOLDER/temp

# If not BED file provided
if [ -z "$BED" ];then
    BED=$TEMP/temp.bed
    $BEDTOOLS bamtobed -i $BAM  | cut -f1,2,3 | sort -k1,1 -k2,2n | $BEDTOOLS merge > $BED
fi


# Sampling bed file. A max of 10000 positions sampled to avoid time
#### Translate barcode groups to the bariables required in bam correct python script
BED2=$TEMP/random.bed
REF_genome=$TEMP/reference.genome

awk -v OFS='\t' '{print $1, $2}' $REF.fai > $REF_genome

# Getting covered size of the original bed file
A=$(awk -F'\t' '{total+=($3-$2)} END {print total}' $BED)

if ([ "$SAMPLING" -eq 1 ] && [ "$A" -gt 100000 ]);then
    # We aprox want to sample X sites
    SITES=100000
    
    # Intersecting new fasta file
    $BEDTOOLS getfasta -fi $REF -bed $BED > ref.intersected.fasta

    # Indexix new fasta
    $SAMTOOLS faidx ref.intersected.fasta

    # Calculating genome file
    awk -v OFS='\t' '{print $1, $2}' ref.intersected.fasta.fai > $REF_genome
    
    # Generating random simulated bed file with X SNPs. At the end, as bedtools uses 0-based bed files, it has to be transformed to 1-based bed file
    $BEDTOOLS random -g $REF_genome -l 1 -n $SITES -seed 1234 | awk -F'\t' -v OFS='\t' '{print $1, $3}' |  sed 's/:\|-/\t/g' | awk -F'\t' -v OFS='\t' '{print $1, $2+$4, $2+$4+1}' | sort -k1,1 -k2,2n | $BEDTOOLS merge -i stdin > $BED2

else
    cp $BED $BED2
fi

##### RUNNING PROGRAM

echo "Step 1 running. Preparing files"
echo
date
echo

# Spliting the condenset bam file into multiples depending on the duplicate count
BAM2=$MAIN_FOLDER/temp/intersected.bam
$SAMTOOLS sort -n $BAM -o $MAIN_FOLDER/temp/sortedn.bam
$PAIRTOBED -abam $MAIN_FOLDER/temp/sortedn.bam -f 0.1 -b $BED2 | $SAMTOOLS sort -o $BAM2 - 

$SAMTOOLS index $BAM2

# Split bam file by barcode
echo python $SCRIPTS/split_bam_step1.py -i $BAM2 -o $MAIN_FOLDER/temp/DEDUPLICATED 
python $SCRIPTS/split_bam_step1.py -i $BAM2 -o $MAIN_FOLDER/temp/DEDUPLICATED 

# Index new split bam files
for i in $MAIN_FOLDER/temp/*DP*.bam; do
    $SAMTOOLS index $i
done

# Getting TSV
for i in $MAIN_FOLDER/temp/*DP*.bam; do
    NAME2=$(echo $i | sed 's/bam$/tsv/')
    $SAMTOOLS mpileup -d 9999999 -f $REF \
        -l $BED2 \
        -q $MQ \
        -Q 1 \
        -x \
        $i | python $SCRIPTS/pileup2tsv_08_01_19.py --pileup - -o $NAME2 --minBQ $BQ -m 10 &
done
wait


echo
echo "Step 2 running. Ploting stats"
echo
date
echo

# Ploting barcode distributions
mkdir -p $MAIN_FOLDER/PLOTS
echo Rscript $SCRIPTS/Barcode_group_distribution.r -bd $MAIN_FOLDER/temp/DEDUPLICATED.txt \
    -out $MAIN_FOLDER/PLOTS
    
Rscript $SCRIPTS/Barcode_group_distribution.r -bd $MAIN_FOLDER/temp/DEDUPLICATED.txt \
    -out $MAIN_FOLDER/PLOTS &

# Error rates
echo Rscript $SCRIPTS/Plot_error_rate_22_11_2018.R -bc1 $MAIN_FOLDER/temp/DEDUPLICATED_DP1.tsv \
    -bc2 $MAIN_FOLDER/temp/DEDUPLICATED_DP2.tsv \
    -bc3 $MAIN_FOLDER/temp/DEDUPLICATED_DP3.tsv \
    -bc4 $MAIN_FOLDER/temp/DEDUPLICATED_DP4.tsv \
    -out $MAIN_FOLDER/PLOTS
Rscript $SCRIPTS/Plot_error_rate_22_11_2018.R -bc1 $MAIN_FOLDER/temp/DEDUPLICATED_DP1.tsv \
    -bc2 $MAIN_FOLDER/temp/DEDUPLICATED_DP2.tsv \
    -bc3 $MAIN_FOLDER/temp/DEDUPLICATED_DP3.tsv \
    -bc4 $MAIN_FOLDER/temp/DEDUPLICATED_DP4.tsv \
    -out $MAIN_FOLDER/PLOTS &

wait 
##rm -f -r $MAIN_FOLDER/temp

echo "END, bye!!!"
date