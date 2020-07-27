#!/bin/bash

# cfDNA variant caller

set -e -u -o pipefail

# script directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# arguments
ARGPARSE_DESCRIPTION="cfDNA variant caller"
ARGPARSE=$DIR/argparse.bash
source $ARGPARSE || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-bam','--bam', required=True, help='Tumor bam file')
parser.add_argument('-nbam','--nbam', default='', required=False, help='Normal bam file')
parser.add_argument('-b','--bed', default='', help='Bed file of the targeted regions. O-based')
parser.add_argument('-r','--ref', required=True, help='Reference genome - fasta')
parser.add_argument('-o', '--out_file', default='$PWD/Sample.vcf', help='Out vcf file')
parser.add_argument('-p', '--param', default='', help='Beta-binomial parameters table')
parser.add_argument('-mq','--mq', default=30, help='Minimum mapping quality')
parser.add_argument('-bq','--bq', default=20, help='Minimum base quality')
parser.add_argument('-d','--dist', default=5, help='Minimum distance allowed between variants')
parser.add_argument('-ac','--ac', default=3, help='Minimum number of reads supporting a variant')
parser.add_argument('-ns','--num_sites', default=-1, help='Number of sites to be analysed')
parser.add_argument('-str','--strand', default=0, choices=['0','1'], help='Strand filter activation. 0 for desactivating the filter. Default [1]')
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

tumor_id=$(basename $BAM .bam)

out_dir=$(dirname $OUT_FILE)
mkdir -p $out_dir
temp=$(mktemp -p $out_dir -d cfdna_tmp.XXXX)

python3 $DIR/split_bam.py --infile $BAM --outprefix "${temp}/dedup"
for f in ${temp}/dedup*.bam; do
    samtools index $f
    f_tsv="${f%.bam}.tsv"
    # TODO handle -l BED optional
    samtools mpileup \
        -d 9999999 -f $REF -l $BED -Q 1 -x $f 2> /dev/null | \
        python3 $DIR/pileup2tsv.py --pileup - --variantf $f_tsv --minBQ $BQ --mindepth 10
    rm $f ${f}.bai
done

if [[ -z "$PARAM" ]]; then
    # TODO check for sufficient bases in target file
    f_params="${temp}/parameters.txt"

    Rscript $DIR/Parameters_step1.r \
        -t1 "${temp}/dedup_DP1.tsv" \
        -t2 "${temp}/dedup_DP2.tsv" \
        -t3 "${temp}/dedup_DP3.tsv" \
        -t4 "${temp}/dedup_DP4.tsv" \
        -o $f_params
else
    f_params=$PARAM
fi

Rscript $DIR/Combine_functions_step1.r \
    -t1 "${temp}/dedup_DP1.tsv" \
    -t2 "${temp}/dedup_DP2.tsv" \
    -t3 "${temp}/dedup_DP3.tsv" \
    -t4 "${temp}/dedup_DP4.tsv" \
    -params $f_params \
    -num $NUM_SITES \
    -o "${temp}/stats.tsv"

python3 $DIR/TSV2VCF.py \
    -i "${temp}/stats.tsv" \
    -tID $tumor_id \
    -ref $REF \
    -o $OUT_FILE \
    -cov 10 \
    -ac $AC \
    --strand $STRAND \
    -variant_dist $DIST \
    -tmpdir $temp
