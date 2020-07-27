#!/usr/bin/env python3

import argparse
import collections

import pysam

parser = argparse.ArgumentParser(
    description='Split BAM file by DP tag and output DP frequency table.')
parser.add_argument('-i', '--infile', required=True,
                    dest='infile', help='Input BAM file.')
parser.add_argument('-o', '--outprefix', required=True,
                    dest='out', help='Prefix of output files.')
parser.add_argument('-f', '--freq', required=False,
                    dest='freq', help='Frequency table output file.')

args = parser.parse_args()

# output BAM files
out1_fp = args.out+'_DP1.bam'
out2_fp = args.out+'_DP2.bam'
out3_fp = args.out+'_DP3.bam'
out4_fp = args.out+'_DP4.bam'

# frequency table
dp_counter = collections.Counter()

# input BAM file
with pysam.AlignmentFile(args.infile, 'rb') as infile, \
        pysam.AlignmentFile(out1_fp, 'wb', template=infile) as out1, \
        pysam.AlignmentFile(out2_fp, 'wb', template=infile) as out2, \
        pysam.AlignmentFile(out3_fp, 'wb', template=infile) as out3, \
        pysam.AlignmentFile(out4_fp, 'wb', template=infile) as out4:

    # iterate over BAM file
    for aln in infile:
        # read DP tag
        try:
            dp = aln.get_tag('DP')
        except KeyError:
            continue

        # add to counter
        dp_counter.update([dp])

        # split by DP
        if (dp >= 4):
            out4.write(aln)
        elif (dp == 3):
            out3.write(aln)
        elif (dp == 2):
            out2.write(aln)
        elif(dp <= 1):
            out1.write(aln)

if args.freq:
    with open(args.freq, 'w') as f:
        f.write('dp\tcount\n')
        for (dp, count) in sorted(dp_counter.items()):
            f.write('{}\t{}\n'.format(dp, count))
