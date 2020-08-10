#! /usr/bin/python

# IMPORTS
from os import getcwd, path, makedirs
from argparse import ArgumentParser

working_directory = getcwd()

# Create argument parser
argument_parser = ArgumentParser("cfDNA variant caller")

# And add arguments
argument_parser.add_argument('-bam','--bam', required=True, help='Tumor bam file')
argument_parser.add_argument('-nbam','--nbam', default='', required=False, help='Normal bam file')
argument_parser.add_argument('-b','--bed', default='', help='Bed file of the targeted regions. O-based')
argument_parser.add_argument('-r','--ref', required=True, help='Reference genome - fasta')
argument_parser.add_argument('-o', '--out_file', default='$PWD/Sample.vcf', help='Out vcf file')
argument_parser.add_argument('-p', '--param', default='', help='Beta-binomial parameters table')
argument_parser.add_argument('-mq','--mq', default=30, help='Minimum mapping quality')
argument_parser.add_argument('-bq','--bq', default=20, help='Minimum base quality')
argument_parser.add_argument('-d','--dist', default=5, help='Minimum distance allowed between variants')
argument_parser.add_argument('-ac','--ac', default=3, help='Minimum number of reads supporting a variant')
argument_parser.add_argument('-ns','--num_sites', default=-1, help='Number of sites to be analysed')
argument_parser.add_argument('-str','--strand', default=0, choices=['0','1'], help='Strand filter activation. 0 for desactivating the filter. Default [1]')
argument_parser.add_argument('-t', '--temp_dir', default='$PWD', help='Temporary directory')

# Retrieve arguments from stio
arguments = argument_parser.parse_args()

# print(arguments)

tumor_id = arguments.bam.replace(".bam", "")

out_directory = path.split(arguments.out_file)[0]
if not path.exists(out_directory):
    makedirs(out_directory)

