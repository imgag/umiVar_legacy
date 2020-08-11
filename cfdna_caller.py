#! /usr/bin/python

# LIBRARY IMPORTS
import collections
import glob
import os
import argparse
import tempfile
import pysam
import subprocess


# Create argument parser
argument_parser = argparse.ArgumentParser("cfDNA variant caller")

# And add arguments
argument_parser.add_argument('-bam','--bam', required=True, help='Tumor bam file')
argument_parser.add_argument('-nbam','--nbam', default='', required=False, help='Normal bam file')
argument_parser.add_argument('-b','--bed', default='', help='Bed file of the targeted regions. O-based')
argument_parser.add_argument('-r','--ref', required=True, help='Reference genome - fasta')
argument_parser.add_argument('-o', '--out_file', default='./Sample.vcf', help='Out vcf file')
argument_parser.add_argument('-p', '--param', default='', help='Beta-binomial parameters table')
argument_parser.add_argument('-mq','--mq', default=30, help='Minimum mapping quality')
argument_parser.add_argument('-bq','--bq', default=20, help='Minimum base quality')
argument_parser.add_argument('-d','--dist', default=5, help='Minimum distance allowed between variants')
argument_parser.add_argument('-ac','--ac', default=3, help='Minimum number of reads supporting a variant')
argument_parser.add_argument('-ns','--num_sites', default=-1, help='Number of sites to be analysed')
argument_parser.add_argument('-str','--strand', default=0, choices=['0','1'], help='Strand filter activation. 0 for desactivating the filter. Default [1]')
argument_parser.add_argument('-t', '--temp_dir', default='.', help='Temporary directory')

# Retrieve arguments from stio
arguments = argument_parser.parse_args()

tumor_id = arguments.bam.replace('.bam', '')

# Get the current working directory
working_directory = "../"

# Retrieve and create the given output directory
out_directory = os.path.dirname(arguments.out_file)
if (out_directory == ''):
    print("out directory not given, using working directory")
    out_directory = working_directory

# If the output directory does not exist, create it with all subdirectories
if not os.path.exists(out_directory):
    os.makedirs(out_directory)

# Create a temporary working directory
temp_directory_path = arguments.temp_dir
temp_directory = tempfile.mkdtemp(dir=arguments.temp_dir, prefix='cfdna_tmp.')


def split_bam(infile, outprefix, frequency):
    '''Split BAM files...'''

    # output BAM files
    out1_fp = outprefix + '_DP1.bam'
    out2_fp = outprefix + '_DP2.bam'
    out3_fp = outprefix + '_DP3.bam'
    out4_fp = outprefix + '_DP4.bam'

    dp_counter = collections.Counter()

    # input BAM file
    with pysam.AlignmentFile(infile, 'rb') as infile, \
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

    if frequency:
        with open(frequency, 'w') as f:
            f.write('dp\tcount\n')
            for (dp, count) in sorted(dp_counter.items()):
                f.write('{}\t{}\n'.format(dp, count))


out_prefix = temp_directory + '/dedup'
frequency = temp_directory + '/dp_freq.tsv'

# Split BAM file into 4 different BAM files
split_bam(infile=arguments.bam, outprefix=out_prefix, frequency=frequency)

# Go over all of the split BAM files
for file in sorted(glob.glob(temp_directory + '/dedup*.bam')):
    print(file)

    # Index them
    try:
        subprocess.run("samtools index " + file, shell=True)
    except subprocess.CalledProcessError as error:
        print(error)

    # And replace the BAM file extension with tsv
    f_tsv = file.replace('.bam', '.tsv')

    # TODO handle -l BED as optional
    mpileup_command = 'samtools mpileup' + \
                      ' -d 0 ' + \
                      ' -f ' + arguments.ref + \
                      ' -l ' + arguments.bed + \
                      ' -Q 1 ' + \
                      ' -x ' + file + \
                      ' | ' + \
                      ' python3 ' + working_directory + '/pileup2tsv.py' + \
                      ' --pileup - ' + \
                      ' --variantf ' + f_tsv + \
                      ' --minBQ ' + str(arguments.bq) + \
                      ' --mindepth 10'

    try:
        subprocess.run(mpileup_command, shell=True)
    except subprocess.CalledProcessError as error:
        print(error)

    # Remove the temporary bam and bam.bai files
    os.remove(file)
    os.remove(file + '.bai')

f_params = None

# If no parameter files is given, estimate the optimal parameters
if arguments.param == '':
    f_params = temp_directory + '/parameters.txt'

    parameters_step1_command = 'Rscript ' + working_directory + '/Parameters_step1.R' + \
                              ' -t1 ' + temp_directory + '/dedup_DP1.tsv' + \
                              ' -t2 ' + temp_directory + '/dedup_DP2.tsv' + \
                              ' -t3 ' + temp_directory + '/dedup_DP3.tsv' + \
                              ' -t4 ' + temp_directory + '/dedup_DP4.tsv' + \
                              ' -o '  + f_params

    try:
        subprocess.run(parameters_step1_command, shell=True)
    except subprocess.CalledProcessError as error:
        print(error)
# Else set use the given parameters
else:
    f_params = arguments.param

print("Worked 1")

combine_functions_step1_command = 'Rscript ' + working_directory + '/Combine_functions_step1.R' + \
                                  ' -t1 ' + temp_directory + '/dedup_DP1.tsv' + \
                                  ' -t2 ' + temp_directory + '/dedup_DP2.tsv' + \
                                  ' -t3 ' + temp_directory + '/dedup_DP3.tsv' + \
                                  ' -t4 ' + temp_directory + '/dedup_DP4.tsv' + \
                                  ' -params ' + f_params + \
                                  ' -num ' + str(arguments.num_sites) + \
                                  ' -o ' + temp_directory + '/stats.tsv'

try:
    subprocess.run(combine_functions_step1_command, shell=True)
except subprocess.CalledProcessError as error:
    print(error)

print("Worked 2")

tsv2vcf_command = 'python3 ' + working_directory + '/TSV2VCF.py' + \
                  ' -i ' + temp_directory + '/stats.tsv' + \
                  ' -tID ' + tumor_id + \
                  ' -ref ' + arguments.ref + \
                  ' -o ' + arguments.out_file + \
                  ' -cov 10' + \
                  ' -ac ' + str(arguments.ac) + \
                  ' --strand ' + str(arguments.strand) + \
                  ' -variant_dist ' + str(arguments.dist) + \
                  ' -tmpdir ' + temp_directory

try:
    subprocess.run(tsv2vcf_command, shell=True)
except subprocess.CalledProcessError as error:
    print(error)

print('worked 3')

plot_error_rate_command = 'Rscript ' + working_directory + '/Plot_error_rate_22_11_2018.R' + \
                          ' -bc1 ' + temp_directory + '/dedup_DP1.tsv' + \
                          ' -bc2 ' + temp_directory + '/dedup_DP2.tsv' + \
                          ' -bc3 ' + temp_directory + '/dedup_DP3.tsv' + \
                          ' -bc4 ' + temp_directory + '/dedup_DP4.tsv' + \
                          ' -out ' + temp_directory

try:
    subprocess.run(plot_error_rate_command, shell=True)
except subprocess.CalledProcessError as error:
    print(error)
