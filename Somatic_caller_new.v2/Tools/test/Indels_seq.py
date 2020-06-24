#!/usr/bin/env python

import pybedtools

infile = '/mnt/users/ahmuyaf1/projects/ctDNA/Muyas_et_al/results/test/GRCh37.fa'


def longestRun(s):
    if len(s) == 0: return 0
    runs = ''.join('*' if x == y else ' ' for x,y in zip(s,s[1:]))
    starStrings = runs.split()
    if len(starStrings) == 0: return 1
    return 1 + max(len(stars) for stars in starStrings)
    
def FrequentBase(s):
    L = list(s)
    
    # Get the most common element and the percentage respect the length of the sequence
    MAX = max([L.count(base) for base in set(L)])
    PERC = round(float(MAX)/len(s), 2)
    
    return(PERC)

def Up_Down_sequence(CHROM, START, N, infile):
    # As it works with 0-based coordinates, we must subtract 1 base to our start coordinate
    START = START-1
    
    ## Upstream
    start = str(START)
    end = str(START - N + 1)
    
    ## Getting sequence downstream
    a = pybedtools.BedTool("\t".join([CHROM, end, start]), from_string=True)
    a = a.sequence(fi=infile)
    
    J = open(a.seqfn).read()
    J = J.rstrip('\n')
    
    # Sequence
    SEQ_up = J.split('\n')[1]
    
    ## Downstream
    start = str(START + 1)
    end = str(START + N)
    
    # Getting sequence upstream
    a = pybedtools.BedTool("\t".join([CHROM, start, end]), from_string=True)
    a = a.sequence(fi=infile)
    
    J = open(a.seqfn).read()
    J = J.rstrip('\n')
    SEQ_down = J.split('\n')[1]
    
    # List of sequences
    LIST = [SEQ_up, SEQ_down]
    
    return(LIST)
    

CHROM='chr19'
START=10909154
N=6



Seq_down, Seq_up = Up_Down_sequence(CHROM, START, 6, infile)

# Down
LONG = longestRun(Seq_down)
FREQ = FrequentBase(Seq_down)

print Seq_down, LONG, FREQ


# Up
LONG = longestRun(Seq_up)
FREQ = FrequentBase(Seq_up)

print Seq_up, LONG, FREQ


