#!/usr/bin/env python
import argparse
import time
import scipy.stats
import numpy
import pybedtools

numpy.seterr(divide = 'ignore')

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
    START = int(START)
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
    
    
def vcalling(VAR, MM, Seq_Upstream, Seq_Downstream):
    GT, SIG, ALT_COUNT, AB, ALT_COUNT_p, ALT_COUNT_padj, STRAND, FISHER, ALT_COUNT_o, OUT_adj =  VAR.split(":")
    
    # Analysis of flanking sequences
    Seq_Upstream_l = longestRun(Seq_Upstream)
    Seq_Upstream_freq = FrequentBase(Seq_Upstream)
    
    Seq_Downstream_l = longestRun(Seq_Downstream)
    Seq_Downstream_freq = FrequentBase(Seq_Downstream)   
    
    # Reference and alternative for this variant
    ref, alt = GT.split(">")
    
    Alt_fwd, Alt_rev, Ref_fwd, Ref_rev = STRAND.split("-")
    if float(ALT_COUNT) > 0 and MM > 0:
        filter_criteria = []
           
        if (float(ALT_COUNT_padj) > 0.1):
            filter_criteria.append("Error")
            
        if ((Seq_Upstream_l >= 3 or Seq_Upstream_freq >= 0.8) and len(GT) > 3): # Len(GT) > 3 to mark INDELs only. Example: A>TC (len = 4)
            filter_criteria.append("LC_Upstream")
            
        if ((Seq_Downstream_l >= 3 or Seq_Downstream_freq >= 0.8) and len(GT) > 3): # Len(GT) > 3 to mark INDELs only. Example: A>TC (len = 4)
            filter_criteria.append("LC_Downstream")                                            
        
        if (float(AB) < AF):
            filter_criteria.append("Low_AF")
            
        if (min(int(Alt_fwd), int(Alt_rev)) == 0 and min(int(Ref_fwd),int(Ref_rev)) != 0 and strand == 1 and int(ALT_COUNT) > 3):
            filter_criteria.append("Strand_imbalanced")
        
        if (int(DP_HQ) < min_COV):
            filter_criteria.append("Low_Cov")
            
        if (int(ALT_COUNT) < int(min_AC)):
            filter_criteria.append("Low_AC")
            
        if (DIST != "Inf" and int(DIST) < min_DIST):
            filter_criteria.append("Clustered_Variant")
    
        if (int(DP_HQ)/float(DP) < 0.75):
            filter_criteria.append("Low_qual_pos")
    
        if (int(ALT_COUNT_o)/(float(MM)) > 0.3 or float(OUT_adj) < 0.1):
            filter_criteria.append("Variant_contamination")
            
        if (float(FISHER) < 0.01 and strand == 1):
            filter_criteria.append("Fisher_Strand")
                      
        if (len(filter_criteria) == 0):
            FILTER = "PASS"
        else:
            FILTER = ';'.join(filter_criteria)
    else:
        alt = '.'
        FILTER = '.'
    
    TSV = [ALT_COUNT, DP_HQ, AB, REFt, ALT_COUNT_p, ALT_COUNT_padj, STRAND, FISHER, ALT_COUNT_o, OUT_adj, str(Seq_Upstream_freq), str(Seq_Downstream_freq)]
    
    CALL = [ref, alt, TSV, FILTER]
    
    return(CALL)
    



parser = argparse.ArgumentParser(description='Getting barcodes in fastq file and labels')
parser.add_argument('-i', '--infile', type=str, help='Tsv table', required= True)
parser.add_argument('-tID', '--tumorid', type=str, default='Tumor', help='Tumor sample id', required= False)
parser.add_argument('-ref', '--reference', type=str, help='Reference fasta file which table was build', required= True)
parser.add_argument('-o', '--outfile', type=str, help='Vcf output file', required= True)
parser.add_argument('-cov', '--min_COV', type=int, default=10, help='Minimum Coverage', required= False)
parser.add_argument('-ac', '--min_AC', type=int, default=3, help='Minimum reads supporting alternative allele', required= False)
parser.add_argument('-variant_dist', '--min_DIST', type=int, default=20, help='Minimum distance allowed between variants (to avoid clustered errors)', required= False)
parser.add_argument('-str', '--strand', type=int, choices = [0,1], default=1, help='Strand bias test (Fisher test). 0 for turn it off', required= False)
parser.add_argument('-af', '--min_AF', type=float, default=0, help='Minimum allele frequency allowed', required= False)
parser.add_argument('-mrd', '--mrd', type=int, choices = [0,1], default=1, help='1 if MRD must be printed. 0 for not [Default = 1]', required= False)
parser.add_argument('-tmpdir', '--tmpdir', default=None, help='Folder for temp files', required= False)

args = parser.parse_args()


# Check temp folder
if args.tmpdir != None:
    pybedtools.set_tempdir(args.tmpdir)

FILE = args.infile
SAMPLE = args.tumorid
REF = args.reference
out = args.outfile

out2 = out + '.tsv'

min_COV = args.min_COV
min_AC = args.min_AC

min_DIST = args.min_DIST
strand = args.strand
AF = args.min_AF

OUT_vcf = open(out,'w')
OUT_tsv = open(out2,'w')

###### GETTING HEADER
date = time.strftime("%d/%m/%Y")## dd/mm/yyyy format

VCF_format="##fileformat=VCFv4.1"
DATE="##fileDate=%s" % date
source="##source=CRG_UKT_somatic_variant_calling"
reference="##reference=%s" % REF
CONCEPTS="""##INFO=<ID=Variant_Dist,Number=1,Type=Integer,Description="Distance to the closest short variant">
##INFO=<ID=Upstream,Type=String,Description="Upstream sequence (5 nucleotides)">
##INFO=<ID=Downstream,Type=String,Description="Downstream sequence (5 nucleotides)">
##FILTER=<ID=PASS,Description="Passed filter">
##FILTER=<ID=Low_COV,Description="Low coverage">
##FILTER=<ID=Strand_imbalanced,Description="All alternative reads found in only one strand">
##FILTER=<ID=Low_AC,Description="Less than defined minimum of alternative counts">
##FILTER=<ID=Clustered_Variant,Description="Clustered variants">
##FILTER=<ID=LC_Upstream,Description="Low complexity region (5bps) upstream. >= 80% of bases show the same nucleotide or tandem of >= 3 equal nucleotides in a row">
##FILTER=<ID=LC_Downstream,Description="Low complexity region (5bps) downstream. >= 80% of bases show the same nucleotide  or tandem of >= 3 equal nucleotides in a row">
##FILTER=<ID=Error,Description="Alternative counts inside the expected error rate distribution">
##FILTER=<ID=Fisher_Strand,Description="Strand bias based on fisher test">
##FILTER=<ID=Low_qual_pos,Description="Position enriched with too many low quality bases">
##FILTER=<ID=Variant_contamination,Description="Reads supporting other alleles outsite of the error rate distribution">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AC,Number=.,Type=Integer,Description="Allele read counts"
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AB,Number=1,Type=Float,Description="Allele balance">
##FORMAT=<ID=Strand,Number=2,Type=String,Description="Alleles in strands: Alternative forward, Alternative reverseReference forward, Reference reverse">
##FORMAT=<ID=FS,Number=1,Type=Float,Description="Fisher strand test (q-value)">
##FORMAT=<ID=VCB,Number=1,Type=Integer,Description="Variant Count bias, number of other different alternative alleles found">
##FORMAT=<ID=Perror,Number=1,Type=Float,Description="Q-value to belong to the error rate distribution (beta binomial distribution) - After FDR correction">"""
INFILE="##Input file:%s" % FILE
sample_name="##Tumor sample:%s" % SAMPLE
#Parameters="##Parameters for filtering = min_COV: %s; min_AC: %s; min_DIST: %s; strand: %s; min_AF: %s; end_read_filter: %s; somatic_type: %s"  % (min_COV, min_AC, min_DIST, strand, AF, end, args.somatic_type)




## "##" VCF header
OUT_vcf.write(VCF_format + '\n')
OUT_vcf.write(DATE + '\n')
OUT_vcf.write(source + '\n')
OUT_vcf.write(reference + '\n')
OUT_vcf.write(CONCEPTS + '\n')
#OUT_vcf.write(INFILE + '\n')
#OUT_vcf.write(sample_name + '\n')
#OUT_vcf.write(Parameters + '\n')

## TSV header and VCF header
TSV_HEADER = ['CHROM', 'POS', 'REF', 'ALT', 'Upstream_5', 'Downstream_5', 'DP_HQ', 'REFt', 'ALT_COUNT', 'AB', 'P_VAL', 'P_VAL_adj', 'STRAND', 'FISHER', 'ALT_COUNT_o', 'P_VALo_adj', 'LC_Upstream', 'LC_Downstream', 'FILTER']
TSV_HEADER = '\t'.join(TSV_HEADER)
OUT_tsv.write(TSV_HEADER + '\n')

VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', SAMPLE]

MRD = []

with open(FILE) as f1:
    for i in f1:
        line  =  i.rstrip('\n')
        
        if line.startswith('##'):
            continue
        elif line.startswith('CHROM'):
            header_file = line
            form = header_file.split('\t')
            
            if ("CHROM" in form):
                CHROMi=[i for i, x in enumerate(form) if x == "CHROM"][0]
            else:
                print("CHROM not present in vcf")
                break

            if ("POS" in form):
                POSi=[i for i, x in enumerate(form) if x == "POS"][0]
            else:
                print("POS not present in vcf")
                break

            if ("REF" in form):
                REFi=[i for i, x in enumerate(form) if x == "REF"][0]
            else:
                print("REF not present in vcf")
                break

            # Coverage 
            if ("DP" in form):
                DPi=[i for i, x in enumerate(form) if x == "DP"][0]
            else:
                print("DP not present in vcf")
                break
            
            if ("DP_HQ" in form):
                DP_HQi=[i for i, x in enumerate(form) if x == "DP_HQ"][0]
            else:
                print("DP_HQ not present in vcf")
                break
            
            if ("REFf" in form):
                REFfi=[i for i, x in enumerate(form) if x == "REFf"][0]
            else:
                print("REFf not present in vcf")
                break
 
            if ("REFr" in form):
                REFri=[i for i, x in enumerate(form) if x == "REFr"][0]
            else:
                print("REFr not present in vcf")
                break
            
            if ("DIST" in form):
                DISTi=[i for i, x in enumerate(form) if x == "DIST"][0]
            else:
                print("DIST not present in vcf")
                break

            if ("MM" in form):
                MMi=[i for i, x in enumerate(form) if x == "MM"][0]
            else:
                print("MM not present in vcf")
                break

            if ("CALL" in form):
                CALLi=[i for i, x in enumerate(form) if x == "CALL"][0]
            else:
                print("CALL not present in vcf")
                break

            VCF_HEADER_l = "#"+'\t'.join(VCF_HEADER)            
            
            OUT_vcf.write(VCF_HEADER_l + '\n')
            
        else:
            info = line.split("\t")
            ## COMMON VCF columns (FIRST COLUMNS)
            CHROM = info[CHROMi]
            POS = info[POSi]
            ID = '.'
            REF = info[REFi]
            REFf = int(info[REFfi])
            REFr = int(info[REFri])
            REFt = REFf + REFr
            MM = int(info[MMi])
            DP = int(info[DPi])
            DP_HQ = int(info[DP_HQi])
            QUAL = '.'
            DIST = info[DISTi]
            
            # Getting 5 bases up and downstream
            Seq_up, Seq_down = Up_Down_sequence(CHROM, POS, 6, args.reference)
            
            # Read count info
            call = info[CALLi]
            CALL = call.split("|")
            
            #print CHROM, POS
            ref = []
            alt = []
            TSV = []
            filter = []
            for VAR in CALL:
                refVAR, altVAR, tsvVAR, filterVAR = vcalling(VAR, MM, Seq_up, Seq_down)
                
                ## Append variants
                ref.append(refVAR)
                alt.append(altVAR)
                TSV.append(tsvVAR)
                filter.append(filterVAR)
            
            if (len(ref) > 1):

                # Getting ref. For both SNP, Deletion and Insertion, the longer reference will represent the original reference
                REF = max(ref, key=len)
                
                FILTER = []
                ALT = []
                ALT_COUNT_all = []
                sample = []
                
                n = len(ref)
                
                P_val_list = []
                for count in range(0,n):
                    refi = ref[count]
                    if len(refi) != len(REF):
                         # Sample info
                        ALT_COUNT, DP_HQ, AB, REFt, ALT_COUNT_p, ALT_COUNT_padj, STRAND, FISHER, ALT_COUNT_o, OUT_adj, Seq_Upstream_freq, Seq_Downstream_freq = TSV[count]
                        
                        # Getting normalized alternative allele
                        ALTi = alt[count]
                        ALTtemp = list(REF)
                        ALTtemp[0] = ALTi
                        ALTi = ''.join(ALTtemp)
                        ALT.append(ALTi)
                        ALT_COUNT_all.append(ALT_COUNT)
                        
                        if (int(ALT_COUNT) > 0):
                            GT = ALTi
                        else:
                            GT = REF
                        
                        FILTER.append(filter[count])
                        
                        samplei = [GT, str(ALT_COUNT), str(DP_HQ), str(AB), STRAND, str(FISHER), str(ALT_COUNT_o), str(OUT_adj), str(ALT_COUNT_padj)]
                        sample.append(':'.join(samplei))
                        
                        P_val_list.append(float(ALT_COUNT_p))
                        
                    else:
                        
                        # Sample info
                        ALT_COUNT, DP_HQ, AB, REFt, ALT_COUNT_p, ALT_COUNT_padj, STRAND, FISHER, ALT_COUNT_o, OUT_adj, Seq_Upstream_freq, Seq_Downstream_freq = TSV[count]
                        
                        ALTi = alt[count]
                        ALT.append(ALTi)
                        ALT_COUNT_all.append(ALT_COUNT)
                        
                        if (int(ALT_COUNT) > 0):
                            GT = ALTi
                        else:
                            GT = REF
                        
                        FILTER.append(filter[count])
                        
                        samplei = [GT, str(ALT_COUNT), str(DP_HQ), str(AB), STRAND, str(FISHER), str(ALT_COUNT_o), str(OUT_adj), str(ALT_COUNT_padj)]
                        sample.append(':'.join(samplei))
                        
                        P_val_list.append(float(ALT_COUNT_p))

                # Common vcf columns
                ALT = ','.join(ALT)
                COMMON = [CHROM, POS, ID, REF, ALT, QUAL]
                
                # Filter column
                FILTER = '|'.join(FILTER)
                
                # Format column
                FORMAT = ["GT","Alt_Count","DP", "AB", "Strand", "FS", "VCB", "Pvcb", "Perror"]
                
                # Info column
                INFO = ['Variant_dist='+str(DIST), 'Upstream='+str(Seq_up), 'Downstream='+str(Seq_down)]
                
                # Combine p-values of this site
                P_merge = scipy.stats.combine_pvalues(P_val_list)[1]
                MRD.append(float(P_merge))
                
                # VCF variant line
                sample = '|'.join(sample)
                VCF_LINE = ['\t'.join(COMMON), FILTER, ';'.join(INFO), ':'.join(FORMAT), sample]
                VCF_LINE = '\t'.join(VCF_LINE)
                
                # TSV line
                ALT_COUNT_all = ','.join(ALT_COUNT_all)
                TSV_LINE = [str(CHROM), str(POS), str(REF), str(ALT), Seq_up, Seq_down, str(DP_HQ), str(REFt), str(ALT_COUNT_all), str(AB), ALT_COUNT_p, ALT_COUNT_padj, STRAND, FISHER, ALT_COUNT_o, OUT_adj, Seq_Upstream_freq, Seq_Downstream_freq, FILTER]
                TSV_LINE = '\t'.join(TSV_LINE)

                if (int(ALT_COUNT) > 0):
                    OUT_vcf.write(VCF_LINE + '\n')
                    OUT_tsv.write(TSV_LINE + '\n')
                
            else:
                
                # Common vcf columns
                REF = ref[0]
                ALT = alt[0]
                COMMON = [CHROM, POS, ID, REF, ALT, QUAL]
                
                # Filter column
                FILTER = filter[0]
                
                # Info column
                INFO = ['Variant_dist='+str(DIST), 'Upstream='+str(Seq_up), 'Downstream='+str(Seq_down)]
                
                # Format column
                FORMAT = ["GT","Alt_Count","DP", "AB", "Strand", "FS", "VCB", "Pvcb", "Perror"]
                
                # Sample info
                ALT_COUNT, DP_HQ, AB, REFt, ALT_COUNT_p, ALT_COUNT_padj, STRAND, FISHER, ALT_COUNT_o, OUT_adj, Seq_Upstream_freq, Seq_Downstream_freq = TSV[0]
                
                # Append p-val for MRD calculation
                MRD.append(float(ALT_COUNT_p))
                
                # Get genotype
                if (int(ALT_COUNT) > 0):
                    GT = ALT
                else:
                    GT = REF
                
                sample = [GT, str(ALT_COUNT), str(DP_HQ), str(AB), STRAND, str(FISHER), str(ALT_COUNT_o), str(OUT_adj), str(ALT_COUNT_padj)]
                
                # VCF variant line
                VCF_LINE = ['\t'.join(COMMON), FILTER, ';'.join(INFO), ':'.join(FORMAT), ':'.join(sample)]
                VCF_LINE = '\t'.join(VCF_LINE)
                
                # TSV line
                TSV_LINE = [str(CHROM), str(POS), str(ref[0]), str(alt[0]), Seq_up, Seq_down, str(DP_HQ), str(REFt), str(ALT_COUNT), str(AB), ALT_COUNT_p, ALT_COUNT_padj, STRAND, FISHER, ALT_COUNT_o, OUT_adj, Seq_Upstream_freq, Seq_Downstream_freq, FILTER]
                TSV_LINE = '\t'.join(TSV_LINE)
                if (int(ALT_COUNT) > 0):
                    OUT_vcf.write(VCF_LINE + '\n')
                    OUT_tsv.write(TSV_LINE + '\n')

OUT_vcf.close()
OUT_tsv.close()

## MRD printing
if (args.mrd == 1):
    out3 = out + '.mrd'
    OUT_mrd = open(out3,'w')
    
    # Header
    MRD_HEADER = ['#MRD_log10', 'MRD_pval']
    MRD_HEADER = '\t'.join(MRD_HEADER)
    OUT_mrd.write(MRD_HEADER + '\n')
    
    # Getting MRD values
    numpy.seterr(divide = 'ignore') 
    VAL_MRD_temp = scipy.stats.combine_pvalues(MRD)[1]
    VAL_MRD_temp2 = max(1e-20, VAL_MRD_temp)
    VAL_MRD = str(format(VAL_MRD_temp2, "5.2e"))
    VAL_MRD_log = str(round(numpy.log10(float(VAL_MRD))*-1, 4))
    
    # Printing line
    LINE = [VAL_MRD_log, VAL_MRD]
    LINE = '\t'.join(LINE)
    OUT_mrd.write(LINE + '\n')
    
    OUT_mrd.close()
    
exit(0)
