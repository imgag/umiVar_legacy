import pysam
import argparse

### Read parameters
parser = argparse.ArgumentParser(description='Creates a statistic on mismatch and gap frequencies in BAM files.')
parser.add_argument('-i', '--infile', required=True, dest='infile', help='Input BAM file.')
parser.add_argument('-o', '--outprefix', required=True, dest='out', help='Prefix of the out files.')

args = ''
try:
	args = parser.parse_args()
except IOError as io:
	print io
	sys.exit('Error reading parameters.')


### Input BAM
try:
	infile = pysam.Samfile( args.infile, "rb" )
except:
	exit("Cannot open input file.")


### Output BAMs
outfile0=args.out+".txt"
outfile1=args.out+"_DP1.bam"
outfile2=args.out+"_DP2.bam"
outfile3=args.out+"_DP3.bam"
outfile4=args.out+"_DP4.bam"

try:
	Outfile0=open(outfile0,'w')
	Outfile1 = pysam.Samfile(outfile1, mode="wb", template = infile)
	Outfile2 = pysam.Samfile(outfile2, mode="wb", template = infile)
	Outfile3 = pysam.Samfile(outfile3, mode="wb", template = infile)
	Outfile4 = pysam.Samfile(outfile4, mode="wb", template = infile)
except:
	exit("Cannot open output file.")



### Parse BAM file.
while 1:

	# Get read
	try:
		read = infile.next()
	except StopIteration:
		break
	
	try:
		TEST = read.opt("DP")
	except:
		continue
	
	
	# Write duplicates in a single file
	Outfile0.write(str(read.opt("DP"))+"\n")
	
	# Split by the amount of duplicates
	if( read.opt("DP") >= 4 ):
		Outfile4.write(read)
		
	elif( read.opt("DP") == 3 ):
		Outfile3.write(read)
	
	elif( read.opt("DP") == 2 ):
		Outfile2.write(read)

	elif( read.opt("DP") <= 1 ):
		Outfile1.write(read)


### Finish
infile.close()
Outfile0.close()
Outfile1.close()
Outfile2.close()
Outfile3.close()
Outfile4.close()
exit(0)


### pysam documentation
# http://pysam.readthedocs.org/en/latest/
