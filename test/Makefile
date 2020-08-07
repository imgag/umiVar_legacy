all: clean run check

clean:
	rm -rf cfdna_tmp.* test_out.*

run:
	../cfdna_caller.sh -bam CFDNA_01.bam -b target.bed -r /mnt/share/data/genomes/GRCh37.fa -o test_out.vcf -p p.txt

check:
	diff -u test_out.vcf out.vcf
	diff -u test_out.tsv out.tsv
	diff -u test_out.mrd out.mrd
