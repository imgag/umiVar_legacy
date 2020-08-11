all: clean run new check

clean:
	rm -rf cfdna_tmp.* test_out.*

run:
	../cfdna_caller.sh -bam CFDNA_01.bam -b target.bed -r ./GRCh37.fa -o test_out.vcf -p p.txt -t .

new:
	python3 ../cfdna_caller.py -bam CFDNA_01.bam -b target.bed -r GRCh37.fa -o test_out_new.vcf -p p.txt -t .

check:
	-diff -u test_out.vcf out.vcf
	-diff -u test_out.tsv out.tsv
	-diff -u test_out.mrd out.mrd
	-diff -u test_out_new.vcf out.vcf
	-diff -u test_out_new.tsv out.tsv
	-diff -u test_out_new.mrd out.mrd
