#!/usr/bin/make -f
#
# Integrate annotation, transcriptome-based efforts to identify TRGs
# Daniel Standage, 2014
SHELL=bash
PdomData=/iplant/home/standage/Polistes_dominula

pdom-loci-r1.2-pot-trgs.txt:		pdom-r1.2-no-anm-hits-iloci.txt pdom-tsa-aligned-unmatched-loci.txt
					comm -12 $^ > $@

pdom-r1.2-no-anm-hits-iloci.txt:	pdom-r1.2-no-anm-hits.ids pdom-loci-r1.2-mrnamap.txt
					./selex --out=1 $^ | sort | uniq > $@

pdom-r1.2-no-anm-hits.ids:		
					iget $(PdomData)/r1.2/genome-annotation/$@

pdom-loci-r1.2-mrnamap.txt:		
					iget $(PdomData)/r1.2/interval-loci/$@

pdom-tsa-aligned-unmatched-loci.txt:	pdom-loci-r1.2.gff3 pdom-tsa-aligned-unmatched.gff3
					bedtools intersect -a <(grep $$'\tlocus\t' pdom-loci-r1.2.gff3) -b pdom-tsa-aligned-unmatched.gff3 -wa -u \
					    | perl -ne 'm/(PdomILCr1.2-\d+)/ and print "$$1\n"' \
					    | sort | uniq > $@

pdom-loci-r1.2.gff3:		
					iget $(PdomData)/r1.2/interval-loci/$@

pdom-tsa-aligned-unmatched.gff3:	pdom-tsa-r1.2-unmatched-pep.txt pdom-tsa-masked.gff3
					python align-keep.py $^ > $@

pdom-tsa-masked.gff3:			
					iget $(PdomData)/r1.2/transcript-alignment/$@

pdom-tsa-r1.2.fa:			
					iget $(PdomData)/r1.2/transcript-assembly/$@.gz
					gunzip $@.gz

pdom-tsa-r1.2-unmatched-pep.txt:	
					iget $(PdomData)/r1.2/transcript-assembly/$@

clean:				
				rm -f pdom-loci-r1.2*.txt pdom-tsa*.txt pdom-tsa*.gff3 pdom-r1.2-no-anm-hits* pdom-loci-r1.2.fa pdom-loci-r1.2.gff3 pdom-tsa-r1.2.fa
