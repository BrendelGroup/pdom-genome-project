#!/usr/bin/make -f
SHELL=bash
PdomData=/iplant/home/standage/Polistes_dominula

all:					pdom-annot-r1.2-iloci.gff3 pdom-annot-r1.2-iloci-sanstrna.gff3 pdom-loci-r1.2.fa pdom-loci-r1.2.gff3 amel-pdom-iloci-gc.png

amel-pdom-iloci-gc.png:			hymenoptera-ilocus-stats-2014.tsv
					./locus-composition.R

pdom-annot-r1.2-iloci-sanstrna.gff3:	pdom-loci-r1.2.gff3
					perl glb2lcl.pl < $^ \
					    | grep -v -e trnascan -e PdomTRNA \
					    | gt gff3 -retainids -sort -tidy -o $@

pdom-annot-r1.2-iloci.gff3:		pdom-loci-r1.2.gff3
					perl glb2lcl.pl < $^ | gt gff3 -retainids -sort -tidy -o $@

pdom-loci-r1.2.fa:			pdom-loci-r1.2.gff3 pdom-scaffolds-unmasked-r1.2.fa
					xtractore --type=locus --outfile=$@ $^

pdom-loci-r1.2.gff3:			pdom-annot-r1.2.gff3
					locuspocus --intloci --skipends \
					           --delta=500 --verbose \
					           --outfile=$@ \
					           $^

pdom-scaffolds-unmasked-r1.2.fa:	
					iget $(PdomData)/genome-assembly/$@.gz
					gunzip $@.gz

pdom-annot-r1.2.gff3:			
					iget $(PdomData)/genome-annotation/$@

clean:					
					rm -rf pdom-annot-r1.2-iloci.gff3 pdom-annot-r1.2-iloci-sanstrna.gff3 pdom-loci-r1.2.fa pdom-loci-r1.2.gff3 pdom-annot-r1.2.gff3 pdom-scaffolds-unmasked-r1.2.fa
