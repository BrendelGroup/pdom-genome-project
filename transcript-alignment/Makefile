#!/usr/bin/make -f
SHELL=bash
GSQDir=/usr/local/src/GENESEQER
PdomData=/iplant/home/standage/Polistes_dominula

all:					pdom-tsa-r1.2-masked.gsq pcan-tsa-r1.2-masked.gsq pmet-tsa-r1.2-masked.gsq \
					pdom-tsa-r1.2-masked.gff3 pcan-tsa-r1.2-masked.gff3 pmet-tsa-r1.2-masked.gff3 \
					pdom-tsa-r1.2.fa pcan-tsa.fa pmet-tsa-r1.2.fa

%-tsa-r1.2-masked.gff3:			%-tsa-r1.2-masked.gsq
					./gsq2makergff3.py < $^ > $@

pcan-tsa-r1.2-masked.gsq:		pdom-scaffolds-masked-r1.2.fa pcan-tsa.fa
					MakeArray pcan-tsa.fa
					GeneSeqerL -s Arabidopsis \
					           -L pdom-scaffolds-masked-r1.2.fa \
					           -d pcan-tsa.fa \
					           -O $@ \
					           -p $GSQDIR/data/prmfile \
					           -x 16 -y 24 -z 48 -w 0.8 -m 1000000 \
					           > pcan-tsa-masked-gsq.log 2>&1

pmet-tsa-r1.2-masked.gsq:		pdom-scaffolds-masked-r1.2.fa pmet-tsa-r1.2.fa
					MakeArray pmet-tsa-r1.2.fa
					GeneSeqerL -s Arabidopsis \
					           -L pdom-scaffolds-masked-r1.2.fa \
					           -d pmet-tsa-r1.2.fa \
					           -O $@ \
					           -p $GSQDIR/data/prmfile \
					           -x 16 -y 24 -z 48 -w 0.8 -m 1000000 \
					           > pmet-tsa-r1.2-masked-gsq.log 2>&1

pdom-tsa-r1.2-masked.gsq:		pdom-scaffolds-masked-r1.2.fa pdom-tsa-r1.2.fa
					MakeArray pdom-tsa-r1.2.fa
					GeneSeqerL -s Arabidopsis \
					           -L pdom-scaffolds-masked-r1.2.fa \
					           -D pdom-tsa-r1.2.fa \
					           -O $@ \
					           -p $GSQDIR/data/prmfileHQ \
					           -x 30 -y 45 -z 60 -w 0.8 -m 1000000 \
					           > pdom-tsa-r1.2-masked-gsq.log 2>&1

%.fa:			
					iget $(PdomData)/r1.2/transcript-assembly/$@.gz
					gunzip $@.gz

pdom-scaffolds-masked-r1.2.fa:		
					iget $(PdomData)/r1.2/genome-assembly/$@.gz
					gunzip $@.gz
