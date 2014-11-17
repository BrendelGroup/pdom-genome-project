#!/usr/bin/make -f
SHELL=bash
GSQDir=/usr/local/src/GENESEQER
PdomData=/iplant/home/standage/Polistes_dominula

all:					pdom-tsa-r1.2-masked.gsq pcan-tsa-r1.2-masked.gsq pmet-tsa-r1.2-masked.gsq

%-tsa-r1.2-masked.gsq:			pdom-scaffolds-masked-r1.2.fa %-tsa-r1.2.fa
					MakeArray $*-tsa-r1.2.fa
					GeneSeqerL -s Arabidopsis \
					           -L pdom-scaffolds-masked-r1.2.fa \
					           -d $*-tsa-r1.2.fa \
					           -O $@ \
					           -p $GSQDIR/data/prmfile \
					           -x 16 -y 24 -z 48 -w 0.8 -m 1000000 \
					           > $*-tsa-r1.2-masked-gsq.log 2>&1

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
