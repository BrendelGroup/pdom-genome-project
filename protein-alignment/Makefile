#!/usr/bin/make -f
SHELL=bash
PdomData=/iplant/home/standage/Polistes_dominula
GTHbin=/usr/local/src/GTH
BSSMDIR=$(GTHbin)/bssm
GTHDATADIR=$(GTHbin)/gthdata

export BSSMDIR := $(BSSMDIR)
export GTHDATADIR := $(GTHDATADIR)

all:						amel-ogs-prot-masked.gth amel-ncbi-prot-masked.gth dmel-flybase-prot-masked.gth

%-prot-masked.gth:				pdom-scaffolds-masked-r1.2.fa %-prot.fa
						gth -genomic pdom-scaffolds-masked-r1.2.fa \
						    -protein $*-prot.fa \
						    -species arabidopsis \
						    -gcmaxgapwidth 20000 \
						    -gcmincoverage 25 \
						    -prhdist 6 \
						    -prminmatchlen 18 \
						    -prseedlength 6 \
						    -o $@ \
						    -force \
						    > $*-prot-masked.log 2>&1

dmel-flybase-prot.fa:				
						curl ftp://ftp.flybase.net/releases/FB2014_01/dmel_r5.55/fasta/dmel-all-translation-r5.55.fasta.gz \
						    | zcat > $@

amel-ncbi-prot.fa:				
						curl ftp://ftp.ncbi.nih.gov/genomes/Apis_mellifera/protein/protein.fa.gz \
						    | zcat | perl -ne 's/>gi\|(\d+)\|ref\|([^|]+)\|\S*/>$1 $2/; print' > $@

amel-ogs-prot.fa:				
						curl http://hymenopteragenome.org/beebase/sites/hymenopteragenome.org.beebase/files/data/consortium_data/amel_OGSv3.2_pep.fa.gz \
						    | zcat | sed 's/gnl|Amel_4.5|//g' > $@

pdom-scaffolds-masked-r1.2.fa:		
						iget $(PdomDATA)/r1.2/genome-assembly/$@.gz
						gunzip $@.gz
