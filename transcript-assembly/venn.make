#!/usr/bin/make -f
#
# Generating a Venn diagram of sequence overlap.
# Based on Volker's script xme-VENN
# Volker Brendel,  2013
# Daniel Standage, 2014
SHELL=bash
BLASTTHREADS=1

all:		pdom-venn.png pdom-tsa-r1.2-unmatched-pep.fa

clean:		
		rm -rf blastp-P* vID-P* ID-P* msb-P* *.log *MRNA* *PRT* Pc* Pd* Pm*

#-------------------------------------------------------------------------------
# Download transcript assemblies
#-------------------------------------------------------------------------------
IRODS_DIR=/iplant/home/standage/Polistes_dominula/r1.2/transcript-assembly

Pc:		
		test -f Pc-MRNAMARKUP.tar.gz || iget $(IRODS_DIR)/Pc-MRNAMARKUP.tar.gz
		test -d Pc-MRNAMARKUP || tar xzf Pc-MRNAMARKUP.tar.gz
		cp Pc-MRNAMARKUP/Pc-out/unmatched-PcMRNA.fas Pc

Pd:		
		test -f Pd-MRNAMARKUP.tar.gz || iget $(IRODS_DIR)/Pd-MRNAMARKUP.tar.gz
		test -d Pd-MRNAMARKUP || tar xzf Pd-MRNAMARKUP.tar.gz
		cp Pd-MRNAMARKUP/unmatched-0.1.3.Trinity.fasta.fas Pd

Pm:		
		test -f Pm-MRNAMARKUP.tar.gz || iget $(IRODS_DIR)/Pm-MRNAMARKUP.tar.gz
		test -d Pm-MRNAMARKUP || tar xzf Pm-MRNAMARKUP.tar.gz
		cp Pm-MRNAMARKUP/Pm-out/unmatched-PmTSA.fas Pm

#-------------------------------------------------------------------------------
# Translate TSAs, index protein databases for blast search
#-------------------------------------------------------------------------------
allPRT:		PcPRT PdPRT PmPRT
allID:		ID-PcPRT ID-PdPRT ID-PmPRT

%PRT:		%
		dnatopro -F 80 -l $* -o $@ 2> dnatoprot-$*.log
		./xedf$* $@
		makeblastdb -in $@ -dbtype prot -parse_seqids > makeblastdb-$*.log 2>&1

ID-%PRT:	%PRT
		egrep '^>' $< | cut -d ' ' -f 1 | cut -d '|' -f 2 > $@

#-------------------------------------------------------------------------------
# Blast searches and post-processing
#-------------------------------------------------------------------------------
allBLASTP:		blastp-Pc-vs-Pd blastp-Pc-vs-Pm blastp-Pd-vs-Pc blastp-Pd-vs-Pm blastp-Pm-vs-Pc blastp-Pm-vs-Pd
allMSB:			msb-Pc-vs-Pd msb-Pc-vs-Pm msb-Pd-vs-Pc msb-Pd-vs-Pm msb-Pm-vs-Pc msb-Pm-vs-Pd

blastp-Pc-vs-Pd:	PcPRT PdPRT
			blastp -num_threads $(BLASTTHREADS) -query PcPRT -db PdPRT -evalue 1e-20 -num_descriptions 10 -num_alignments 10 -out $@
blastp-Pc-vs-Pm:	PcPRT PmPRT
			blastp -num_threads $(BLASTTHREADS) -query PcPRT -db PmPRT -evalue 1e-20 -num_descriptions 10 -num_alignments 10 -out $@
blastp-Pd-vs-Pc:	PdPRT PcPRT
			blastp -num_threads $(BLASTTHREADS) -query PdPRT -db PcPRT -evalue 1e-20 -num_descriptions 10 -num_alignments 10 -out $@
blastp-Pd-vs-Pm:	PdPRT PmPRT
			blastp -num_threads $(BLASTTHREADS) -query PdPRT -db PmPRT -evalue 1e-20 -num_descriptions 10 -num_alignments 10 -out $@
blastp-Pm-vs-Pc:	PmPRT PcPRT
			blastp -num_threads $(BLASTTHREADS) -query PmPRT -db PcPRT -evalue 1e-20 -num_descriptions 10 -num_alignments 10 -out $@
blastp-Pm-vs-Pd:	PmPRT PdPRT
			blastp -num_threads $(BLASTTHREADS) -query PmPRT -db PdPRT -evalue 1e-20 -num_descriptions 10 -num_alignments 10 -out $@

msb-%:			blastp-%
			./run-msb.sh $*


#-------------------------------------------------------------------------------
# IDs of proteins in different sectors of Venn diagram
#-------------------------------------------------------------------------------
allVID:		vID-Pc vID-Pd vID-Pm

vID-Pc:		allMSB ID-PcPRT
		./run-diff.sh Pc Pd Pm
vID-Pd:		allMSB ID-PdPRT
		./run-diff.sh Pd Pc Pm
vID-Pm:		allMSB ID-PmPRT
		./run-diff.sh Pm Pc Pd

pdom-tsa-r1.2.fa:	
			iget /iplant/home/standage/Polistes_dominula/r1.2/transcript-assembly/$@.gz                                
			gunzip $@.gz

pdom-tsa-r1.2-unmatched-pep.txt:	vID-Pd PdPRT pdom-tsa-r1.2.fa
					./id-resolve.py vID-PdPcPm \
					    PdPRT \
					    pdom-tsa-r1.2.fa \
					    > $@

pdom-tsa-r1.2-unmatched-pep.fa:		pdom-tsa-r1.2-unmatched-pep.txt pdom-tsa-r1.2.fa
					./select-seq.pl $^ > $@

#-------------------------------------------------------------------------------
# Make the Venn diagram graphic
#-------------------------------------------------------------------------------

pdom-venn.png:	allVID make-rvenn.py
		./make-rvenn.py | Rscript -
