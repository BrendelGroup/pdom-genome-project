#!/usr/bin/make -f
# Copyright (c) 2015, Daniel S. Standage

SHELL := bash
SRAURL=ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR

.SECONDARY:

all: pdom-gdnaseq-200bp-1.fq pdom-gdnaseq-200bp-2.fq \
     pdom-gdnaseq-500bp-1.fq pdom-gdnaseq-500bp-2.fq \
     pdom-gdnaseq-1kb-1.fq   pdom-gdnaseq-1kb-2.fq   \
     pdom-gdnaseq-3kb-1.fq   pdom-gdnaseq-3kb-2.fq   \
     pdom-gdnaseq-8kb-1.fq   pdom-gdnaseq-8kb-2.fq

pdom-gdnaseq-200bp-1.fq: SRR1393722.fastq
	@ ln -s $< $@
pdom-gdnaseq-200bp-2.fq: SRR1393725.fastq
	@ ln -s $< $@
pdom-gdnaseq-500bp-1.fq: SRR1409963.fastq
	@ ln -s $< $@
pdom-gdnaseq-500bp-2.fq: SRR1409964.fastq
	@ ln -s $< $@
pdom-gdnaseq-1kb-1.fq: SRR1409967.fastq
	@ ln -s $< $@
pdom-gdnaseq-1kb-2.fq: SRR1409968.fastq
	@ ln -s $< $@
pdom-gdnaseq-3kb-1.fq: SRR1409969.fastq
	@ ln -s $< $@
pdom-gdnaseq-3kb-2.fq: SRR1409970.fastq
	@ ln -s $< $@
pdom-gdnaseq-8kb-1.fq: SRR1409971.fastq
	@ ln -s $< $@
pdom-gdnaseq-8kb-2.fq: SRR1409972.fastq
	@ ln -s $< $@

%.fastq: %.sra
	@ echo [Convert accession $* from .sra to .fastq]
	@ which fastq-dump > /dev/null
	@ fastq-dump $<

%.sra:
	@ echo [Download accession $* from SRA]
	@ PREFIX=$$(echo $* | cut -c 1-6); wget ${SRAURL}/$$PREFIX/$*/$*.sra
