#!/usr/bin/make -f
# Copyright (c) 2015, Daniel S. Standage

SHELL=bash
NumThreads=1
FastqFiles=pdom-gdnaseq-1kb-1.fq   \
           pdom-gdnaseq-1kb-2.fq   \
           pdom-gdnaseq-200bp-1.fq \
           pdom-gdnaseq-200bp-2.fq \
           pdom-gdnaseq-3kb-1.fq   \
           pdom-gdnaseq-3kb-2.fq   \
           pdom-gdnaseq-500bp-1.fq \
           pdom-gdnaseq-500bp-2.fq \
           pdom-gdnaseq-8kb-1.fq   \
           pdom-gdnaseq-8kb-2.fq
AllJf=pdom-17mers.jf pdom-21mers.jf pdom-25mers.jf pdom-29mers.jf
AllHistos=pdom-17mers.hist pdom-21mers.hist pdom-25mers.hist pdom-29mers.hist

all:			pdom-size-kmers.png $(FastqFiles) $(AllJf) $(AllHistos)
			

pdom-size-kmers.png:	$(AllHistos)
			which Rscript
			./size-coverage-estimate.R

pdom-%mers.hist:	pdom-%mers.jf
			which jellyfish
			jellyfish histo $< > $@

pdom-%mers.jf:		$(FastqFiles)
			which jellyfish
			jellyfish count -m $* -s 100M -t $(NumThreads) -C -o $@ $(FastqFiles)

clean:			
			rm -rf  $(AllJf)

clean-all:		clean
			rm -rf $(AllHistos)

