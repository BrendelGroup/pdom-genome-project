#!/usr/bin/make -f
SHELL=bash
NumThreads=1
TrimJar=/usr/local/src/Trimmomatic-0.22/trimmomatic-0.22.jar
CWD=$(shell pwd)
DataDir=$(CWD)/assembly/Polistes_dominula/data-trim

RawFastqs=pdom-gdnaseq-200bp-1.fq.gz \
          pdom-gdnaseq-200bp-2.fq.gz \
          pdom-gdnaseq-500bp-1.fq.gz \
          pdom-gdnaseq-500bp-2.fq.gz \
          pdom-gdnaseq-1kb-1.fq.gz \
          pdom-gdnaseq-1kb-2.fq.gz \
          pdom-gdnaseq-3kb-1.fq.gz \
          pdom-gdnaseq-3kb-2.fq.gz \
          pdom-gdnaseq-8kb-1.fq.gz \
          pdom-gdnaseq-8kb-2.fq.gz

AllFastqs=$(DataDir)/pdom-gdnaseq-200bp-trim-1.fq \
          $(DataDir)/pdom-gdnaseq-500bp-trim-1.fq \
          $(DataDir)/pdom-gdnaseq-1kb-trim-1.fq \
          $(DataDir)/pdom-gdnaseq-3kb-trim-1.fq \
          $(DataDir)/pdom-gdnaseq-8kb-trim-1.fq

all:						pdom-scaffolds-unmasked-r1.2.fa.gz $(RawFastqs)

pdom-scaffolds-unmasked-r1.2.fa.gz:		AssemblyComplete scaff-ids.pl
						./scaff-ids.pl PdomSCFr1.2- \
						    < $(DataDir)/run01/ASSEMBLIES/test/final.assembly.fasta \
						    > pdom-scaffolds-unmasked-r1.2.fa

AssemblyComplete:				PrepComplete in_libs.csv in_groups.csv
						RunAllPathsLG PRE=$(CWD)/assembly \
						    REFERENCE_NAME=Polistes_dominula \
						    DATA_SUBDIR=data-trim \
						    RUN=run01 \
						    TARGETS=standard \
						    THREADS=32
						date > $@

PrepComplete:					$(AllFastqs) in_libs.csv in_groups.csv
						PrepareAllPathsInputs.pl \
						    DATA_DIR=$(DataDir) \
						    PLOIDY=2
						date > $@

$(DataDir)/pdom-gdnaseq-%-trim-1.fq:		pdom-gdnaseq-%-1.fq.gz pdom-gdnaseq-%-2.fq.gz
						./run-trim.sh $* $(TrimJar) $(NumThreads)
						mkdir -p $(DataDir)
						mv pdom-gdnaseq-$*-trim-?.fq $(DataDir)

pdom-gdnaseq-%.fq.gz:				
						which iget
						iget /iplant/home/standage/Polistes_dominula/sequence/genome/$@

clean:						
						rm -rf assembly pdom-gdnaseq*.fq*
