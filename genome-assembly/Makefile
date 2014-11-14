#!/usr/bin/make -f
SHELL=bash
NumThreads=1
TrimJar=/usr/local/src/Trimmomatic-0.22/trimmomatic-0.22.jar
DataDir=assembly/Polistes_dominula/data-trim

pdom-scaffolds-unmasked-r1.2.fa.gz:		AssemblyComplete scaff-ids.pl
						./scaff-ids.pl PdomSCFr1.2- \
						    < $(DataDir)/run01/ASSEMBLIES/test/final.assembly.fasta \
						    > pdom-scaffolds-unmasked-r1.2.fa

AssemblyComplete:				PrepComplete in_libs.csv in_groups.csv
						RunAllPathsLG PRE=assembly \
						    REFERENCE_NAME=Polistes_dominula \
						    DATA_SUBDIR=data-trim \
						    RUN=run01 \
						    TARGETS=standard
						date > $@

PrepComplete:					$(AllFastqs) in_libs.csv in_groups.csv
						PrepareAllPathsInputs.pl \
						    DATA_DIR=$(DataDir) \
						    PLOIDY=2
						date > $@

$(DataDir)/pdom-gdnaseq-%-trim-paired-1.fq:	pdom-gdnaseq-%-1.fq.gz pdom-gdnaseq-%-2.fq.gz
						./run-trim.sh $* $(TrimJar) $(NumThreads)
						mkdir -p $(DataDir)
						pdom-gdnaseq-$*-trim-paired-*.fq $(DataDir)

pdom-gdnaseq-%.fq.gz:				
						which iget
						iget /iplant/home/standage/Polistes_dominula/sequence/genome/$@.gz

clean:						
						rm -rf assembly pdom-gdnaseq*.fq*
