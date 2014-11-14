#!/usr/bin/make -f
SHELL=bash
NumThreads=1
TrimJar=/usr/local/src/Trimmomatic-0.22/trimmomatic-0.22.jar
PdomData=/iplant/home/standage/Polistes_dominula

LeftReads=pdom-rnaseq-q1-trim-1.fq \
          pdom-rnaseq-q2-trim-1.fq \
          pdom-rnaseq-q3-trim-1.fq \
          pdom-rnaseq-q4-trim-1.fq \
          pdom-rnaseq-q5-trim-1.fq \
          pdom-rnaseq-q6-trim-1.fq \
          pdom-rnaseq-w1-trim-1.fq \
          pdom-rnaseq-w2-trim-1.fq \
          pdom-rnaseq-w3-trim-1.fq \
          pdom-rnaseq-w4-trim-1.fq \
          pdom-rnaseq-w5-trim-1.fq \
          pdom-rnaseq-w6-trim-1.fq 

RawFastqs=pdom-rnaseq-q1-1.fq.gz pdom-rnaseq-q1-2.fq.gz \
          pdom-rnaseq-q2-1.fq.gz pdom-rnaseq-q2-2.fq.gz \
          pdom-rnaseq-q3-1.fq.gz pdom-rnaseq-q3-2.fq.gz \
          pdom-rnaseq-q4-1.fq.gz pdom-rnaseq-q4-2.fq.gz \
          pdom-rnaseq-q5-1.fq.gz pdom-rnaseq-q5-2.fq.gz \
          pdom-rnaseq-q6-1.fq.gz pdom-rnaseq-q6-2.fq.gz \
          pdom-rnaseq-w1-1.fq.gz pdom-rnaseq-w1-2.fq.gz \
          pdom-rnaseq-w2-1.fq.gz pdom-rnaseq-w2-2.fq.gz \
          pdom-rnaseq-w3-1.fq.gz pdom-rnaseq-w3-2.fq.gz \
          pdom-rnaseq-w4-1.fq.gz pdom-rnaseq-w4-2.fq.gz \
          pdom-rnaseq-w5-1.fq.gz pdom-rnaseq-w5-2.fq.gz \
          pdom-rnaseq-w6-1.fq.gz pdom-rnaseq-w6-2.fq.gz

all:				pdom-trinity/Trinity.fasta $(LeftReads) $(RawFastqs) pdom-tsa-r1.2-unmatched-pep.fa

# Downloads precomputed mRNAmarkup results; too much of a hassle to automate here
pdom-tsa-r1.2-unmatched-pep.fa:	venn.make
				./venn.make BLASTTHREADS=$(NumThreads) all

pdom-trinity/Trinity.fasta:	pdom-rnaseq-all-trim-1.fq
				Trinity.pl --seqType fq \
				           --JM 100G \
				           --bflyHeapSpaceMax 50G \
				           --output pdom-trinity \
				           --CPU $(NumThreads) \
				           --left  pdom-rnaseq-all-trim-paired-1.fq \
				           --right pdom-rnaseq-all-trim-paired-2.fq \
				           --full_cleanup \
				           --jaccard_clip \
				           --CuffFly

pdom-rnaseq-all-trim-1.fq:	$(LeftReads)
				cat pdom-rnaseq-*-trim-1.fq > pdom-rnaseq-all-trim-1.fq
				cat pdom-rnaseq-*-trim-2.fq > pdom-rnaseq-all-trim-2.fq

pdom-rnaseq-%-trim-1.fq:	pdom-rnaseq-%-1.fq.gz pdom-rnaseq-%-2.fq.gz
				./run-trim $* $(TrimJar) $(NumThreads)

pdom-rnaseq-%.fq.gz:		
				iget $(PdomData)/sequence/transcriptome/$@
