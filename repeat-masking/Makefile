#!/usr/bin/make -f
SHELL=bash
NumThreads=1
GCContent=30.77
PdomData=/iplant/home/standage/Polistes_dominula

all:					pdom-scaffolds-masked-r1.2.fa.gz

pdom-scaffolds-masked-r1.2.fa.gz:	pdom-scaffolds-masked-r1.2.fa
					gzip pdom-scaffolds-masked-r1.2.fa

pdom-scaffolds-masked-r1.2.fa:		pdom.repeats.19.50.gff3 pdom.repeats.19.50.hexapodhits.txt pdom-rm-masked.fa
					./mask.pl $^ > $@

pdom.repeats.19.50.hexapodhits.txt:	pdom.repeats.19.50.blastx
					MuSeqBox -i pdom.repeats.19.50.blastx -L 100 \
					    | cut -f 1 -d ' ' | sort | uniq \
					    | perl -ne 'm/(PdomSCFr1.2-\d+)-\d+\/(\d+)-(\d+)/ and print "$$1\t$$2\t$$3\n"' \
					    > $@

pdom.repeats.19.50.blastx:		pdom.repeats.19.50.fa hexapoda.fa hexapoda.fa.pin
					blastx -query pdom.repeats.19.50.fa -db hexapoda.fa \
					       -num_alignments 10 -evalue 1e-5 -num_threads 64 \
					       -out $@ \
					    > pdom.repeats.19.50.log 2>&1

hexapoda.fa.pin:			hexapoda.fa
					makeblastdb -in hexapoda.fa -dbtype prot -parse_seqids

hexapoda.fa:				
					curl 'http://www.uniprot.org/uniprot/?query=taxonomy%3a6960&force=yes&format=fasta' \
					    > $@

pdom.repeats.19.50.fa:			pdom.repeats.19.50.gff3 pdom-rm-masked.fa
					gff2fasta.plx -s pdom-rm-masked.fa \
					              -f pdom.repeats.19.50.gff3 \
					    > $@

pdom.repeats.19.50.gff3:		pdom.repeats.19.50.tmer
					./tallymer2gff3.plx -k 19 -s $IDX \
					                    pdom.repeats.19.50.tmer \
					    > $@

pdom.repeats.19.50.tmer:		pdom.idx.19.50
					gt tallymer search -v \
					            -output qseqnum qpos counts \
					            -tyr pdom.idx.19.50 \
					            -q $IDX \
					    > $@ 2> tallymer.search.log

pdom.idx.19.50:				pdom.occratio.10.45.dump
					gt tallymer mkindex -v \
					            -mersize 19 \
					            -minocc 50 \
					            -esa $IDX \
					            -counts -pl \
					            -indexname $@ \
					    > mkindex.log 2>&1

pdom.occratio.10.45.dump:		suffixerator.log
					gt tallymer occratio -v \
					                     -minmersize 10 \
					                     -maxmersize 45 \
					                     -output unique nonunique nonuniquemulti total relative \
					                     -esa $IDX \
					    > $@

suffixerator.log:			pdom-rm-masked.fa
					gt suffixerator -v \
					                -db $IDX \
					                -indexname $IDX \
					                -tis -suf -lcp -des -ssp -sds -dna \
					    > suffixerator.log 2>&1

pdom-rm-masked.fa:			pdom-scaffolds-unmasked-r1.2.fa.masked
					./lc2n.py < $^ > $@

pdom-scaffolds-unmasked-r1.2.fa.masked:	pdom-scaffolds-unmasked-r1.2.fa
					RepeatMasker -species insects -parallel ${NumThreads} -gc ${GCContent} \
					             -frag 4000000 -lcambig -xsmall -gff \
					             pdom-scaffolds-unmasked-r1.2.fa \
					    > rm.log 2>&1

pdom-scaffolds-unmasked-r1.2.fa:		
					iget ${PdomData}/r1.2/genome-assembly/$@.gz
					gunzip $@.gz

clean:					
					echo FIXME
