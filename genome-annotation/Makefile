#!/usr/bin/make -f
SHELL=bash
NumThreads=1
PdomData=/iplant/home/standage/Polistes_dominula
AUGUSTUS_CONFIG_DIR=/usr/local/src/augustus/config


all:					pdom-annot-r1.2-long-unmatched.fa pdom-hym-hits-tax.csv pdom-nothym-hits-tax.csv

#--------------------------------------------------------------------------------
# BLASTp search
#--------------------------------------------------------------------------------
pdom-annot-r1.2-long-unmatched.fa:	pdom-annot-r1.2-proteins.fa pdom-r1.2-no-anm-hits.fa
					bash fasta-median-plus.sh $^ > $@

pdom-r1.2-no-anm-hits.fa:		pdom-r1.2-no-anm-hits.ids pdom-annot-r1.2-proteins.fa
					select-seq $^ > $@

pdom-r1.2-no-anm-hits.ids:		pdom-annot-r1.2-besthit-gis.txt pdom-annot-r1.2-proteins.fa
					comm -13 <(cut -f 1 pdom-annot-r1.2-besthit-gis.txt | sort) \
					         <(perl -ne 'm/>(\S+)/ and print "$1\n"' < pdom-annot-r1.2-proteins.fa) \
					    > $@

pdom-nothym-hits-tax.csv:		pdom-annot-r1.2-besthit-gis.txt taxtrav
					./taxtrav --rank family \
					          --filter 'order=Hymenoptera' \
					          --nodes taxonomy/nodes.dmp \
					          --names taxonomy/names.dmp \
					          --complement \
					          <(cut -f 2 pdom-annot-r1.2-besthit-gis.txt)
					    > pdom-nothym-hits-tax.csv

pdom-hym-hits-tax.csv:			pdom-annot-r1.2-besthit-gis.txt taxtrav
					taxtrav --rank family \
					        --filter 'order=Hymenoptera' \
					        --nodes taxonomy/nodes.dmp \
					        --names taxonomy/names.dmp \
					        <(cut -f 2 pdom-annot-r1.2-besthit-gis.txt)
					    > pdom-hym-hits-tax.csv

taxtrav:				taxtrav.c
					gcc -Wall -O3 -o taxtrav taxtrav.c

pdom-annot-r1.2-besthit-gis.txt:	taxonomy/names.dmp pdom-r1.2-vs-anm-blastp.xml
					perl blastxml-to-besthit.pl taxonomy/names.dmp \
					    < pdom-r1.2-vs-anm-blastp.xml \
					    > $@

taxonomy/names.dmp:			
					mkdir taxonomy
					cd taxonomy
					curl -O ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
					tar xzf taxdump.tar.gz
					cd ..

pdom-r1.2-vs-anm-blastp.xml:		anm.pal pdom-annot-r1.2-proteins.fa
					blastp -db anm \
					       -query pdom-annot-r1.2-proteins.fa \
					       -evalue 1e-4 \
					       -num_threads $NumThreads \
					       -outfmt 5 \
					       -out $@

anm.pal:				nr anm.p.gil
					blastdb_aliastool -gilist anm.p.gil -dbtype prot -db nr -out anm

#--------------------------------------------------------------------------------
# Run Maker
#--------------------------------------------------------------------------------
pdom-annot-r1.2-proteins.fa:		pdom-annot-r1.2.gff3
					cat pdom.maker.output/pdom_datastore/PdomSCFr1.2-*/PdomSCFr1%2E2-*.maker.proteins.fasta \
					    | python seq-ids.py rnaids.txt > $@

pdom-annot-r1.2-transcripts.fa:		pdom-annot-r1.2.gff3
					cat pdom.maker.output/pdom_datastore/PdomSCFr1.2-*/PdomSCFr1%2E2-*.maker.transcripts.fasta \
					    | python seq-ids.py rnaids.txt > $@

pdom-annot-r1.2.gff3:			pdom-annot-p1.2-raw.gff3
					bash clean.sh pdom-annot-p1.2-raw.gff3 \
					    | perl -n fix.pl > pdom \
					    > pdom-annot-p1.2-fixed.gff3
					gt gff3 -retainids -sort -tidy pdom-annot-p1.2-fixed.gff3 2> gt.log \
					    | python annot-ids.py --idfmt='Pdom%sr1.2-%05lu' -n --rnamap=rnaids.txt --dbxref=MAKER - \
					    | python fix-regions.py pdom-scaffolds-unmasked-r1.2.fa > pdom-annot-r1.2.gff3

pdom-annot-p1.2-raw.gff3:		pdom-scaffolds-masked-r1.2.fa pdom-tsa-r1.2-masked.gff3 \
					amel-ncbi-prot.gff3 pdom-viga-improved.gff3 \
					pdom-yrgate.gff3 pdom.snap.hmm
					maker -CTL
					cp pdom_maker_opts.ctl maker_opts.ctl
					maker -genome pdom-scaffolds-masked-r1.2.fa \
					      -fix_nucleotides \
					      -nodatastore \
					      -RM_off \
					      -base pdom \
					      > pdom.maker.log 2>&1
					gff3_merge -o pdom-annot-p1.2-raw.gff3 \
					              pdom.maker.output/pdom_master_datastore_index.log

#--------------------------------------------------------------------------------
# Download data files
#--------------------------------------------------------------------------------

pdom.snap.hmm:				
					iget $(PdomData)/r1.2/genome-annotation/models/pdom.snap.hmm
					iget $(PdomData)/r1.2/genome-annotation/models/pdom.genemark.mod
					iget $(PdomData)/r1.2/genome-annotation/models/pdom.augustus.tar.gz
					tar -xzf pdom.augustus.tar.gz
					mv pdom $(AUGUSTUS_CONFIG_DIR)/species

pdom-yrgate.gff3:		
					$(PdomData)/r1.2/genome-annotation/$@.gz
					gunzip $@.gz

pdom-viga-improved.gff3:		
					$(PdomData)/r1.2/genome-annotation/$@.gz
					gunzip $@.gz

amel-ncbi-prot.gff3:			
					iget $(PdomData)/r1.2/protein-alignment/pdom-r1.2-refrprot-align.tar.gz
					tar -xzf pdom-r1.2-refrprot-align.tar.gz
					mv pdom-r1.2-refrprot-align/*.gff3 .

pdom-tsa-r1.2-masked.gff3:		
					iget $(PdomData)/r1.2/transcript-alignment/pdom-r1.2-trans-align.tar.gz
					tar -xzf pdom-r1.2-trans-align.tar.gz
					mv pdom-r1.2-trans-align/*.gff3 .

pdom-scaffolds-masked-r1.2.fa:		
					iget $(PdomData)/genome-assembly/$@.gz
					gunzip $@.gz
