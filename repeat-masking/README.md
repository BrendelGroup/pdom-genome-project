# Repeat masking

The genome assembly was screened for known repetitive elements using [RepeatMasker][] version [open-4.0.5][] and [Repbase][] version 20140131.
After masking repeats identified by RepeatMasker, the assembly was screened for additional repeats using [Tallymer] version [1.5.2].
To discriminate _bona fide_ repetitive elements from genes occurring in high copy number in the genome, all repeats identified by Tallymer were subjected to a BLASTX search against a database of Hexapod proteins.
Any repeats with matches in the database and e-values < 1e-5 were discarded as probable high copy number genes, while the rest were used to mask the genome.
The final masked sequence has been deposited in the Pdom data store at `r1.2/genome-assembly/pdom-scaffolds-masked-r1.2.fa.gz`.

## Procedure (interactive)

First, download the unmasked genome sequence.

```bash
PdomData=/iplant/home/standage/Polistes_dominula
iget -V $PdomData/r1.2/genome-assembly/pdom-scaffolds-unmasked-r1.2.fa.gz
gunzip pdom-scaffolds-unmasked-r1.2.fa.gz
```

### Screening with RepeatMasker

Next, identify known repeats with RepeatMasker.
By default, RepeatMasker produces soft-masked (lower-case) sequences, so we need to post-process the output to hard mask (N) the sequence.

```bash
NumThreads=16
GCContent=30.77
RepeatMasker -species insects -parallel $NumThreads -gc $GCContent \
             -frag 4000000 -lcambig -xsmall -gff \
             pdom-scaffolds-unmasked-r1.2.fa \
             > rm.log 2>&1
python lc2n.py < pdom-scaffolds-unmasked-r1.2.fa.masked \
               > pdom-rm-masked.fa
```

### Screening with Tallymer

The, do additional *k*-mer based screening for repetitive elements using Tallymer (procedure published by [Dan Bolser][]).

```bash
gt suffixerator -v \
                -db $IDX \
                -indexname $IDX \
                -tis -suf -lcp -des -ssp -sds -dna \
                > suffixerator.log 2>&1

gt tallymer occratio -v \
                     -minmersize 10 \
                     -maxmersize 45 \
                     -output unique nonunique nonuniquemulti total relative \
                     -esa $IDX \
                     > pdom.occratio.10.45.dump

gt tallymer mkindex -v \
            -mersize 19 \
            -minocc 50 \
            -esa $IDX \
            -counts -pl \
            -indexname pdom.idx.19.50 \
            > mkindex.log 2>&1

gt tallymer search -v \
            -output qseqnum qpos counts \
            -tyr pdom.idx.19.50 \
            -q $IDX \
            > pdom.repeats.19.50.tmer \
            2> tallymer.search.log

tallymer2gff3.plx -k 19 -s $IDX \
                  pdom.repeats.19.50.tmer \
                  > pdom.repeats.19.50.gff3

gff2fasta.plx -s pdom-rm-masked.fa \
              -f pdom.repeats.19.50.gff3 \
              > pdom.repeats.19.50.fa
```

Do a BLASTx search of repeats found by Tallymer vs known hexapod proteins, and parse out those with hits using [MuSeqBox][].

```bash
curl 'http://www.uniprot.org/uniprot/?query=taxonomy%3a6960&force=yes&format=fasta' \
     > hexapoda.fa
makeblastdb -in hexapoda.fa -dbtype prot -parse_seqids
blastx -query pdom.repeats.19.50.fa -db hexapoda.fa \
       -num_alignments 10 -evalue 1e-5 -num_threads 64 \
       -out pdom.repeats.19.50.blastx \
       > pdom.repeats.19.50.log 2>&1

MuSeqBox -i pdom.repeats.19.50.blastx -L 100 \
    | cut -f 1 -d ' ' | sort | uniq \
    | perl -ne 'm/(PdomSCFr1.2-\d+)-\d+\/(\d+)-(\d+)/ and print "$1\t$2\t$3\n"' \
    > pdom.repeats.19.50.hexapodhits.txt
```

Finally, hard mask the Tallymer repeats, excluding any that match Hexapod proteins as probably high-copy-number genes.

```bash
mask.pl pdom.repeats.19.50.gff3 \
        pdom.repeats.19.50.hexapodhits.txt \
        pdom-rm-masked.fa \
        > pdom-scaffolds-masked-r1.2.fa
gzip pdom-scaffolds-masked-r1.2.fa
```

## References

- **Smit AFA, Hubley R, Green P**. RepeatMasker Open-3.0. 1996-2010 <http://www.repeatmasker.org>. 
- **Jurka J, Kapitonov VV, Pavlicek A, Klonowski P, Kohany O, Walichiewicz J** (2005) Repbase Update, a database of eukaryotic repetitive elements. _Cytogentic and Genome Research_, **110**(1-4):462-467, [doi:10.1159/000084979](http://dx.doi.org/10.1159/000084979).
- **Kurtz S, Narechania A, Stein JC, Ware D** (2009) A new method to compute _K_-mer frequencies and its application to annotate large repetitive plant genomes. _BMC Genomics_, **9**:517, [doi:10.1186/1471-2164-9-517](http://dx.doi.ogr/10.1186/1471-2164-9-517).

[RepeatMasker]: http://www.repeatmasker.org/
[open-4.0.5]: http://www.repeatmasker.org/RepeatMasker-open-4-0-5.tar.gz
[Repbase]: http://www.girinst.org/server/RepBase/index.php
[Tallymer]: http://www.zbh.uni-hamburg.de/?id=211
[1.5.2]: http://genometools.org/pub/genometools-1.5.2.tar.gz
[Dan Bolser]: https://github.com/dbolser/PGSC/tree/master/kmer-filter
[MuSeqBox]: http://brendelgroup.org/bioinformatics2go/MuSeqBox.php
