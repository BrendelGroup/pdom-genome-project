# Transcript alignment

Assembled and groomed *P. dominula* transcripts were spliced aligned at high stringency to the repeat-masked genome sequence using [GeneSeqer][] version [2-26-2014][].
*Polistes canadensis* and *P. metricus* transcripts were also aligned with GeneSeqer using less stringent parameters.
The alignments have been deposited in the Pdom Data Store at `r1.2/transcript-alignment`.

## Procedure (interactive)

First, set the `GSQDir` variable to the path of the GeneSeqer source code.

```bash
GSQDir=/usr/local/src/GENESEQER
```

Next, align *P. dominula* TSAs with stringent parameters.

```bash
# Download the assembled and masked genome sequence
PdomData=/iplant/home/standage/Polistes_dominula
iget ${PdomData}/r1.2/genome-assembly/pdom-scaffolds-masked-r1.2.fa.gz
gunzip pdom-scaffolds-masked-r1.2.fa.gz

# P. dominula alignments
iget ${PdomData}/r1.2/transcript-assembly/pdom-tsa-r1.2.fa.gz
gunzip pdom-tsa-r1.2.fa.gz
MakeArray pdom-tsa-r1.2.fa
GeneSeqerL -s Arabidopsis \
           -L pdom-scaffolds-masked-r1.2.fa \
           -D pdom-tsa-r1.2.fa.gz \
           -O pdom-tsa-r1.2-masked.gsq \
           -p $GSQDIR/data/prmfileHQ \
           -x 30 -y 45 -z 60 -w 0.8 -m 1000000 \
           > pdom-tsa-r1.2-masked-gsq.log 2>&1
```

Now the *P. canadensis* and *P, metricus* TSAs with less stringent parameters.

```bash
# P. metricus alignments
iget ${PdomData}/r1.2/transcript-assembly/pmet-tsa-r1.2.fa.gz
gunzip pmet-tsa-r1.2.fa.gz
MakeArray pmet-tsa-r1.2.fa
GeneSeqerL -s Arabidopsis \
           -L pdom-scaffolds-masked-r1.2.fa \
           -d pmet-tsa-r1.2.fa \
           -O pmet-tsa-r1.2-masked.gsq \
           -p $GSQDIR/data/prmfile \
           -x 16 -y 24 -z 48 -w 0.8 -m 1000000 \
           > pmet-tsa-r1.2-masked-gsq.log 2>&1

# P. canadensis alignments
iget ${PdomData}/r1.2/transcript-assembly/pcan-tsa.fa.gz
gunzip pcan-tsa.fa.gz
MakeArray pcan-tsa.fa
GeneSeqerL -s Arabidopsis \
           -L pdom-scaffolds-masked-r1.2.fa \
           -d pcan-tsa.fa \
           -O pcan-tsa-masked.gsq \
           -p $GSQDIR/data/prmfile \
           -x 16 -y 24 -z 48 -w 0.8 -m 1000000 \
           > pcan-tsa-masked-gsq.log 2>&1
```

Lastly, convert all of the alignments into GFF3 format, filtering out alignments with GeneSeqer similarity or coverage scores < 0.8.

```bash
./gsq2makergff3.py < pdom-tsa-r1.2-masked.gsq > pdom-tsa-r1.2-masked.gff3
./gsq2makergff3.py < pmet-tsa-r1.2-masked.gsq > pmet-tsa-r1.2-masked.gff3
./gsq2makergff3.py < pcan-tsa-r1.2-masked.gsq > pcan-tsa-masked.gff3
```

## Procedure (automated)

The same procedure can also be run in batch mode using the following commands (in the `transcript-alignment` directory).

```bash
make GSQDir=/usr/local/src/GENESEQER
```

## References

- **Brendel V, Xing L, Zhu W** (2004) Gene structure prediction from consensus spliced alignment of multiple ESTs matching the same genomic locus. *Bioinformatics*, **20**:1157-1169, [doi:10.1093/bioinformatics/bth058](http://dx.doi.org/10.1093/bioinformatics/bth058).
- **Berens _et al_**, in preparation.
- **Ferreira P, Patalano S, Chauhan R, Ffrench-Constant R, Gabaldon T, Guigo R, Sumner S** (2013) Transcriptome analyses of primitively eusocial wasps reveal novel insights into the evolution of sociality and the origin of alternative phenotypes. *Genome Biology*, **14**(2):R20, [doi:10.1186/gb-2013-14-2-r20](http://dx.doi.org/10.1186/gb-2013-14-2-r20).

[GeneSeqer]: http://brendelgroup.org/bioinformatics2go/GeneSeqer.php
[2-26-2014]: http://www.brendelgroup.org/bioinformatics2go/Download/GeneSeqer-2-26-2014.tar.gz
