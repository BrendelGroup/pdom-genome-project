# Novel genes

Two previous analyses produced potential sources of novel genes warranting further investigation.

- The first potential source of novel genes are gene models predicted by Maker that have no matches against animal proteins in the NCBI nr database.
  These are described in the genome annotation documentation and are referred to hereafter as *unmatched gene models*.
- The second potential source of novel genes are assembled transcripts that have no matches in protein databases but are conserved among the 3 *Polistes* species.
  These are described in the transcript assembly documentation and are referred to hereafter as *unmatched TSAs*.

Here, we look for interval loci (iLoci) that have support from both of these sources of evidence.
Relevant data has been uploaded to the Pdom Data Store at `r1.2/novel-genes`.

## Procedure (interactive)

First we need to download the relevant data.

```bash
PdomData=/iplant/home/standage/Polistes_dominula/r1.2
iget ${PdomData}/transcript-assembly/pdom-tsa-r1.2-unmatched-pep.txt
iget ${PdomData}/transcript-alignment/pdom-tsa-masked-filtered.gff3
iget ${PdomData}/interval-loci/pdom-loci-r1.2.gff3
iget ${PdomData}/interval-loci/pdom-loci-r1.2-mrnamap.txt
iget ${PdomData}/genome-annotation/pdom-r1.2-no-anm-hits.ids
iget ${PdomData}/transcript-assembly/pdom-tsa-r1.2.fa.gz
gunzip pdom-tsa-r1.2.fa.gz
```

Next, identify iLoci to which unmatched TSAs align.

```bash
# Pull out alignment coordinates for unmatched TSAs
python align-keep.py pdom-tsa-r1.2-unmatched-pep.txt \
                     pdom-tsa-masked-filtered.gff3 \
    > pdom-tsa-aligned-unmatched.gff3

# Find the corresponding iLoci
bedtools intersect -a <(grep $'\tlocus\t' pdom-loci-r1.2.gff3) \
                   -b pdom-tsa-aligned-unmatched.gff3 -wa -u \
    | perl -ne 'm/(PdomILCr1.2-\d+)/ and print "$1\n"' \
    | sort | uniq > pdom-tsa-aligned-unmatched-loci.txt
```

Then look at unmatched gene models and determine the iLoci to which they correspond.

```bash
./selex --out=1 pdom-r1.2-no-anm-hits.ids pdom-loci-r1.2-mrnamap.txt \
    | sort | uniq > pdom-r1.2-no-anm-hits-iloci.txt
```

Finally, determine which iLoci appear in both lists (iLoci to which unmatched TSAs align and iLoci containing unmatched gene models).
```bash
comm -12 pdom-r1.2-no-anm-hits-iloci.txt pdom-tsa-aligned-unmatched-loci.txt \
    > pdom-loci-r1.2-pot-trgs.txt
```

## Procedure (automated)

The same procedure can also be run in batch mode using the following commands (in the `novel-genes` directory).

```bash
make
make clean
```
