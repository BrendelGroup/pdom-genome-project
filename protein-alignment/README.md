# Reference protein alignment

Reference protein sequences from Apis mellifera ([OGS 3.2] & [NCBI GNOMON][]) and Drosophila melanogaster ([FlyBase r5.55]) were splice-aligned to the genome using [GenomeThreader][] version 1.6.0.
The alignments have been deposited in the Pdom Data Store at `r1.2/protein-alignment`.

## Procedure (interactive)

First, download the assembled and masked genome sequence.
```bash
PdomData=/iplant/home/standage/Polistes_dominula
iget ${PdomData}/r1.2/genome-assembly/pdom-scaffolds-masked-r1.2.fa.gz
gunzip pdom-scaffolds-masked-r1.2.fa.gz
```

Next, download the reference proteins.

```bash
BeeData=http://hymenopteragenome.org/beebase/sites/hymenopteragenome.org.beebase
curl ${BeeData}/files/data/consortium_data/amel_OGSv3.2_pep.fa.gz \
    | zcat \
    | sed 's/gnl|Amel_4.5|//g' \
    > amel-ogs-prot.fa
curl ftp://ftp.ncbi.nih.gov/genomes/Apis_mellifera/protein/protein.fa.gz \
    | zcat \
    | perl -ne 's/>gi\|(\d+)\|ref\|([^|]+)\|\S*/>$1 $2/; print' \
    > amel-ncbi-prot.fa
FlyData=ftp://ftp.flybase.net/releases/FB2014_01/dmel_r5.55
curl ${FlyData}/fasta/dmel-all-translation-r5.55.fasta.gz \
    | zcat \
    > dmel-flybase-prot.fa
```

Then perform the alignments.

```bash
# The $GTHBIN variable should point to the bin directory in the GenomeThreader
# distribution.
export BSSMDIR=$GTHBIN/bssm
export GTHDATADIR=$GTHBIN/gthdata

for prot in amel-ogs amel-ncbi dmel-flybase
do
  $GTHBIN/gth -genomic pdom-scaffolds-masked-r1.2.fa \
              -protein ${prot}-prot.fa \
              -species arabidopsis \
              -gcmaxgapwidth 20000 \
              -gcmincoverage 25 \
              -prhdist 6 \
              -prminmatchlen 18 \
              -prseedlength 6 \
              -o ${prot}-prot-masked.gth \
              -force \
              > ${prot}-prot-masked.log 2>&1
done
```

## Procedure (automated)

The same procedure can also be run in batch mode using the following commands (in the `protein-alignment` directory).

```bash
make GTHbin=/usr/local/bin/GTH
```

## References

- **Gremme G, Brendel V, Sparks ME, Kurtz S** (2005) Engineering a software tool for gene structure prediction in higher organisms. *Information and Software Technology*, **47**(15):965-978, [doi:10.1016/j.infsof.2005.09.005](http://dx.doi.org/10.1016/j.infsof.2005.09.005).

[OGS 3.2]: http://hymenopteragenome.org/beebase/?q=download_sequences
[NCBI GNOMON]: ftp://ftp.ncbi.nih.gov/genomes/Apis_mellifera/protein
[FlyBase r5.55]: ftp://ftp.flybase.net/releases/FB2014_01/dmel_r5.55
[GenomeThreader]: http://genomethreader.org
