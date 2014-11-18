# Interval loci

The LocusPocus program ([AEGeAn Toolkit][] version [53d092091c][]) was used to compute interval loci (iLoci) from the genome annotation.
The xtractore program was used to extract the iLocus sequences from the genome sequence.
The iLocus sequences and annotations have been deposited in the Pdom Data Store at `/r1.2/interval-loci`.

## Procedure (interactive)

We first need the genome sequence and corresponding annotations.

```bash
PdomData=/iplant/home/standage/Polistes_dominula/r1.2
iget ${PdomData}/genome-annotation/pdom-annot-r1.2.gff3
iget ${PdomData}/genome-assembly/pdom-scaffolds-unmasked-r1.2.fa.gz
gunzip pdom-scaffolds-unmasked-r1.2.fa.gz
```

Then we compute the iLocus coordinates, excluding gene-less iLoci if they are at the end of a scaffold.

```bash
locuspocus --intloci --skipends \
           --delta=500 --verbose \
           --outfile=pdom-loci-r1.2.gff3 \
           pdom-annot-r1.2.gff3
```

The xtractore program gives us the iLocus sequences.

```bash
xtractore --type=locus --outfile=pdom-loci-r1.2.fa \
          pdom-loci-r1.2.gff3 \
          pdom-scaffolds-unmasked-r1.2.fa
```

For some analyses, we want gene coordinates expressed in reference to the iLocus to which they belong, as opposed to the entire scaffold.
We transformed the annotations to iLocus-based coordinates and created two files, one including tRNA genes and one without.

```bash
perl glb2lcl.pl < pdom-loci-r1.2.gff3 \
    | gt gff3 -retainids -sort -tidy -o pdom-annot-r1.2-iloci.gff3
    
perl glb2lcl.pl < pdom-loci-r1.2.gff3 \
    | grep -v -e trnascan -e PdomTRNA \
    | gt gff3 -retainids -sort -tidy -o pdom-annot-r1.2-iloci-sanstrna.gff3
```

## Procedure (automated)

The same procedure can also be run in batch mode using the following commands (in the `interval-loci` directory).

```bash
make
make clean
```

[AEGeAn Toolkit]: http://standage.github.io/AEGeAn
[53d092091c]: https://github.com/standage/AEGeAn/tree/53d092091c928136d5ab2d031dcd32f293ce3a4f
[d72e59a]: https://github.com/standage/AEGeAn/tree/d72e59ad0012f69fa2fd036f4d715af3ca72d1ab
