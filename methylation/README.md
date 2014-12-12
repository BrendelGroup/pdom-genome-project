# Analysis of DNA methyltransferases

## Selecting probes

Repeated the following steps for Dnmt1, Dnmt2, and Dnmt3, MDB (methyl-DNA binding domain), and TET2 (methylcytosine dioxygenase) proteins.

- Search OrthoDB for dnmt\* (example: [`dnmt1`](http://orthodb.org/orthodb7/results?tree=Arth&searchtext=dnmt1&level=Arthropoda&swaptree=))
- Retrieve GB number for longest *dnmt*\* gene in *A. mellifera*: GB48403; retrieve sequence from [amel_OGSv3.2_pep.fa.gz](http://hymenopteragenome.org/beebase/?q=download_sequences)
  - Dnmt1: GB48403
  - Dnmt2: GB54141
  - Dnmt3: GB55485
  - MBD: GB41522
  - TET2: gi|571563963|ref|XP_006561260.1 (could not find in OrthoDB; found with NCBI keyword search
- Do BLASTP search at NCBI: *A. mellifera* protein vs Hymenoptera entries in NR
- Retrieve proteins from as many species as possible where length and sequence are similar
  - see [dnmt1-hym.faa](dnmt1-hym.faa)
  - see [dnmt2-hym.faa](dnmt2-hym.faa)
  - see [dnmt3-hym.faa](dnmt3-hym.faa)
  - see [mbd-hym.faa](mdb-hym.faa)
- Do BLASTP search at PdomGDB to retrieve *P. dominula* homolog
- Go back to OrthoDB, select the "Vertebrate" tab, find the link for the Dnmt\* gene in *Mus musculus*; download from Ensembl, create a new file with mouse as outgroup
  - see [dnmt1-vert.faa](dnmt1-vert.faa)
  - see [dnmt2-vert.faa](dnmt2-vert.faa)
  - see [dnmt3-vert.faa](dnmt3-vert.faa)
  - see [mbd-vert.faa](mbd-vert.faa)

### Notes

- The *A. mellifera* Dnmt3 entry in OrthoDB is substantially longer than most other Arthropod Dnmt3 entries.
  One of the first BLASTP hits of this protein against NR is [a DNMT3 protein which lists Wang et al paper as reference](http://www.ncbi.nlm.nih.gov/protein/NP_001177350.1); its length is much closer to the consensus length of the OrthoDB entries.
  This sequence therefore replaced the original *A. mellifera* probe and was used to reinitiate the BLASTP search for Dnmt3.
- *P. dominula* does not have a Dnmt3 gene, so it is not included in that file.
- The *P. dominula* TET2 gene was incorrectly annotated by Maker, requiring manual refinement with yrGATE.
- The Hymenopteran TET2 gene does not appear to be conserved in mammals.

## Preliminary assessment of conservation

Used CLUSTAL to align each set of DNMTs and get a preliminary assessment their conservation.

```bash
clustalo --in dnmt1-hym.faa --outfmt clu --wrap 100 > dnmt1-hym.clu
clustalo --in dnmt2-hym.faa --outfmt clu --wrap 100 > dnmt2-hym.clu
clustalo --in dnmt3-hym.faa --outfmt clu --wrap 100 > dnmt3-hym.clu
clustalo --in mbd-hym.faa --outfmt clu --wrap 100 > mdb-hym.clu
clustalo --in tet2-hym.faa --outfmt clu --wrap 100 > tet2-hym.clu
clustalo --in dnmt1-vert.faa --outfmt clu --wrap 100 > dnmt1-vert.clu
clustalo --in dnmt2-vert.faa --outfmt clu --wrap 100 > dnmt2-vert.clu
clustalo --in dnmt3-vert.faa --outfmt clu --wrap 100 > dnmt3-vert.clu
clustalo --in mbd-vert.faa --outfmt clu --wrap 100 > mbd-vert.clu
```

Or just run `make`.

See the following alignment files.

- just Hymenoptera
  - [dnmt1-hym.clu](dnmt1-hym.clu)
  - [dnmt2-hym.clu](dnmt2-hym.clu)
  - [dnmt3-hym.clu](dnmt3-hym.clu)
  - [mbd-hym.clu](mbd-hym.clu)
  - [tet2-hym.clu](tet2-hym.clu)
- Hymenoptera plus vertebrate outgroup
  - [dnmt1-vert.clu](dnmt1-vert.clu)
  - [dnmt2-vert.clu](dnmt2-vert.clu)
  - [dnmt3-vert.clu](dnmt3-vert.clu)
  - [mbd-vert.clu](mbd-vert.clu)

## Next step: re-annotation

Use representative probe(s) from each set of DNMTs as input for a new VIGA run.

- splice align the probe to the genome
- cut out the genomic region containing the gene
- annotate the gene using CpGAT
- fill in a CGD database: conserved intron plots, alignments, gene trees
