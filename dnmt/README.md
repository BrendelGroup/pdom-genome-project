# Analysis of DNA methyltransferases

## Selecting probes

### DNMT1

- Search OrthoDB for [`dnmt1`](http://orthodb.org/orthodb7/results?tree=Arth&searchtext=dnmt1&level=Arthropoda&swaptree=)
- Retrieve GB number for longest *dnmt1* gene in *A. mellifera*: GB48403; retrieve sequence from [amel_OGSv3.2_pep.fa.gz](http://hymenopteragenome.org/beebase/?q=download_sequences)
- Do BLASTP search at NCBI: GB48403 vs Hymenoptera entries in NR
- Retrieve proteins from as many species as possible where length and sequence are similar: see [dnmt1-hym.faa](dnmt1-hym.faa)

### DNMT2

- Search OrthoDB for [`dnmt2`](http://orthodb.org/orthodb7/results?tree=Arth&searchtext=dnmt2&level=Arthropoda&swaptree=)
- Retrieve GB number for longest *dnmt2* gene in *A. mellifera*: GB54141; retrieve sequence from [amel_OGSv3.2_pep.fa.gz](http://hymenopteragenome.org/beebase/?q=download_sequences)
- Do BLASTP search at NCBI: GB54141 vs Hymenoptera entries in NR
- Retrieve proteins from as many species as possible where length and sequence are similar: see [dnmt2-hym.faa](dnmt2-hym.faa)

### DNMT3

- Search OrthoDB for [`dnmt3`](http://orthodb.org/orthodb7/results?tree=Arth&searchtext=dnmt3&level=Arthropoda&swaptree=)
- Retrieve GB number for longest *dnmt3* gene in *A. mellifera*: GB55485; retrieve sequence from [amel_OGSv3.2_pep.fa.gz](http://hymenopteragenome.org/beebase/?q=download_sequences)
  - Note that GB55485 is substantially longer than most other Arthropod DNMT3 entries in OrthoDB
- Do BLASTP search at NCBI: GB54141 vs Hymenoptera entries in NR
  - One of the first hits is [a DNMT3 protein which lists Wang et al paper as reference](http://www.ncbi.nlm.nih.gov/protein/NP_001177350.1); its length is much closer to the consensus length of the OrthoDB entries
  - reinitiate BLASTP search with this sequence as a query
- Retrieve proteins from as many species as possible where length and sequence are similar: see [dnmt3-hym.faa](dnmt3-hym.faa)

## Preliminary assessment of conservation

Used CLUSTAL to align each set of DNMTs and get a preliminary assessment their conservation.

```bash
clustalo --in dnmt1-hym.faa --outfmt clu --wrap 100 > dnmt1-hym.clu
clustalo --in dnmt2-hym.faa --outfmt clu --wrap 100 > dnmt2-hym.clu
clustalo --in dnmt3-hym.faa --outfmt clu --wrap 100 > dnmt3-hym.clu
```

See [dnmt1-hym.clu](dnmt1-hym.clu), [dnmt2-hym.clu](dnmt2-hym.clu), and [dnmt3-hym.clu](dnmt3-hym.clu).

## Next step: re-annotation

Use representative probe(s) from each set of DNMTs as input for a new VIGA run.

- splice align the probe to the genome
- cut out the genomic region containing the gene
- annotate the gene using CpGAT
- fill in a CGD database: conserved intron plots, alignments, gene trees
