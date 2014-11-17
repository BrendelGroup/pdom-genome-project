# Genome annotation

The [Maker annotation pipeline][] version 2.31.6 was used to annotate protein-coding genes and tRNA genes in the *P. dominula* genome.
The final, cleaned up annotation file and corresponding sequence files have been deposited in the Pdom Data Store at `r1.2/genome-annotation/`.

## Gene predictor training

Approximately 3000 genes were selected from annotations whose reliability was confirmed by strict measures of conservation with several other Hymenopteran species.
These genes were then used to train 3 ab initio gene predictors (SNAP, Augustus, and GeneMark) for use with the Maker pipeline.
For Augustus and SNAP, the procedures documented in their respective source code distributions were followed.
For GeneMark, a model built by the self-training version of GeneMark was used with Maker.
All models have been deposited in the Pdom Data Store at `r1.2/genome-annotation/models`.

## Procedure (interactive)

### Data files

First, there are *a lot* of data files to download.
Make sure to set the `AUGUSTUS_CONFIG_DIR` variable appropriately.

```bash
PdomData=/iplant/home/standage/Polistes_dominula
AUGUSTUS_CONFIG_DIR=/usr/local/src/augustus/config
export AUGUSTUS_CONFIG_DIR

# Genome sequence
iget ${PdomData}/genome-assembly/pdom-scaffolds-masked-r1.2.fa.gz
gunzip pdom-scaffolds-masked-r1.2.fa.gz

# Transcript alignments
iget ${PdomData}/r1.2/transcript-alignment/pdom-r1.2-trans-align.tar.gz
tar -xzf pdom-r1.2-trans-align.tar.gz
mv pdom-r1.2-trans-align/*.gff3 .

# Protein alignments
iget ${PdomData}/r1.2/protein-alignment/pdom-r1.2-refrprot-align.tar.gz
tar -xzf pdom-r1.2-refrprot-align.tar.gz
mv pdom-r1.2-refrprot-align/*.gff3 .

# Highly conserved and manually curated annotations
iget ${PdomData}/r1.2/genome-annotation/pdom-viga-improved.gff3.gz
iget ${PdomData}/r1.2/genome-annotation/pdom-yrgate.gff3.gz
gunzip pdom-viga-improved.gff3.gz
gunzip pdom-yrgate.gff3.gz

# Download species-specific parameter settings for ab initio gene predictors
iget ${PdomData}/r1.2/genome-annotation/models/pdom.snap.hmm
iget ${PdomData}/r1.2/genome-annotation/models/pdom.genemark.mod
iget ${PdomData}/r1.2/genome-annotation/models/pdom.augustus.tar.gz
tar -xzf pdom.augustus.tar.gz
mv pdom ${AUGUSTUS_CONFIG_DIR}/species
```

### Configuration files

Next, we need to create the Maker control files.
If all of the supplementary programs are the in the path, the ``maker_exe.ctl`` file will be populated automatically.
If not, you will need to manually fill it in with the location of all the programs.

```bash
maker -CTL
# Replace the empty `maker_opts.ctl` file with our file,
# which has all the values filled in.
cp pdom_maker_opts.ctl maker_opts.ctl
```

### Annotation with Maker

After all of the data files and control files are in place, running Maker is trivial.

```bash
NumThreads=16
maker -genome pdom-scaffolds-masked-r1.2.fa \
      -fix_nucleotides \
      -nodatastore \
      -RM_off \
      -cpus $NumThreads \
      -base pdom \
      > pdom.maker.log 2>&1
```

### Post-processing

The feature IDs created internally by Maker are pretty unwieldy, so we use a few scripts to clean them up.

```bash
# Merge and clean up data
gff3_merge -o pdom-annot-p1.2-raw.gff3 pdom.maker.output/pdom_master_datastore_index.log
bash clean.sh pdom-annot-p1.2-raw.gff3 > pdom-annot-p1.2-cleaned.gff3
# Fixes strange output artifact for some VIGA-based predictions
perl -n fix.pl < pdom-annot-p1.2-cleaned.gff3 > pdom-annot-p1.2-fixed.gff3

# Create a polished GFF3 file with official IDs for the annotations and proper
# ##sequence-region pragmas.
gt gff3 -retainids -sort -tidy pdom-annot-p1.2-fixed.gff3 2> gt.log \
    | python annot-ids.py --idfmt='Pdom%sr1.2-%05lu' -n --rnamap=rnaids.txt --dbxref=MAKER - \
    | python fix-regions.py pdom-scaffolds-masked-r1.2.fa > pdom-annot-r1.2.gff3

# Create transcript and protein sequence files with proper feature IDs
cat pdom.maker.output/pdom_datastore/PdomSCFr1.2-*/PdomSCFr1%2E2-*.maker.transcripts.fasta \
    | python seq-ids.py rnaids.txt > pdom-annot-r1.2-transcripts.fasta
cat pdom.maker.output/pdom_datastore/PdomSCFr1.2-*/PdomSCFr1%2E2-*.maker.proteins.fasta \
    | python seq-ids.py rnaids.txt > pdom-annot-r1.2-proteins.fasta
```

### Functional annotation by BLASTp search

Translation products of all gene models predicted by Maker were searched against all animal proteins in the NCBI nr database.
This BLAST search provides a preliminary functional annotation for the gene models, but also identified gene models without any matches to known proteins.

```bash
# NCBI nr database installed using update_blastdb.pl script
# distributed with BLAST.
blastdb_aliastool -gilist anm.p.gil -dbtype prot -db nr -out anm

# BLAST search
blastp -db anm \
       -query pdom-annot-r1.2-proteins.fa \
       -evalue 1e-4 \
       -num_threads $NumThreads \
       -outfmt 5 \
       -out pdom-r1.2-vs-anm-blastp.xml
```

Using the NCBI Taxonomy data, we can examine the best BLAST hit for each gene model and tally these up according to the family of the best hit.

```bash
# Download NCBI taxonomy data
mkdir taxonomy
cd taxonomy
curl -O ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xzf taxdump.tar.gz
cd -

perl blastxml-to-besthit.pl taxonomy/names.dmp \
    < pdom-r1.2-vs-anm-blastp.xml \
    > pdom-annot-r1.2-besthit-gis.txt

# Breakdown of hits to Hymenoptera proteins, by family
./taxtrav --rank family \
        --filter 'order=Hymenoptera' \
        --nodes taxonomy/nodes.dmp \
        --names taxonomy/names.dmp \
        <(cut -f 2 pdom-annot-r1.2-besthit-gis.txt)
    > pdom-hym-hits-tax.csv
cut -f 3 -d , pdom-hym-hits-tax.csv | sort | uniq -c | sort -rn

# Number of hits to non-Hymenoptera proteins
./taxtrav --rank family \
        --filter 'order=Hymenoptera' \
        --nodes taxonomy/nodes.dmp \
        --names taxonomy/names.dmp \
        --complement \
        <(cut -f 2 pdom-annot-r1.2-besthit-gis.txt)
    > pdom-nothym-hits-tax.csv
wc -l pdom-nothym-hits-tax.csv
```

Finally, let's pull out the gene models with no matches in the database.
Unmatched gene models tend to be short, but we can  also pull out gene models whose length (in amino acids) exceeds the median length for all gene models and set these aside as for further examination.

```bash
comm -13 <(cut -f 1 pdom-annot-r1.2-besthit-gis.txt | sort) \
         <(perl -ne 'm/>(\S+)/ and print "$1\n"' < pdom-annot-r1.2-proteins.fa) \
    > pdom-r1.2-no-anm-hits.ids
./select-seq pdom-r1.2-no-anm-hits.ids \
             pdom-annot-r1.2-proteins.fa \
    > pdom-r1.2-no-anm-hits.fa
bash fasta-median-plus.sh pdom-annot-r1.2-proteins.fa \
                          pdom-r1.2-no-anm-hits.fa \
    > pdom-annot-r1.2-long-unmatched.fa
```

## Procedure (automated)

The same procedure can also be run in batch mode using the following commands (in the `genome-annotation` directory).
Keep in mind that unless all of the supplemental programs are in your path, you will still need to run `maker -CTL` and edit the `maker_exe.ctl` file manually.

```bash
make NumThreads=16
```

## References

- **Cantarel BL, Korf I, Robb SMC, Parra G, Ross E, Moore B, Holt C, Alvarado AS, Yandell M** (2008) MAKER: An easy-to-use annotation pipeline designed for emerging model organism genomes. *Genome Research*,  18:88-196, [doi:10.1101/gr.6743907](http://dx.doi.org/10.1101/gr.6743907).

[Maker annotation pipeline]: http://www.yandell-lab.org/software/maker.html
