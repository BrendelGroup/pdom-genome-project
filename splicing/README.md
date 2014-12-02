# Analysis of alternative and differential splicing

The Tuxedo suite and the AEGeAn Toolkit were used to catalog alternative splicing events and identify genes that are differentially spliced between queens and workers.
RNA-seq reads (trimmed for quality control\*) were mapped to the filtered iLoci\* using [Tophat][] version [2.0.12][]  and [Bowtie2][] version [2.2.3][].
[Cufflinks][] version [2.2.1][] was then used to assemble the aligned reads perform differential splicing tests, while the ``asinspect`` program ([AEGeAn][] version [4fc8909d9b][]) was used to construct a catalog of alternative splicing events.
The output has been deposited in the Pdom Data Store at ``r1.2/splicing``.

\*Described in the section on differential expression.

## Procedure

### RNA-Seq alignments

First we align the RNA-Seq reads to the genome (iLoci) sample-by-sample.
We start by downloading and indexing the iLocus sequences.

```bash
mkdir -p genome
PdomData=/iplant/home/standage/Polistes_dominula/
IlocusData=${PdomData}/r1.2/interval-loci
iget ${IlocusData}/pdom-loci-r1.2-filtered.fa genome/
iget ${IlocusData}/pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 genome/

cd genome
bowtie2-build pdom-loci-r1.2-filtered.fa pdom-loci-r1.2
cd -
```

Then we perform the alignments.
Make sure the groomed RNA-Seq data\* is in the `rnaseq-clean` directory.

```bash
NumThreads=16
mkdir -p alignments
for caste in q w
do
  for rep in {1..6}
  do
    sample=${caste}${rep}
    tophat --output-dir alignments/${sample} \
           --GTF genome/pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 \
           --mate-inner-dist 70 \
           --min-anchor 6 \
           --num-threads $NumThreads \
           genome/pdom-loci-r1.2 \
           rnaseq-clean/pdom-rnaseq-${sample}-trim-paired-1.fq \
           rnaseq-clean/pdom-rnaseq-${sample}-trim-paired-2.fq
  done
done
```

### Assemblies

Next we need to assemble mapped reads to construct transcripts.
We first assemble each sample individually.

```bash
mkdir -p assemblies
rm -f assemblies/all.txt
for caste in q w
do
  for rep in {1..6}
  do
    sample=${caste}${rep}
    echo assemblies/${sample}/transcripts.gtf >> assemblies/all.txt
    cufflinks --output-dir assemblies/${sample} \
              --num-threads $NumThreads \
              --multi-read-correct \
              --GTF-guide genome/pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 \
              alignments/${sample}/accepted_hits.bam
  done
done
```

Then we combine all assemblies into a single aggregate assembly.

```bash
cuffmerge -o assemblies/merged \
          --num-threads $NumThreads \
          --ref-sequence genome/pdom-loci-r1.2-filtered.fa \
          --ref-gtf genome/pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 \
          assemblies/all.txt
```

For tools that expect transcript structures in GFF3 format, we convert the final assembly from GTF to GFF3 (removing strandless single exon features).

```bash
gt gtf_to_gff3 -tidy assemblies/merged/merged.gtf \
    | perl -ne '@f = split(/\t/); next if($f[6] eq "."); print' \
    | gt gff3 -sort -tidy -o pdom-annot-r1.2-tuxedo.gff3 -force
```

### Differential splicing

With a combined assembly, we can analyze the queen and worker data and look for genes showing caste-dependent alternative splicing patterns.

```bash
qaligns=$(ls alignments/q?/accepted_hits.bam | tr '\n' ',')
waligns=$(ls alignments/w?/accepted_hits.bam | tr '\n' ',')
qaligns=${qaligns%?}
waligns=${waligns%?}
cuffdiff --output-dir diff \
         --labels Queen,Worker \
         --num-threads $NP \
         --multi-read-correct \
         assemblies/merged/merged.gtf \
         $qaligns $waligns
```

### Catalog alternative splicing events

Cuffdiff examines differential splicing patterns, but we also want to compose a catalog of alternative splicing events, ignoring whether the castes exhibit any difference in isoform preference.

```bash
asinspect --refr genome/pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 \
          --gtf \
          --out pdom-as.tsv \
          assemblies/merged/merged.gtf \
```

The ``asinspect`` program reports the following classes of alternative splicing.

  - *Cassette exons*: an exon from the reference annotation is designated a cassette exon (CE) if its flanking exons are spliced together in any of the isoforms present in the Cufflinks assembly.
    The program also reports novel exons from the Cufflinks assembly that are present between 2 exons included in the reference assembly.
  - *Retained introns*: an intron from the reference annotation is designated a retained intron (RI) if a single exon from the Cufflinks assembly matches the outer boundaries of the intron's flanking exons.
    The program also reports novel introns from the Cufflinks assembly that split a single exon from the Maker annotation into two new exons sharing the same outer boundaries.

We are particularly interested in further investigating conserved CE events.
First we identify cassette exons that are entirely coding exons (i.e. contain no start or stop codon).

```bash
python coding-cassette-exons.py \
       pdom-as.tsv pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 \
       > pdom-as-ce-coding.tsv
```

Next we create a GFF3 and Fasta file containing both the exon-skipped isoform and the exon-contained isoform.

```bash
python ce-isoforms.py \
       pdom-as-ce-coding.tsv pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 \
    | gt gff3 -retainids -sort -tidy \
    > pdom-annot-r1.2-withceisoforms.gff3

gt extractfeat -matchdescstart -join -retainids -translate \
               -seqfile pdom-loci-r1.2.fa -type CDS \
               pdom-annot-r1.2-withceisoforms.gff3 \
    | perl get_ce_iso.pl \
    > pdom-annot-r1.2-ce-isoforms.fa
```

Finally, we do a BLASTp search of these isoforms against *Apis mellifera* proteins.
To identify high-confidence conserved CE events, we used MuSeqBox to filter the BLAST output and select the genes satisfying the following criteria.

  - both isoforms have a match with expect value < 1e-20
  - $>=$ 90% reciprocal coverage between each query and match
  - the CE-containing isoform has a different best match than the CE-skipped isoform

```bash
curl -O ftp://ftp.ncbi.nih.gov/genomes/Apis_mellifera/protein/protein.fa.gz
gunzip protein.fa.gz
mv protein.fa amel-ncbi-prot.fa
makeblastdb -in amel-ncbi-prot.fa -dbytype prot -parse_seqids
blastp -db amel-ncbi-prot.fa -query pdom-annot-r1.2-ce-isoforms.fa \
       -evalue 1e-6 -num_threads $NumThreads -out pdom-vs-amel.blastp

MuSeqBox -i pdom-vs-amel.blastp -n 1 -s 1 -l 25 -d 4 -L 128 -c ce-crtfile \
    > pdom-vs-amel-blastp.msb

python ce-filter.py < pdom-vs-amel-blastp.msb > pdom-vs-amel-conserved.txt
cut -f 1 -d '=' pdom-vs-amel-conserved.txt | sort | uniq \
    > pdom-annot-r1.2-ce-list.txt
```

## Procedure (automated)

The same procedure can also be run in batch mode using the following commands (in the `splicing` directory).

```bash
bash procedure.sh . 16
```

## References

- **Kim D, Pertea G, Trapnell C, Pimentel H, Kelley R, Salzberg SL** (2013) TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions. *Genome Biology*, **14**:R36, [doi:10.1186/gb-2013-14-4-r36](http://dx.doi.org/10.1186/gb-2013-14-4-r36).
- **Langmead B, Salzberg SL** (2012) Fast gapped-read alignment with Bowtie 2. *Nature Methods*, **9**:357-359, [doi:10.1038/nmeth.1923](http://dx.doi.org/10.1038/nmeth.1923).
- **Trapnell C, Hendrickson D, Sauvageau S, Goff L, Rinn JL, Pachter L** (2013) Differential analysis of gene regulation at transcript resolution with RNA-seq. *Nature Biotechnology*, **31**:46-53, [doi:10.1038/nbt.2450](http://dx.doi.org/10.1038/nbt.2450).

[Tophat]: http://ccb.jhu.edu/software/tophat/index.shtml
[2.0.12]: http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.12.tar.gz
[Bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
[2.2.3]: http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.3/bowtie2-2.2.3-source.zip
[Cufflinks]: http://cufflinks.cbcb.umd.edu/
[2.2.1]: http://cufflinks.cbcb.umd.edu/downloads/cufflinks-2.2.1.tar.gz
[AEGeAn]: http://standage.github.io/AEGeAn
[4fc8909d9b]: https://github.com/standage/AEGeAn/tree/4fc8909d9b0fd73f7bfa5275c7cbe5e8000c558e
