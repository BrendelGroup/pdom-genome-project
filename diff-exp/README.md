# Differential expression analysis

iLocus expression levels were determined for each replicate using [RSEM][] version [1.2.7][], which in turn utilizes [Bowtie][] version [1.0.0][] to perform read alignments.
[EBSeq][] version 1.1.5 (bundled with RSEM) was used to identify genes that are differentially expressed between queens and workers.
The [tk-rnaseq toolkit][] version [5ce0f8c075][] was used at various stages to filter, process, and visualize the data.
Various data tables, sequence files, and graphics described here have been deposited in the Pdom Data Store at ``r1.2/diff-exp``.

## Procedure (interactive)

The differential expression analysis took place in 3 stages.
In the first stage we conducted the full analysis to collect information about each iLocus and each sample.

```bash
# Ensure the RSEM directory is in the $PATH variable.
# Download and trim RNA-Seq data as described in the
# transcriptome assembly.

# Download the iLoci
mkdir iloci
PdomData=/iplant/home/standage/Polistes_dominula
iget ${PdomData}/r.1/interval-loci/pdom-loci-r1.2.fa iloci/

# Prepare/index the iLocus sequences
mkdir -p logs/unfiltered
rsem-prepare-reference iloci/pdom-loci-r1.2.fa \
                       iloci/pdom-loci-r1.2 \
                       > logs/unfiltered/rsem-prep.log 2>&1

# Align the RNA-Seq reads, calculate expression estimates
mkdir -p abundances/unfiltered tmp/unfiltered
for caste in q w
do
  for rep in {1..6}
  do
    sample=${caste}${rep}
    rsem-calculate-expression --paired-end \
                              --temporary-folder=tmp/unfiltered/rsem-calc-${sample}-temp \
                              --time \
                              --num-threads=32 \
                              rnaseq-clean/pdom-rnaseq-${sample}-trim-paired-1.fq \
                              rnaseq-clean/pdom-rnaseq-${sample}-trim-paired-2.fq \
                              iloci/pdom-loci-r1.2 \
                              abundances/unfiltered/pdom-${sample} \
                              > logs/unfiltered/rsem-calc-${sample}.log 2>&1
  done
  wait
done

# Compile expression data matrix: 1 row per iLocus, 1 column per sample
rsem-generate-data-matrix abundances/unfiltered/pdom-q1.genes.results \
                          abundances/unfiltered/pdom-q2.genes.results \
                          abundances/unfiltered/pdom-q3.genes.results \
                          abundances/unfiltered/pdom-q4.genes.results \
                          abundances/unfiltered/pdom-q5.genes.results \
                          abundances/unfiltered/pdom-q6.genes.results \
                          abundances/unfiltered/pdom-w1.genes.results \
                          abundances/unfiltered/pdom-w2.genes.results \
                          abundances/unfiltered/pdom-w3.genes.results \
                          abundances/unfiltered/pdom-w4.genes.results \
                          abundances/unfiltered/pdom-w5.genes.results \
                          abundances/unfiltered/pdom-w6.genes.results \
                          > abundances/unfiltered/pdom.genes.results

# Determine differentially expressed genes, FDR < 0.05
rsem-run-ebseq abundances/unfiltered/pdom.genes.results \
               6,6 pdom.genes.diffexp > logs/unfiltered/ebseq.log 2>&1
rsem-control-fdr pdom.genes.diffexp \
                 0.05 pdom.genes.diffexp.sig.05

# Generate a heatmap to visualize DEGs
git clone https://github.com/standage/tk-rnaseq.git
tk-rnaseq/de-viz ebseqinfile=abundances/unfiltered/pdom.genes.results \
                 ebseqoutfile=pdom.genes.diffexp \
                 workdir=expression-data \
                 samples="Q1,Q2,Q3,Q4,Q5,Q6,W1,W2,W3,W4,W5,W6" \
                 all
```

In the second stage, we filtered the iLoci to remove outliers and sequences that exhibited too much inconsistency across replicated, and then re-ran the entire analysis workflow.

```bash
# Build a table containing various data about each iLocus, most notably:
# the fold change; the number of raw reads mapping for each sample; and
# the unnormalized and normalized expression estimates for each sample.
make -f tk-rnaseq/build-table -j 16 \
     ebseqinfile=abundances/unfiltered/pdom.genes.results \
     ebseqoutfile=pdom.genes.diffexp \
     bamfilepattern=abundances/unfiltered/pdom-*.transcript.sorted.bam \
     fasta=iloci/pdom-loci-r1.2.fa \
     workdir=expression-data \
     all

# Filter the iLoci based on the specified criteria.
tk-rnaseq/de-filter --density=0.01,10 \
                    --nonzeros=5 \
                    --varcoef=1.0 \
                    --noheader \
                    --numsamples=6,6 \
                    < expression-data/expression-data.txt \
                    > expression-data/expression-data-filtered.txt

# Create a new fasta file containing only the iLoci that satisfy
# the filtering criteria.
tk-rnaseq/select-seq <(cut -f 1 expression-data/expression-data-filtered.txt) \
                     iloci/pdom-loci-r1.2.fa \
                     > iloci/pdom-loci-r1.2-filtered.fa

# Prepare/index the filtered iLocus sequences
mkdir -p logs/filtered
rsem-prepare-reference iloci/pdom-loci-r1.2-filtered.fa \
                       iloci/pdom-loci-r1.2-filtered \
                       > logs/filtered/rsem-prep.log 2>&1

# Align the RNA-Seq reads, calculate expression estimates
mkdir -p abundances/filtered tmp/filtered
for caste in q w
do
  for rep in {1..6}
  do
    sample=${caste}${rep}
    rsem-calculate-expression --paired-end \
                              --temporary-folder=tmp/filtered/rsem-calc-${sample}-temp \
                              --time \
                              --num-threads=48 \
                              rnaseq-clean/pdom-rnaseq-${sample}-trim-paired-1.fq \
                              rnaseq-clean/pdom-rnaseq-${sample}-trim-paired-2.fq \
                              iloci/pdom-loci-r1.2-filtered \
                              abundances/filtered/pdom-${sample} \
                              > logs/filtered/rsem-calc-${sample}.log 2>&1
  done
  wait
done

# Compile expression data matrix
rsem-generate-data-matrix abundances/filtered/pdom-q1.genes.results \
                          abundances/filtered/pdom-q2.genes.results \
                          abundances/filtered/pdom-q3.genes.results \
                          abundances/filtered/pdom-q4.genes.results \
                          abundances/filtered/pdom-q5.genes.results \
                          abundances/filtered/pdom-q6.genes.results \
                          abundances/filtered/pdom-w1.genes.results \
                          abundances/filtered/pdom-w2.genes.results \
                          abundances/filtered/pdom-w3.genes.results \
                          abundances/filtered/pdom-w4.genes.results \
                          abundances/filtered/pdom-w5.genes.results \
                          abundances/filtered/pdom-w6.genes.results \
                          > abundances/filtered/pdom.filtered.genes.alllibs.results

# Determine DEGs
rsem-run-ebseq abundances/filtered/pdom.filtered.genes.alllibs.results \
               6,6 pdom.filtered.genes.alllibs.diffexp > logs/filtered/ebseq-alllibs.log 2>&1
rsem-control-fdr pdom.filtered.genes.alllibs.diffexp \
                 0.05 pdom.filtered.genes.alllibs.diffexp.sig.05

# Visualize DEGs
tk-rnaseq/de-viz ebseqinfile=abundances/filtered/pdom.filtered.genes.alllibs.results \
                 ebseqoutfile=pdom.filtered.genes.alllibs.diffexp \
                 workdir=expression-data-filtered \
                 samples="Q1,Q2,Q3,Q4,Q5,Q6,W1,W2,W3,W4,W5,W6" \
                 all
```

In the final stage, we identified the least characteristic replicate from each condition (caste), removed those replicates, and re-ran the final step of the analysis workflow.

```bash
# Determine which replicate from each caste is the most similar to the other caste
tk-rnaseq/bully --samples="Q1,Q2,Q3,Q4,Q5,Q6,W1,W2,W3,W4,W5,W6" --numreps=6,6 \
                --norm=abundances/filtered/pdom.filtered.genes.alllibs.results \
                < expression-data-filtered/de.expr.dat.sorted

# Compile expression data matrix, this time sans outlier replicates
rsem-generate-data-matrix abundances/filtered/pdom-q1.genes.results \
                          abundances/filtered/pdom-q2.genes.results \
                          abundances/filtered/pdom-q3.genes.results \
                          abundances/filtered/pdom-q5.genes.results \
                          abundances/filtered/pdom-q6.genes.results \
                          abundances/filtered/pdom-w1.genes.results \
                          abundances/filtered/pdom-w2.genes.results \
                          abundances/filtered/pdom-w3.genes.results \
                          abundances/filtered/pdom-w4.genes.results \
                          abundances/filtered/pdom-w5.genes.results \
                          > abundances/filtered/pdom.filtered.genes.results

# Determine DEGs
rsem-run-ebseq abundances/filtered/pdom.filtered.genes.results \
               5,5 pdom.filtered.genes.diffexp > logs/filtered/ebseq.log 2>&1
rsem-control-fdr pdom.filtered.genes.diffexp \
                 0.05 pdom.filtered.genes.diffexp.sig.05

# Visualize DEGs
tk-rnaseq/de-viz ebseqinfile=abundances/filtered/pdom.filtered.genes.results \
                 ebseqoutfile=pdom.filtered.genes.diffexp \
                 workdir=expression-data-filtered-final \
                 samples="Q1,Q2,Q3,Q5,Q6,W1,W2,W3,W4,W5" \
                 all
```

## References

- **Li B, Dewey CN** (2011) RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. *BMC Bioinformatics*, **12**:323, [doi:10.1186/1471-2105-12-323](http://dx.doi.org/10.1186/1471-2105-12-323).
- **Langmead B, Trapnell C, Pop M, Salzberg SL** (2009) Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. *Genome Biology*, **10**:R25, [doi:10.1186/gb-2009-10-3-r25](http://dx.doi.org/10.1186/gb-2009-10-3-r25).
- **Leng N, Dawson JA, Thomson JA, Ruotti V, Rissman AI, Smits BMG, Haag JD, Gould MN, Stewart RN, Kendziorski C** (2013) EBSeq: an empirical Bayes hierarchical model for inference in RNA-seq experiments. *Bioinformatics*, **29**(8):1035-1043, [doi:10.1093/bioinformatics/btt087](http://dx.doi.org/10.1093/bioinformatics/btt087).

[RSEM]: http://deweylab.biostat.wisc.edu/rsem
[1.2.7]: https://github.com/bli25wisc/RSEM/archive/v1.2.7.tar.gz
[Bowtie]: http://bowtie-bio.sourceforge.net/index.shtml
[1.0.0]: https://github.com/BenLangmead/bowtie/archive/v1.0.0.tar.gz
[EBSeq]: https://www.biostat.wisc.edu/~kendzior/EBSEQ/
[tk-rnaseq toolkit]: https://github.com/standage/tk-rnaseq
[5ce0f8c075]: https://github.com/standage/tk-rnaseq/tree/5ce0f8c075dcc2971e7b120e7c1886c299291d8c
