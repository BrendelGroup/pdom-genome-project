# Transcript assembly

Raw DNA-Seq reads were groomed using [Trimmomatic][] version [0.22][], and the groomed reads were then assembled using [Trinity][] version [r20131110][].
Then the assembled transcripts were processed with [mRNAmarkup][] version [10-3-2013][] to remove contaminants and correct erroneously assembled chimeric transcripts.
The cleaned and annotated transcripts have been deposited in the Pdom Data Store at `r1.2/transcript-assembly/`.

## Procedure (interactive)

### Short read quality control

First, designate the number of available processors to speed up Trimmomatic's computations.
Also, provide the path of the `trimmomatic-0.22.jar` file contained in the Trimmomatic source code distribution.

```bash
NumThreads=16
TrimJar=/usr/local/src/Trimmomatic-0.22/trimmomatic-0.22.jar
PdomData=/iplant/home/standage/Polistes_dominula
```

Now for the processing.
We apply the following filters to each read pair.

  - remove adapter contamination
  - remove any nucleotides at either end of the read whose quality score is below 3
  - trim the read once the average quality in a 5bp sliding window falls below 20
  - discard any reads which, after processing, fall below the 40bp length threshold

```bash
for caste in q w
do
  for rep in {1..6}
  do
    sample=${caste}${rep}
    iget ${PdomData}/sequence/transcriptome/pdom-rnaseq-${sample}-1.fq.gz
    iget ${PdomData}/sequence/transcriptome/pdom-rnaseq-${sample}-2.fq.gz
    ./run-trim.sh $sample $TrimJar $NumThreads
  done
done
```

### Assembly with Trinity

Trinity requires a single input file---or a pair of input files for paired-end data.
We need to combine all of the data into a single pair of files.

```bash
cat pdom-rnaseq-*-trim-1.fq > pdom-rnaseq-all-trim-1.fq
cat pdom-rnaseq-*-trim-2.fq > pdom-rnaseq-all-trim-2.fq
```

We'll then execute the Trinity assember using the `--CuffFly` reconstruction algorithm.

```bash
Trinity.pl --seqType fq \
           --JM 100G \
           --bflyHeapSpaceMax 50G \
           --output pdom-trinity \
           --CPU $NumThreads \
           --left  pdom-rnaseq-all-trim-1.fq \
           --right pdom-rnaseq-all-trim-2.fq \
           --full_cleanup \
           --jaccard_clip \
           --CuffFly
```

### Post-processing with mRNAmarkup

Contaminant, reference protein, and miRNA databases were collected as described in the mRNAmarkup documentation (`db/0README` and `db/0README-hy`).
The mRNAmarkup procedure was then run on the Trinity output.
Be sure to edit the `mRNAmarkup.conf` file with the correct paths to the databases.

```bash
mRNAmarkup -c mRNAmarkup.conf \
           -i pdom-trinity/Trinity.fasta \
           -o output-mRNAmarkup
```

### Potential *Polistes*-specific genes

*Polistes metricus* and *P. canadensis* transcript assemblies were groomed and annotated using the same mRNAmarkup procedure.
All three transcript sets had a substantial number of TSAs that could not be annotated by mRNAmarkup.
The longest open reading frame translations of at least 80 amino acids derived from these unmatched TSAs were pairwise compared with BLASTp to detect any *Polistes*-conserved and species-specific genes.
This procedure relies on a workflow and several scripts available in the supplemental/ directory of the documentation repository.

The workflow to find the *Polistes*-conserved transcripts with no external protein matches can be executed with this command.
The workflow depends on the `dnatopro` program (part of the VBTools utilities included in the mRNAmarkup distribution) and the NCBI BLAST+ suite.
```bash
./venn.make BLASTTHREADS=$NumThreads all
```

Transcripts associated with potential *Polistes*-specific genes will be placed in the `pdom-tsa-r1.2-unmatched-pep.fa` file.

![Otherwise unmatched TSAs with significant pairwise protein-level similarity among the *Polistes*. For each intersection, only the smallest number is shown from the species numbers for the comparison. For example, the 3-species intersection consists of 144 TSAs from *P. dominula*, 136 TSAs from *P. canadensis*, and 95 TSAs from *P. metricus*. Multiply matched TSA products within a species would include products from alternative transcript isoforms and from duplicated genes.](pdom-venn.png)

## Procedure (automated)

The same procedure can also be run in batch mode using the following commands (in the `transcript-assembly` directory).
Note that this procedure does not automate the mRNAmarkup analysis.
After producing the Trinity assembly, it retrieves previously computed mRNAmarkup results for the clade-specific gene analysis.

```bash
make NumThreads=16 \
     TrimJar=/usr/local/src/Trimmomatic-0.22/trimmomatic-0.22.jar
make clean
```

## References

- **Lohse M, Bolger AM, Nagel A, Fernie AR, Lunn JE, Stitt M, Usadel B** (2012) RobiNA: a user-friendly, integrated software solution for RNA-Seq-based transcriptomics. *Nucleic Acids Research*, **40**:W622-7, [doi:10.1093/nar/gks540](http://dx.doi.org/10.1093/nar/gks540).
- **Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A** (2011) Full-length transcriptome assembly from RNA-seq data without a reference genome. *Nature Biotechnology*, **29**(7):644-52, [doi:10.1038/nbt.1883](http://dx.doi.org/10.1038/nbt.1883).

[Trimmomatic]: http://www.usadellab.org/cms/index.php?page=trimmomatic
[0.22]: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.22.zip
[Trinity]: http://trinityrnaseq.sourceforge.net/
[r20131110]: http://downloads.sourceforge.net/project/trinityrnaseq/trinityrnaseq_r20131110.tar.gz
[mRNAmarkup]: http://brendelgroup.org/bioinformatics2go/mRNAmarkup.php
[10-3-2013]: http://www.brendelgroup.org/bioinformatics2go/Download/mRNAmarkup-10-3-2013.tar.gz
