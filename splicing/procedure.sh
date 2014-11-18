#!/usr/bin/env bash
set -e

prep_genome()
{
  local WD=$1

  mkdir -p $WD/genome
  IPLANT=/iplant/home/standage/Polistes_dominula/r1.2/interval-loci
  iget $IPLANT/pdom-loci-r1.2-filtered.fa $WD/genome/
  iget $IPLANT/pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 $WD/genome/

  cd $WD/genome
  bowtie2-build pdom-loci-r1.2-filtered.fa pdom-loci-r1.2 \
      > $WD/logs/pdom-loci-r1.2.index.log 2>&1
  cd -
}

#-------------------------------------------------------------------------------
# Align RNA-Seq reads using Tophat (2.0.12; Bowtie2 2.2.3)
align_reads()
{
  local WD=$1
  local NP=$2

  mkdir -p $WD/alignments
  for caste in q w
  do
    for rep in {1..6}
    do
      sample=${caste}${rep}
      tophat --output-dir $WD/alignments/${sample} \
             --GTF $WD/genome/pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 \
             --mate-inner-dist 70 \
             --min-anchor 6 \
             --num-threads $NP \
             $WD/genome/pdom-loci-r1.2 \
             $WD/rnaseq-clean/pdom-rnaseq-${sample}-trim-paired-1.fq \
             $WD/rnaseq-clean/pdom-rnaseq-${sample}-trim-paired-2.fq \
             > $WD/logs/${sample}.tophat.log 2>&1
    done
  done
}


#-------------------------------------------------------------------------------
# Assemble read alignments using Cufflinks (2.2.1)
assemble_reads()
{
  local WD=$1
  local NP=$2

  mkdir -p $WD/assemblies
  rm -f $WD/assemblies/all.txt
  for caste in q w
  do
    for rep in {1..6}
    do
      sample=${caste}${rep}
      echo $WD/assemblies/${sample}/transcripts.gtf >> $WD/assemblies/all.txt

      cufflinks --output-dir $WD/assemblies/${sample} \
                --num-threads $NP \
                --multi-read-correct \
                --GTF-guide $WD/genome/pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 \
                $WD/alignments/${sample}/accepted_hits.bam \
                > $WD/logs/${sample}.cufflinks.log 2>&1
    done
  done

  cuffmerge -o $WD/assemblies/merged \
            --num-threads $NP \
            --ref-sequence $WD/genome/pdom-loci-r1.2-filtered.fa \
            --ref-gtf $WD/genome/pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 \
            $WD/assemblies/all.txt \
            > $WD/logs/cuffmerge.log 2>&1

  # Convert final assembly into GFF3 format, remove strandless features
  gt gtf_to_gff3 -tidy assemblies/merged/merged.gtf 2> $WD/logs/gtf2gff.log \
      | perl -ne '@f = split(/\t/); next if($f[6] eq "."); print' \
      | gt gff3 -sort -tidy -o pdom-annot-r1.2-tuxedo.gff3 -force
}


#-------------------------------------------------------------------------------
# Catalog alternative splicing events
alt_splice()
{
  local WD=$1
  asinspect --refr $WD/genome/pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 \
            --gtf \
            --out pdom-as.tsv \
            $WD/assemblies/merged/merged.gtf \

  python coding-cassette-exons.py \
         pdom-as.tsv pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 \
         > pdom-as-ce-coding.tsv

  python ce-isoforms.py \
         pdom-as-ce-coding.tsv pdom-annot-r1.2-iloci-sanstrna-filtered.gff3 \
      | gt gff3 -retainids -sort -tidy \
      > pdom-annot-r1.2-withceisoforms.gff3

  gt extractfeat -matchdescstart -join -retainids -translate \
                 -seqfile pdom-loci-r1.2.fa -type CDS \
                 pdom-annot-r1.2-withceisoforms.gff3 \
      | perl get_ce_iso.pl \
      > pdom-annot-r1.2-ce-isoforms.fa

  curl -O ftp://ftp.ncbi.nih.gov/genomes/Apis_mellifera/protein/protein.fa.gz
  gunzip protein.fa.gz
  mv protein.fa amel-ncbi-prot.fa
  makeblastdb -in amel-ncbi-prot.fa -dbytype prot -parse_seqids
  blastp -db amel-ncbi-prot.fa -query pdom-annot-r1.2-ce-isoforms.fa \
         -evalue 1e-6 -num_threads 32 -out pdom-vs-amel.blastp

  MuSeqBox -i pdom-vs-amel.blastp -n 1 -s 1 -l 25 -d 4 -L 128 -c ce-crtfile \
      > pdom-vs-amel-blastp.msb

  python ce-filter.py < pdom-vs-amel-blastp.msb > pdom-vs-amel-conserved.txt
  cut -f 1 -d '=' pdom-vs-amel-conserved.txt | sort | uniq \
      > pdom-annot-r1.2-ce-list.txt
}

#-------------------------------------------------------------------------------
# Perform differential splicing analysis
diff_splice()
{
  local WD=$1
  local NP=$2

  cuffdiff --output-dir $WD/diff \
           --labels Queen,Worker \
           --num-threads $NP \
           --multi-read-correct \
           $WD/assemblies/merged/merged.gtf \
           $WD/alignments/q1/accepted_hits.bam,$WD/alignments/q2/accepted_hits.bam,$WD/alignments/q3/accepted_hits.bam,$WD/alignments/q4/accepted_hits.bam,$WD/alignments/q5/accepted_hits.bam,$WD/alignments/q6/accepted_hits.bam \
           $WD/alignments/w1/accepted_hits.bam,$WD/alignments/w2/accepted_hits.bam,$WD/alignments/w3/accepted_hits.bam,$WD/alignments/w4/accepted_hits.bam,$WD/alignments/w5/accepted_hits.bam,$WD/alignments/w6/accepted_hits.bam \
           > $WD/logs/cuffdiff.log 2>&1
}


#-------------------------------------------------------------------------------
# Main method
main()
{
  local WD=$1
  if [ -z "$WD" ]; then
    WD=.
  fi
  local NP=$2
  if [ -z "$NP" ]; then
    NP=32
  fi

  mkdir -p $WD/logs

  prep_genome    $WD
  align_reads    $WD $NP
  assemble_reads $WD $NP
  alt_splice     $WD
  diff_splice    $WD $NP
}

main "$@"
