#!/usr/bin/env bash
set -euo pipefail

Sample=$1
TrimJar=$2
NumThreads=$3

# The '|| true' at the end of this line is simply because awk gives a non-zero return status on success, causing the script to prematurely terminate.
ReadLength=$(zcat pdom-gdnaseq-${Sample}-1.fq.gz | head -n 2 | tail -n 1 | awk '{ print length($0) }' || true)
if [ $ReadLength == 100 ]; then
  MinLength=40
elif [ $ReadLength == 35 ]; then
  MinLength=26
else
  echo "Error: read length $ReadLength"
  exit 1
fi

java -classpath $TrimJar org.usadellab.trimmomatic.TrimmomaticPE \
     -threads $NumThreads \
     -phred33 \
     pdom-gdnaseq-${Sample}-1.fq.gz \
     pdom-gdnaseq-${Sample}-2.fq.gz \
     pdom-gdnaseq-${Sample}-trim-1.fq \
     pdom-gdnaseq-${Sample}-trim-unpaired-1.fq \
     pdom-gdnaseq-${Sample}-trim-2.fq \
     pdom-gdnaseq-${Sample}-trim-unpaired-2.fq \
     ILLUMINACLIP:illumina-std-adapters.fa:2:40:15 \
     LEADING:3 \
     TRAILING:3 \
     SLIDINGWINDOW:5:20 \
     MINLEN:$MinLength
