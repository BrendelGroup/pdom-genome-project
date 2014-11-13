#!/usr/bin/env bash
set -euo pipefail

if [ -z "${1:-}" ]; then
  NumThreads=1
else
  NumThreads=$1
fi

iget -Vr /iplant/home/standage/Polistes_dominula/sequence/genome
ls genome/*.gz | parallel --gnu --jobs $NumThreads gunzip

FastqFiles=$(ls genome/*.fq)
for k in  17 21 25 29
do
  jellyfish count -m $k -s 100M -t $NumThreads -C -o pdom-${k}mers.jf $FastqFiles
  jellyfish histo pdom-${k}mers.jf > pdom-${k}mers.hist
done

./size-coverage-estimate.R

if [ ! -z "${2:-}" ]; then
  if [ "$2" != "keep" ]; then
    rm -r genome/ *.jf
  fi
fi

