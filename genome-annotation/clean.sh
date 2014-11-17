#!/usr/bin/env bash
set -euo pipefail
perl -ne 'last if(m/^##FASTA/); print' < $1 \
    | egrep -v $'\t(contig|expressed_sequence_match|match|match_part|protein_match|translated_nucleotide_match)\t'
