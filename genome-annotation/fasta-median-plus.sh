#!/usr/bin/env bash
set -eu
MEDIAN=$(perl -MBio::SeqIO -e '$l = Bio::SeqIO->new(-fh=>\*STDIN, -format=>"Fasta"); while($s = $l->next_seq){ printf "%s\t%lu\n", $s->id, $s->length }' < $1 | ./upper-median.R)
perl -MBio::SeqIO -e "\$l = Bio::SeqIO->new(-fh=>\*STDIN, -format=>'Fasta'); \$w = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta'); while(\$s = \$l->next_seq){ \$w->write_seq(\$s) if(\$s->length >= $MEDIAN) }" < $2
