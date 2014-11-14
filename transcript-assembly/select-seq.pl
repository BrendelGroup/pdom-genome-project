#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

# ------------------------------------------------------------------------------
# select-seq: extract sequences by ID from a Fasta file.
# ------------------------------------------------------------------------------

my $usage = "perl $0 seqidlist.txt seqs.fasta";
my $seqlist = shift(@ARGV) or die("Usage: $usage");
my $seqfile = shift(@ARGV) or die("Usage: $usage");

my %ids;
open(my $ID, "<", $seqlist) or die("error: unable to open $seqlist");
while(my $line = <$ID>)
{
  chomp($line);
  $ids{$line} = 1;
}

my $seqinstream  = Bio::SeqIO->new(-file => $seqfile, -format => "Fasta");
my $seqoutstream = Bio::SeqIO->new(-fh   => \*STDOUT, -format => "Fasta");
while(my $seq = $seqinstream->next_seq)
{
  my $molid = $seq->id;
  $seqoutstream->write_seq($seq) if($ids{$molid});
}

