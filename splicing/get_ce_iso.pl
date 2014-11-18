#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my $loader = Bio::SeqIO->new(-fh => \*STDIN,  -format => "Fasta");
my $writer = Bio::SeqIO->new(-fh => \*STDOUT, -format => "Fasta");
while(my $seq = $loader->next_seq())
{
  if($seq->id() =~ m/PdomMRNAr1.2-\d+\.\d+\.\d+/)
  {
    my $s = $seq->seq();
    if($s =~ m/\*/)
    {
      my @aa = split(/\*+/, $s);
      $seq->seq($aa[0]);
    }
    $writer->write_seq($seq);
  }
}
