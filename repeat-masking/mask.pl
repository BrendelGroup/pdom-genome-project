#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my $usage = "perl mask.pl rep.gff exlusions.txt seq.fa > seq.masked.fa";
my $gff_file = shift(@ARGV) or die("Usage: $usage");
my $exc_file = shift(@ARGV) or die("Usage: $usage");
my $seq_file = shift(@ARGV) or die("Usage: $usage");

my %exclusions;
open(my $EXC, "<", $exc_file) or die("unable to open '$exc_file'");
while(my $line = <$EXC>)
{
  chomp($line);
  my($seqid, $start, $end) = split(/\t/, $line);
  $exclusions{"${seqid}_${start}-${end}"} = 1;
}
close($EXC);

my %reps;
open(my $GFF, "<", $gff_file) or die("unable to open '$gff_file'");
while(my $line = <$GFF>)
{
  my @v = split(/\t/, $line);
  next if(scalar(@v) != 9);
  my($seqid, $start, $end) = ($v[0], $v[3], $v[4]);
  if($exclusions{"${seqid}_${start}-${end}"})
  {
    print STDERR "skip ${seqid}_${start}-${end}\n";
    next;
  }

  $reps{$seqid} = [] unless($reps{$seqid});
  push(@{$reps{$seqid}}, [$start, $end]);
}
close($GFF);

my $reader = Bio::SeqIO->new( -file => $seq_file, -format => "Fasta" );
my $writer = Bio::SeqIO->new( -fh => \*STDOUT,    -format => "Fasta" );
while(my $sequence = $reader->next_seq)
{
  my $seq = $sequence->seq;
  for my $rep(@{$reps{$sequence->id}})
  {
    my($start, $end) = @{$rep};
    my $length = $end - $start + 1;
    substr($seq, $start, $length, "N" x $length);
  }
  $sequence->seq($seq);
  $writer->write_seq($sequence);
}
$reader->close();
$writer->close();

