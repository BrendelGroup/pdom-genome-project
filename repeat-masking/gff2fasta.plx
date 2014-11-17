#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Bio::SeqIO;

## Define some paramters
my $verbose = 0;

my $seq_file;
my $feature_file;



## Get options
GetOptions
  (
   "seq_file=s"     => \$seq_file,
   "feature_file=s" => \$feature_file,
  )
  or die "failed to parse command line\n";

die "pass a sequence file\n"
  unless -s $seq_file;

die "pass a GFF file\n"
  unless -s $feature_file;



## read in gff

warn "reading GFF\n";

my %gff;

open GFF, '<', $feature_file
  or die "fail\n";

while(<GFF>){
  my ($seqid, undef, undef, $start, $end,
      undef, undef, undef, $attrs) = split;
  
  push @{$gff{$seqid}}, [$start, $end, $attrs];
}

warn "OK\n";



## Do the fasta

my $seqio = Bio::SeqIO->
  new( -file => $seq_file,
       -format => 'fasta' )
  or die "double fail\n";

while(my $sobj = $seqio->next_seq){
  my $seqid = $sobj->id;
  
  unless(defined($gff{$seqid})){
    warn "no repeats for $seqid\n"
      if $verbose > 0;
    next;
  }
  
  my @repeats = @{$gff{$seqid}};
  
  my $seq = $sobj->seq;
  
  for(@repeats){
    my ($start, $end, $attrs) = @$_;
    
    warn join("\t", $start, $end, $attrs), "\n"
      if $verbose > 1;
    
    my %attrs = split(/=|;/, $attrs);
    
    print ">$seqid-". $attrs{"ID"}.
      "/$start-$end (". ($end-$start+1). ")\n";
    
    print substr($seq, $start, $end-$start+1), "\n";
  }
  
  #exit;
}

warn "OK\n";
