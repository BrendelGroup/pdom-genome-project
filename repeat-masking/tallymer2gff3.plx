#!/usr/bin/perl -w

use strict;
use Getopt::Long;

## Define some paramters
my $k = 19;
my $min_len_thresh = 50;

my $merge_thresh = $k * 2;

my $seq_file;



## See code to understand
die "some repeats will be missed\n"
    if $min_len_thresh < $k;



## Get options
GetOptions
  (
   "k=i"    => \$k,
   "min=i"  => \$min_len_thresh,
   "file|seq=s" => \$seq_file,
  )
  or die "failed to parse command line\n";

die "pass a sequence file\n"
  unless -s $seq_file;



## Read in sequence names
warn "reading in sequence names from $seq_file\n";
my @seq_names;

open NAMES, '<', $seq_file
  or die "failed to open $seq_file\n";

while(<NAMES>){
  next unless /^>/; chomp;
  push @seq_names, substr($_, 1);
}

warn "read ", scalar(@seq_names), "\n";
warn "OK\n";





## Parse the output of the tallymer search

my ($p_seq, $p_pos, $len, $depth) = (-1, -1, -1);

while(<>){
  # ignore commets
  next if /^#/;
  
  chomp;
  my ($seq, $pos, $x, undef) = split;
  
  ## are we beginning a new sequence?
  if($seq != $p_seq){
    #warn $seq, "\n";
    
    # we should dump the current region
    &print_region( $p_seq, $p_pos, $len, $depth / $len )
      if $len >= $min_len_thresh;
    
    # and reset our counters
    $len = 0;
    $depth = 0;
    $p_pos = $pos;
    $p_seq = $seq;
    next;
  }
  
  ## if not, are we continuing a region?
  if($pos - $p_pos <= $merge_thresh){
    $len += $pos - $p_pos;
    $depth += $x;
    $p_pos = $pos;
    next;
  }
  
  ## if not, we must be leaving a region...
  
  # we should dump the current region
  &print_region( $seq, $p_pos, $len, $depth / $len )
    if $len >= $min_len_thresh;
  
  # and reset our counters
  $len = 0;
  $depth = 0;
  $p_pos = $pos;
}

warn "OK\n";



sub print_region () {
  $p_seq = shift;
  $p_pos = shift;
  $len   = shift;
  $depth = shift;
  
  our $j;
  
  print
    join("\t",
	 $seq_names[$p_seq],
	 'tallymer',
	 'repeat_region',
	 $p_pos - $len,
	 $p_pos + $k,
	 sprintf("%d", $depth),
	 '+',
	 '.',
	 join(';',
	      'ID='. sprintf("%08d", $j++),
	      'dbxref=SO:0000657'
	     ),
	), "\n";
  
  return 1;
}
