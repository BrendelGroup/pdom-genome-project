#!/usr/bin/env perl
use strict;

my %names;
open(my $NAMES, "<", $ARGV[0]) or die("fail");
while(<$NAMES>)
{
  if(m/^(\d+)\t\|\t([^\t]+)\t\|.+\|\tscientific name\t/)
  {
    $names{$2} = $1
  }
}
close($NAMES);

my $findh = 0;
my $qid = "";
while(<STDIN>)
{
  if(m/<Iteration_query-def>(\S+)/)
  {
    $findh = 1;
    $qid = $1;
    next;
  }

  next if(not $findh);
  if(m/<Hit_def>.+\[([^\]]+)\]/)
  {
    $findh = 0;
    my $gi = $names{$1};
#print STDERR "DEBUG species=$1, gi=$gi\n";
    print "$qid\t$gi\t$1\n";
  }
}
