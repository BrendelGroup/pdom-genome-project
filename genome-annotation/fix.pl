#!/usr/bin/env perl
next unless(/^PdomSCFr1.2/);
@v = split(/\t/);
if($v[3] eq "" or $v[4] eq "")
{
  $two = <STDIN>;
  @v2 = split(/\t/, $two);
  $v[3] = $v2[3];
  $v[4] = $v2[4];
  print join("\t", @v), $two;
}
else{ print }
