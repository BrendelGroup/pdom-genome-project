#!/usr/bin/env python
import sys

def parse_fasta(fp):
  """
  Generic fasta parsing function
  Stolen shamelessly from http://stackoverflow.com/a/7655072/459780
  """
  name, seq = None, []
  for line in fp:
    line = line.rstrip()
    if line.startswith(">"):
      if name: yield (name, ''.join(seq))
      name, seq = line, []
    else:
      seq.append(line)
  if name: yield (name, ''.join(seq))

if __name__ == "__main__":
  print "##gff-version   3"
  with open(sys.argv[1], "r") as fp:
    for defline, seq in parse_fasta(fp):
      seqid = defline[1:].split(" ")[0]
      seqlength = len(seq)
      print "##sequence-region   %s 1 %d" % (seqid, seqlength)

  for line in sys.stdin:
    if line.startswith("##gff-version") or line.startswith("##sequence-region"):
      continue
    print line.rstrip()
