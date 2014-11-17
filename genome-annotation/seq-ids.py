#!/usr/bin/env python
import re, sys

if __name__ == "__main__":
  usage = ("\nReplace sequence IDs using the given correspondence of new IDs to old IDs.\n"
           "Usage: python seq-ids.py mapping.txt < seqs.fasta > seqs-fixed.fasta\n")

  if len(sys.argv) != 2:
    print >> sys.stderr, usage

  mapping = {}
  for line in open(sys.argv[1], "r"):
    line = line.rstrip()
    newid, oldid = line.split("\t")
    mapping[oldid] = newid

  for line in sys.stdin:
    line = line.rstrip()
    seqidmatch = re.search("^>(\S+)", line)
    if seqidmatch:
      seqid = seqidmatch.group(1)
      newid = mapping[seqid]
      line = re.sub("^>(\S+)", r">%s \1" % newid, line)
    print line
