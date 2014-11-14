#!/usr/bin/env python
import re, sys

def parse_fasta(fp):
  """
  Stolen shamelessly from http://stackoverflow.com/a/7655072/459780.
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

if len(sys.argv) != 4:
  print "Usage: python id-resolve.py vID-PdPcPm PdPRT pdom-tsa-r1.2.fa"
  exit(0)

vidfile = sys.argv[1]
prtfile = sys.argv[2]
tsafile = sys.argv[3]

vids2keep = {}
with open(vidfile, "r") as fp:
  for line in fp:
    vids2keep[line.rstrip()] = 1

cmpids2prtids = {}
with open(prtfile, "r") as fp:
  for defline, seq in parse_fasta(fp):
    seqidmatch = re.search(">lcl\|(\S+) (\S+)", defline)
    assert seqidmatch, "Issue parsing '%s'" % defline
    prtid = seqidmatch.group(1)
    cmpid = seqidmatch.group(2)
    cmpids2prtids[cmpid] = prtid

with open(tsafile, "r") as fp:
  for line in fp:
    if not line.startswith(">"):
      continue
    idmatch = re.search(">(\S+) (\S+)", line)
    assert idmatch, "Issue parsing '%s'" % line.rstrip()
    tsaid = idmatch.group(1)
    cmpid = idmatch.group(2)
    if cmpid not in cmpids2prtids:
      continue
    prtid = cmpids2prtids[cmpid]
    if prtid in vids2keep:
      print tsaid
