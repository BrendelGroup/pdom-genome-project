#!/usr/bin/env python
import re, sys

idfile  = sys.argv[1]
alnfile = sys.argv[2]

ids2keep = {}
with open(idfile, "r") as fp:
  for line in fp:
    ids2keep[line.rstrip()] = 1

with open(alnfile, "r") as fp:
  for line in fp:
    match = re.search("Name=(PdomTSAr1.2-\d+)", line)
    if match:
      tsaid = match.group(1)
      if tsaid in ids2keep:
        print line.rstrip()
