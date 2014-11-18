#!/usr/bin/env python
import re, sys

isogroups = {}
for line in sys.stdin:
  if not line.startswith("PdomMRNA"):
    continue

  alnmatch = re.search("(\S+)\s+(\d+)", line)
  if not alnmatch:
    continue

  qid, sid = alnmatch.group(1), alnmatch.group(2)
  mid = qid[:-2]
  if mid not in isogroups:
    isogroups[mid] = {}
  isogroups[mid][qid] = sid

for mid in isogroups:
  mid_0 = mid + ".0"
  if mid_0 not in isogroups[mid]:
    continue

  sub_0 = isogroups[mid][mid_0]
  for qid in isogroups[mid]:
    sub_k = isogroups[mid][qid]
    if sub_k != sub_0:
      print "%s=%s\t%s=%s" % (mid_0, sub_0, qid, sub_k)
