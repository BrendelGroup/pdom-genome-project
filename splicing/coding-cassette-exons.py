#!/usr/bin/env python
import sys

if __name__ == "__main__":
  ce_file = sys.argv[1]
  ft_file = sys.argv[2]

  cds = {}
  with open(ft_file) as gff3:
    for line in gff3:
      fields = line.split()
      if len(fields) == 9 and fields[2] == "CDS":
        label = "%s_%s-%s" % (fields[0], fields[3], fields[4])
        cds[label] = 1

  with open(ce_file, "r") as ce_list:
    print ce_list.readline().rstrip()
    for line in ce_list:
      if not line.startswith("CE\t"):
        continue
      fields = line.rstrip().split("\t")
      seqid = fields[4]
      coords = fields[5].split(",")[1]
      ce = "%s_%s" % (seqid, coords)
      if ce in cds:
        print line.rstrip()
