#!/usr/bin/env python
import re, sys

def process_isoforms(mrna_id, mrna_feat, mrna_subfeats, cassette_exons):
  assert mrna_id in cassette_exons
  ce_count = 0
  print mrna_feat.replace(mrna_id, mrna_id + ".0")
  for feat in mrna_subfeats:
    print feat.replace(mrna_id, mrna_id + ".0")

  for ce in cassette_exons[mrna_id]:
    ce_count += 1
    new_id = "%s.%d" % (mrna_id, ce_count)
    print mrna_feat.replace(mrna_id, new_id) + ";skipped_exon=" + ce
    for feat in mrna_subfeats:
      fields = feat.split("\t")
      coords = "%s-%s" % (fields[3], fields[4])
      if coords == ce:
        continue
      print feat.replace(mrna_id, new_id)

if __name__ == "__main__":
  ce_file = sys.argv[1]
  ft_file = sys.argv[2]

  cassette_exons = {}
  with open(ce_file, "r") as fp:
    for line in fp:
      if not line.startswith("CE\t"):
        continue
      fields = line.rstrip().split("\t")
      assert len(fields) == 6, "expected 6 fields, found %d" % len(fields)
      mrna_id = fields[2]
      exon_coords = fields[5].split(",")[1]
      if mrna_id not in cassette_exons:
        cassette_exons[mrna_id] = {}
      cassette_exons[mrna_id][exon_coords] = 1
  
  with open(ft_file, "r") as fp:
    mrna_id = None
    mrna_feat = None
    mrna_subfeats = []
    for line in fp:
      line = line.rstrip()
      fields = line.split("\t")
      if len(fields) != 9:
        if line != "###":
          print line
        continue

      ftype = fields[2]
      assert ftype in ("gene", "mRNA", "exon", "CDS", "five_prime_UTR",
                       "three_prime_UTR"), "unrelated feature type %s" % ftype
      if ftype == "gene":
        if mrna_id:
          process_isoforms(mrna_id, mrna_feat, mrna_subfeats, cassette_exons)
          mrna_id = None
          ce_count = 0
          mrna_subfeats = []
        print line
      elif ftype == "mRNA":
        if mrna_id:
          process_isoforms(mrna_id, mrna_feat, mrna_subfeats, cassette_exons)
          mrna_id = None
          ce_count = 0
          mrna_subfeats = []
        mrna_id = re.search("ID=([^;\n]+)", fields[8]).group(1)
        if mrna_id in cassette_exons:
          mrna_feat = line
        else:
          mrna_id = None
          print line
      elif ftype in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
        parent_id = re.search("Parent=([^;\n]+)", fields[8]).group(1)
        if mrna_id:
          assert parent_id == mrna_id
          mrna_subfeats.append(line)
        else:
          print line
