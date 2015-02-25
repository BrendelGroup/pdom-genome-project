import re, sys


class Entry():
  """
  Very simplistic representation of a GFF3 entry: a single line of a GFF3 file
  which may or may not be a complete representation of a genomic feature.
  """
  def __init__(self, line):
    self.fields = line.rstrip().split("\t")
    if not self.is_feature():
      self.line = line
      return
  
    self.seqid  = self.fields[0]
    self.source = self.fields[1]
    self.ftype  = self.fields[2]
    self.start  = self.fields[3]
    self.end    = self.fields[4]
    self.score  = self.fields[5]
    self.strand = self.fields[6]
    self.phase  = self.fields[7]
    self.attrs  = {}
    for keyvaluepair in self.fields[8].split(";"):
      key, value = keyvaluepair.split("=")
      if key in self.attrs:
        self.attrs[key] += ","+ value
      else:
        self.attrs[key] = value
  
  def is_feature(self):
    return len(self.fields) == 9


class Minter():
  """
  Class for minting IDs for gene and RNA features using the specified format.
  """
  def __init__(self, fp, idfmt = "%s%d"):
    self.oldids = {}
    self.scan_ids(fp)
    self.newids = {}
    self.oldnames = {}
    self.genecount = 0
    self.mint_new_ids(idfmt)

  def scan_ids(self, instream):
    """
    Scan the GFF3-formatted data from the given input stream and store the
    mapping of gene IDs to RNA IDs
    """
    for line in instream:
      entry = Entry(line)
      if not entry.is_feature() or entry.ftype not in ("mRNA", "tRNA", "rRNA"):
        continue
      
      assert "ID" in entry.attrs, "Error: RNA features must have ID attributes"
      rnaid = entry.attrs["ID"]
      assert "Parent" in entry.attrs, ("Error: RNA features must be connected "
                                       "to a gene feature via its Parent "
                                       "attribute")
      geneid = entry.attrs["Parent"]
      
      if geneid not in self.oldids:
        self.oldids[geneid] = []
      self.oldids[geneid].append(entry)

  def mint_new_ids(self, idfmt):
    """
    Create new IDs for each gene and RNA feature.
    """
    for geneid, rnalist in self.oldids.iteritems():
      self.genecount += 1
      rnacount = 0
      newgeneid = idfmt % ("GENE", self.genecount)
      self.newids[geneid] = newgeneid
  
      for rnaentry in rnalist:
        rnacount += 1
        rnaid = rnaentry.attrs["ID"]
        rnaname = rnaentry.attrs["Name"]
        rnatype = rnaentry.ftype.upper()
        rnafmt = idfmt+".%d"
        newrnaid = rnafmt % (rnatype, self.genecount, rnacount)
        self.newids[rnaid] = newrnaid
        self.oldnames[rnaid] = rnaname
  
  def fix_line(self, line, dbxref = None):
    """
    If the given line from a GFF3 file is a feature, replace any gene or RNA IDs
    located therein.
    """
    entry = Entry(line)
    if not entry.is_feature():
      return line

    replaceid = "ID" in entry.attrs and entry.attrs["ID"] in self.newids
    replaceparent = "Parent" in entry.attrs and entry.attrs["Parent"] in self.newids 
    if not replaceid and not replaceparent:
      return line

    if replaceid:
      newid = self.newids[entry.attrs["ID"]]
      line = re.sub("ID=[^;]+", "ID="+newid, line)
      if dbxref:
        line += ";Dbxref=%s:%s" % (dbxref, entry.attrs["ID"])

    if replaceparent:
      newparent = self.newids[entry.attrs["Parent"]]
      line = re.sub("Parent=[^;]+", "Parent="+newparent, line)

    return line

  def write_genemap(self, fp):
    if not fp:
      return
    for geneid in self.oldids.keys():
      print >> fp, "%s\t%s" % (self.newids[geneid], geneid)

  def write_rnamap(self, fp):
    if not fp:
      return
    for rnalist in self.oldids.values():
      for entry in rnalist:
        rnaid = entry.attrs["ID"]
        print >> fp, "%s\t%s\t%s" % (self.newids[rnaid], rnaid, self.oldnames[rnaid])


# Sequence parsing
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


# Assorted GFF3 processing functions
def strip_name(line):
  return re.sub(";*Name=[^;]+", "", line)

def strip_exon_id(line):
  if not re.search("\texon\t", line):
    return line
  return re.sub("ID=[^;]+;*", "", line)

def fix_cds_utr_id(line):
  labels = { "five_prime_UTR": "5putr", "three_prime_UTR": "3putr",
             "UTR": "utr", "CDS": "cds" }
  typematch = re.search("\t([^\t]*UTR|CDS)\t", line)
  if not typematch:
    return line

  ftype = typematch.group(1) 
  assert ftype in labels, "Error: unknown UTR type "+ ftype
  label = labels[typematch.group(1)]
  parentmatch = re.search("Parent=([^;]+)", line)
  assert parentmatch, "Error: CDS or UTR entry must have a Parent attribute"
  newid = "%s.%s" % (parentmatch.group(1), label)
  return re.sub("ID=[^;]+", "ID="+newid, line)
