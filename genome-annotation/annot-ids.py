#!/usr/bin/env python
import getopt, re, sys
from idtidy import *

def print_usage(outstream):
  usage = ("\n"
    "Create new IDs for gene and RNA features using the specified format, tested\n"
    "primarily on output from the Maker genome annotation pipeline.\n\n"
    "Usage: python annot-ids.py [options] annot.gff3\n"
    "  Options:\n"
    "    -f|--idfmt: STRING     printf-style format to use for creating new IDs;\n"
    "                           must accept a string (%s) and long unsigned\n"
    "                           integer (%d or %lu); default is '%s%d'\n"
    "    -g|--genemap: FILE     write the correspondence of new gene IDs to old\n"
    "                           gene IDs to the given file\n"
    "    -h|--help              print this help message and exit\n"
    "    -n|--stripnames        remove names from all features; particularly useful\n"
    "                           when name attributes are uninformative, as is the\n"
    "                           case with Maker\n"
    "    -o|--outfile: FILE     file to which output will be written; default is\n"
    "                           the terminal (standard output)\n"
    "    -r|--rnamap: FILE      write the correspondence of new RNA IDs to old\n"
    "                           RNA IDs to the given file\n"
    "    -x|--dbxref: STRING    use the 'Dbxref' attribute (database cross\n"
    "                           reference) to store a copy of any ID that is\n"
    "                           replaced; the argument provided will serve as the\n"
    "                           Dbxref key\n")
  print >> outstream, usage

def parse_options(argv):
  params = {
    "dbxref": None,
    "genemap": None,
    "idfmt": "%s%d",
    "outfile": sys.stdout,
    "rnamap": None,
    "stripnames": False,
    "infile": None
  }
  optstr = "f:g:hno:r:x:"
  longopts = ["idfmt=", "genemap=", "help", "stripnames", "outfile=", "rnamap=",
              "dbxref="]
  (options, args) = getopt.getopt(argv[1:], optstr, longopts)
  for key, value in options:
    if key in ("-f", "--idfmt"):
      params["idfmt"] = value
    elif key in ("-g", "--genemap"):
      params["genemap"] = open(value, "w")
    elif key in ("-h", "--help"):
      print_usage(sys.stdout)
      sys.exit(0)
    elif key in ("-n", "--stripnames"):
      params["stripnames"] = True
    elif key in ("-o", "--outfile"):
      params["outfile"] = open(value, "w")
    elif key in ("-r", "--rnamap"):
      params["rnamap"] = open(value, "w")
    elif key in ("-x", "--dbxref"):
      params["dbxref"] = value
    else:
      print_usage(sys.stderr)
      assert False, "unsupported option '%s'" % key

  if len(args) != 1:
    print_usage(sys.stderr)
    assert False, "expected 1 argument, %d provided" % len(args)
  if args[0] == "-":
    params["infile"] = sys.stdin
  else:
    params["infile"] = open(args[0])

  return params

if __name__ == "__main__":
  params = parse_options(sys.argv)

  # Load data into memory; increases runtime and memory consumption, but
  # required for correct handling of stdin, named pipes, or process substitutions
  indata = []
  for line in params["infile"]:
    indata.append(line)

  minter = Minter(indata, params["idfmt"])
  for line in indata:
    line = line.rstrip()
    line = minter.fix_line(line, params["dbxref"])
    if params["stripnames"]:
      line = strip_name(line)
    line = strip_exon_id(line)
    line = fix_cds_utr_id(line)
    print >> params["outfile"], line

  minter.write_genemap(params["genemap"])
  minter.write_rnamap(params["rnamap"])
