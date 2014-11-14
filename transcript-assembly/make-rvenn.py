#!/usr/bin/env python
import subprocess, sys

def line_count(filename):
  """
  Stolen shamelessly from http://stackoverflow.com/a/845069/459780
  """
  sp = subprocess.Popen(['wc', '-l', filename],
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
  result, err = sp.communicate()
  if sp.returncode != 0:
    raise IOError(err)
  return int(result.strip().split()[0])

def file_line_counts():
  counts = {}

  counts["Pc"] = line_count("vID-Pc")
  counts["Pd"] = line_count("vID-Pd")
  counts["Pm"] = line_count("vID-Pm")

  counts["PcPd"] = min(line_count("vID-PcPd"),
                       line_count("vID-PdPc"))
  counts["PcPm"] = min(line_count("vID-PcPm"),
                       line_count("vID-PmPc"))
  counts["PdPm"] = min(line_count("vID-PdPm"),
                       line_count("vID-PmPd"))

  counts["PcPdPm"] = min(line_count("vID-PcPdPm"),
                         line_count("vID-PdPcPm"),
                         line_count("vID-PmPcPd"))

  return counts

def venn_ranges(counts):
  r = {}
  i = 0

  r["Pd"] = "%d:%d" % (i+1, i+counts["Pd"])
  i += counts["Pd"]

  r["PcPd"] = "%d:%d" % (i+1, i+counts["PcPd"])
  i += counts["PcPd"]

  r["PcPdPm"] = "%d:%d" % (i+1, i+counts["PcPdPm"])
  i += counts["PcPdPm"]

  r["PdPm"] = "%d:%d" % (i+1, i+counts["PdPm"])
  i += counts["PdPm"]

  r["Pc"] = "%d:%d" % (i+1, i+counts["Pc"])
  i += counts["Pc"]

  r["PcPm"] = "%d:%d" % (i+1, i+counts["PcPm"])
  i += counts["PcPm"]

  r["Pm"] = "%d:%d" % (i+1, i+counts["Pm"])
  i += counts["Pm"]

  return r

if __name__ == "__main__":
  counts = file_line_counts()
  ranges = venn_ranges(counts)

  print ('library("VennDiagram")\n'
'venn.diagram(\n'
'    x = list(\n'
'            Pd = c('+ ranges["Pd"] +', '+ ranges["PcPd"] +', '+ ranges["PcPdPm"] +', '+ ranges["PdPm"] +'),\n'
'            Pc = c('+ ranges["Pc"] +', '+ ranges["PcPd"] +', '+ ranges["PcPdPm"] +', '+ ranges["PcPm"] +'),\n'
'            Pm = c('+ ranges["Pm"] +', '+ ranges["PcPdPm"] +', '+ ranges["PdPm"] +', '+ ranges["PcPm"] +')\n'
'            ),\n'
'    filename = "pdom-venn.png",\n'
'    imagetype = "png",\n'
'    col = "transparent",\n'
'    fill = c("red", "blue", "green"),\n'
'    alpha = 0.5,\n'
'    label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"),\n'
'    cex = 2.5,\n'
'    fontfamily = "serif",\n'
'    fontface = "bold",\n'
'    cat.default.pos = "text",\n'
'    cat.col = c("darkred", "darkblue", "darkgreen"),\n'
'    cat.cex = 2.5,\n'
'    cat.fontfamily = "serif",\n'
'    cat.dist = c(0.06, 0.06, 0.03),\n'
'    cat.pos = 0\n'
');')
