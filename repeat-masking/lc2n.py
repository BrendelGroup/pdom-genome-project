#!/usr/bin/env python
import re, string, sys

for line in sys.stdin:
  line = line.rstrip()
  if not line.startswith(">"):
    line = re.sub('[a-z]', 'N', line)
  print line
