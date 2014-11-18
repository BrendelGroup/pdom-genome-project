#!/usr/bin/make -f
SHELL=bash

READMES=genome-size/README.md genome-assembly/README.md transcript-assembly/README.md repeat-masking/README.md \
	transcript-alignment/README.md protein-alignment/README.md genome-annotation/README.md \
	interval-loci/README.md novel-genes/README.md diff-exp/README.md splicing/README.md

pdom-docs.pdf:	README.md $(READMES)
		which pandoc
		sed -e 's|<sub>|~|g' -e 's|</sub>|~|g' $^ \
		    | grep -v '^\!\[' \
		    | pandoc --number-sections --standalone -f markdown --no-wrap \
		      -V geometry:margin=1in --highlight-style zenburn -o $@

