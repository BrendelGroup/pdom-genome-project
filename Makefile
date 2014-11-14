#!/usr/bin/make -f
SHELL=bash

pdom-docs.pdf:	README.md genome-size/README.md pdom-size-kmers.png genome-assembly/README.md
		which pandoc
		sed -e 's|<sub>|~|g' -e 's|</sub>|~|g' README.md genome-size/README.md genome-assembly/README.md \
		    | pandoc --number-sections --standalone -f markdown --no-wrap \
		      -V geometry:margin=1in --highlight-style zenburn -o $@

pdom-size-kmers.png:	genome-size/pdom-size-kmers.png
			ln -s genome-size/pdom-size-kmers.png

clean:		
		rm -f pdom-size-kmers.png
