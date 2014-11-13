#!/usr/bin/make -f
SHELL=bash

pdom-docs.pdf:	README.md genome-size/README.md
		which pandoc
		test -h pdom-size-kmers.png || ln -s genome-size/pdom-size-kmers.png
		sed -e 's|<sub>|~|g' -e 's|</sub>|~|g' $^ \
		    | pandoc --number-sections --standalone -f markdown --no-wrap \
		      -V geometry:margin=1in --highlight-style zenburn -o $@

clean:		
		rm -f pdom-size-kmers.png
