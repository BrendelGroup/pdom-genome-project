#!/usr/bin/make -f
SHELL=bash

pdom-docs.pdf:	
		which pandoc
		test -h pdom-size-kmers.png || ln -s genome-size/pdom-size-kmers.png
		sed -e 's|<sub>|~|g' -e 's|</sub>|~|g' README.md genome-size/README.md \
		    | pandoc --number-sections --standalone -f markdown --no-wrap \
                      -V geometry:margin=1.15in --highlight-style zenburn -o $@

clean:		
		rm -f pdom-docs.pdf pdom-size-kmers.png
