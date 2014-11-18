# *Polistes dominula* genome project

Daniel Standage  
Volker Brendel  
Amy Toth, PI

As contributors to this project we are committed to open and reproducible research.
The main challenge for documenting our work is how to share the documentation and supplementary code in a way that is accessble both to readers with a conceptual focus (particularly reviewers) as well as readers with a technical focus (those interested in reproducing our work or adapting it for their system).
This repository is our answer to that challenge.

For each analysis (or group of related analyses) described in the paper, this repository contains a directory with the following.

  - A `README.md` markdown file contains a prose description of te analysis, along with the commands someone would type into the terminal to reproduce the analysis.
  - A `Makefile` or batch script that enables someone to reproduce the entire analysis with a single command.
  - Any relevant scripts or configuration files.

The README files can be concatenated and rendered as a single PDF file with [Pandoc](http://johnmacfarlane.net/pandoc/).
They can also be viewed individually on the web at [GitHub](https://github.com/BrendelGroup/pdom-genome-project).
Anyone wanting access to the scripts and Makefiles can simply (and anonymously) clone the repository with git or download an archive file containing all of the repository's contents.

---------

**Disclaimer**: Although we have made every effort to ensure the docs are correct, it's possible errors remain.
Some of the Makefiles were created *post hoc* from lab notebook entries, so there may be typos, or copy/paste errors; sometimes a dependency or technical prerequisite (i.e. *make sure program X is in your path*) may be accidentally ommitted.
That being said, we don't expect any of these issues to be substantial.
If you run into any issues reading or using the documentation or supplementary code, please feel free to open up a ticket in the [Pdom Genome Project issue tracker](https://github.com/BrendelGroup/pdom-genome-project/issues).
