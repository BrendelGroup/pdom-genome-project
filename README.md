# *Polistes dominula* genome project lab notebook

Daniel Standage  
Volker Brendel  
Amy Toth, PI

## Overview

This repository is the electronic lab notebook for the *Polistes dominula* genome project.
The intellectual thrust of the work and its various findings are discussed in *[insert citation here]*.
And although our methods are described in detail in the supplementary **Materials and Methods** section, there are many technical details required to replicate the work (such as configuration files, UNIX commands, and source code) that are not easily communicated in the format of a scientific paper.
This notebook is a complement to the paper and is intended to provide those technical details.

## Disclaimer

Ideally, the code and commands provided in this lab notebook should simply work out-of-the-box on any UNIX system.
Practically, though, the vast majority of the analysis was done on machines running the Fedora and Mac OS X operating systems.
Getting the code to run on other systems may require minor adaptations.

Also, although we have made every effort to ensure the lab notebook content is accurate, it's possible errors remain.
Some of the content was created *post hoc* from personal notes, so there may be minor typos, copy/paste errors, or occasional omitted dependency or technical prerequisite statements (such as *make sure program X is in your path*).
If you run into any issues reading or using the documentation or supplementary code, we'd be happy to help if you open up a ticket in the [lab notebook's issue tracker](https://github.com/BrendelGroup/pdom-genome-project/issues).

## Organization

For each analysis (or group of related analyses) described in the paper, this repository contains a directory with the following.

  - A `README.md` markdown file contains a prose description of te analysis, along with the commands someone would type into the terminal to reproduce the analysis.
  - A `Makefile` or batch script that enables someone to reproduce the entire analysis with a single command.
  - Any relevant scripts or configuration files.

The README files can be concatenated and rendered as a single PDF file with [Pandoc](http://johnmacfarlane.net/pandoc/).
They can also be viewed individually on the web at [GitHub](https://github.com/BrendelGroup/pdom-genome-project).
Anyone wanting access to the scripts and Makefiles can simply (and anonymously) clone the repository with git or download an archive file containing all of the repository's contents.