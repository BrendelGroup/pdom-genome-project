Polistes dominula* genome project

- [Daniel Standage][]
- [Volker Brendel][]
- [Amy Toth][], Principal Investigator

## Overview

This documentation is a record of our work for the *Polistes dominula* genome project.
It was created to 1) serve as full disclosure of all of the methods, commands, and software used to produce the reported results, and 2) facilitate anonymous replication of those results.

### Data access

Raw instrument data and final data outputs are stored in the [iPlant Data Store][] under the path `/iplant/home/standage/Polistes_dominula/`.
All file and directory paths provided in this documentation are relative to that root path, which for the remainder of the documentation will be designated the **Pdom Data Store**.

### Using this documentation

This project is divided into several sections, with each section focusing on a single analysis or small group of related analyses.
Each section has a dedicated directory containing code and documentation specific to that section.
These resources can be browsed or downloaded at [GitHub][].

  - A `README.md` file (in Markdown format) is included for each section, which provides a prose description of what each set of commands is doing.
    This file is intended to facilitate interactive replication of results: typing or pasting the commands into the terminal and executing them manually to produce the output.
    (Note: a single PDF document containing all documentation was produced by concatenating all of the various README files into a single Markdown file and converting to PDF format.)
  - Each section also includes a `run.sh` file (Bash code) which includes the same commands as the README file--sans commentary--and is intended to facilitate replication in batch mode.
  - Most sections also include additional supplementary files, such as source code, graphics, or configuration files necessary for replicating the results. The purpose of each supplemental file should be clear from the documentation.

If you encounter any problems using this documentation or its associated files, please open a ticket with the [Pdom Genome Project issue tracker][].

[Daniel Standage]: http://standage.github.io/
[Volker Brendel]: http://brendelgroup.org/
[Amy Toth]: http://www.public.iastate.edu/~amytoth/Toth_lab/Home.html
[iPlant Data Store]: http://www.iplantcollaborative.org/ci/data-store
[GitHub]: https://github.com/BrendelGroup/pdom-genome-project
[Pdom Genome Project issue tracker]: https://github.com/BrendelGroup/pdom-genome-project/issues

