# Agent-based modelling reveals the role of the tumour microenvironment on the short-term success of combination temozolomide/immune checkpoint blockade to treat glioblastoma
Codes associated with the paper titled "Agent-based modelling reveals the role of the tumour microenvironment on the short-term success of combination temozolomide/immune checkpoint blockade to treat glioblastoma". The code is developed in the C++ based PhysiCell framework, an Open Source Physics-Based Cell Simulator for 3-D Multicellular systems.

Reference: A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellular Systems, PLoS Comput. Biol. 14(2): e1005991, 2018. DOI: 10.1371/journal.pcbi.1005991.

Visit http://physicell.org for the latest tutorials and help.


Makefile rules to compile and run the code:

make GBM-immune-TMZ: populates the PhysiCell environment with the GBM-immune-TMZ project. Use "make" to compile it. After the compilation, type .\project.exe to run the code.

make: compiles the current project. If no project has been defined, it first populates the heterogeneity sample project.

make clean: removes all .o files and the executable, so that the next "make" recompiles the entire project

make data-cleanup: clears all simulation data

make reset: de-populates the sample project and returns to the original PhysiCell state. Use this when switching to a new PhysiCell sample project.
