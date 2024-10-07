# celldev

This repository contains code to reproduce the results of the paper **"Assessing the inference of single-cell phylogenies and population dynamics from genetic lineage tracing data"**.

------------------------------------------------------------------------

We simulate CRISPR-based lineage recordings and apply Bayesian inference in BEAST2 for the joint inference of time-scaled lineage trees and parameters of the cell population processes. We evaluate the inference performance under various conditions and assess how much information on cell development the recordings provide.

We use the following BEAST2 packages:

[TiDeTree](https://github.com/seidels/tidetree): models CRISPR-Cas9 editing on multiple independent target sites; editing outcomes are indels

[SciPhy](https://github.com/azwaans/SciPhy/tree/master): models sequential CRISPR-prime editing on tandem arrays of target sites (tapes); editing outcomes are short template-based insertions

*Note:* SciPhy was initially developed for the analysis of DNA Typewriter recordings (see [Choi et al., 2022](https://www.nature.com/articles/s41586-022-04922-8)). The code in this repository was written with a previous version of the SciPhy package which used the naming 'typewriter' instead of 'sciphy' (see commit [#30 Renaming classes](https://github.com/azwaans/SciPhy/commit/ead14aa57874a6c8157cba155f288ad8bf28707e)).

Further, we use the phylodynamic models [BDSKY](https://github.com/BEAST2-Dev/bdsky/tree/master/src/bdsky/evolution/speciation) and [BDMM-Prime](https://github.com/tgvaughan/BDMM-Prime/tree/master){.uri}, and the add-on [feast](https://github.com/tgvaughan/feast).

------------------------------------------------------------------------

In our study, we consider the development of **homogeneous** (single-type) and **heterogeneous** (multi-type) cell populations. For both, the workflow consists of the following steps:

1.  Simulation of lineage trees

2.  Simulation of CRISPR-based lineage recordings (resulting in 'barcodes')

3.  Bayesian inference of cell phylogenies and phylodynamic parameters

Analysis and evaluation of the inference performance

------------------------------------------------------------------------

*TODO: add installation tips (java, .jar file - in bin?, snakemake, HPC with Slurm?, R)*
