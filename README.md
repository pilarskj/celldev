# Assessing the inference of single-cell phylogenies and population dynamics from CRISPR lineage recordings

Julia Pilarski<sup>1,2</sup>, Tanja Stadler<sup>1,2</sup>, Sophie Seidel<sup>1,2</sup>

<sup>1</sup>Department of Biosystems Science and Engineering, ETH Zürich, Basel, Switzerland 
<sup>2</sup>Swiss Institute of Bioinformatics (SIB), Lausanne, Switzerland

PLOS Computational Biology, 2026

This repository contains code to reproduce the analyses and figures associated with the [paper](<https://doi.org/10.1371/journal.pcbi.1014370>).

------------------------------------------------------------------------

We simulated CRISPR lineage recordings and applied the Bayesian framework in [BEAST 2](https://www.beast2.org) for the joint inference of time-scaled cell lineage trees and parameters of the cell population processes, i.e., cell division, death, and differentiation. 
We evaluated the inference performance under various conditions and assessed how much information on cell development the recordings provide.

We used the following BEAST 2 packages:

[TiDeTree](https://github.com/seidels/tidetree): models CRISPR-Cas9 editing on multiple independent target sites; editing outcomes are random indels.

[SciPhy](https://github.com/azwaans/SciPhy): models sequential CRISPR-prime editing on tandem arrays of target sites (tapes); editing outcomes are short template-based insertions.
*Note:* SciPhy was initially developed for the analysis of DNA Typewriter recordings (see [Choi et al., 2022](https://www.nature.com/articles/s41586-022-04922-8)), hence, 'typewriter' is referred to in many file names and functions.

Further, we used the birth-death models from [BDSKY](https://github.com/BEAST2-Dev/bdsky/tree/master/src/bdsky/evolution/speciation) and [BDMM-Prime](https://github.com/tgvaughan/BDMM-Prime/tree/master) for phylodynamic inference.

------------------------------------------------------------------------

The repository adheres to the following structure: 
We considered the development of **homogeneous** (single-type) and **heterogeneous** (multi-type) cell populations. 
For both, our workflow consisted of four steps:

1.  Simulation of lineage trees

2.  Simulation of CRISPR lineage recordings (resulting in 'barcodes')

3.  Bayesian inference of cell phylogenies and phylodynamic parameters

4.  Analysis and evaluation of the inference performance

Additionally, the directory 'figures' contains code to reproduce the graphs from the manuscript.

------------------------------------------------------------------------

System requirements: We run the inference on a HPC cluster with the Slurm batch system and a Java module provided. We automatized the simulations and inference using Snakemake. We analyzed the data in R.

BEAST2 package versions: BEAST v2.7.7, TiDeTree v1.0.0, SciPhy v0.0.1, BDSKY v1.5.1, BDMM-Prime v1.0.0, feast v10.4.0.
