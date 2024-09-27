# Assessing the inference of single-cell phylogenies and population dynamics from genetic lineage tracing data

In this repository, we simulate CRISPR-based lineage recordings and apply Bayesian inference in BEAST2 for the joint inference of time-scaled lineage trees and parameters of the cell population processes. We evaluate the inference performance under various conditions and assess how much information on cell development the recordings provide.

We use the following BEAST2 packages:

[TiDeTree](https://github.com/seidels/tidetree): models CRISPR-Cas9 editing on multiple independent target sites; editing outcomes are indels

[SciPhy](https://github.com/azwaans/SciPhy/tree/master): models sequential CRISPR-prime editing on tandem arrays of target sites (tapes); editing outcomes are short template-based insertions

Further, we use the phylodynamic models [BDSKY](https://github.com/BEAST2-Dev/bdsky/tree/master/src/bdsky/evolution/speciation) and [BDMM-Prime](https://github.com/tgvaughan/BDMM-Prime/tree/master){.uri}.

*Note:* SciPhy was initially developed for the analysis of DNA Typewriter recordings (see [Choi et al., 2022](https://www.nature.com/articles/s41586-022-04922-8)). The code in this repository was written with a previous version of the SciPhy package which used the naming 'typewriter' instead of 'sciphy' (see commit [#30 Renaming classes](https://github.com/azwaans/SciPhy/commit/ead14aa57874a6c8157cba155f288ad8bf28707e)). **TODO: mention which files should be adapted.**

**TODO: describe how to replicate the study(?)**
