## atlasqtl – fAsT gLobal-locAl hotSpot QTL detection

## Overview

**atlasqtl** is an R package implementing a flexible hierarchical modelling of 
hotspot predictors in multiple-response regression, i.e., predictors which are
associated with several responses simultaneously. The *hotspot propensity* of 
each candidate predictor is modelled using a horseshoe distribution (Carvalho 
et al. 2009), whose local scales flexibly capture the large hotspot effects 
and whose global scale adapts to the overall sparsity level. Inference is 
performed using highly-scalable variational inference updates which are coupled 
with simulated annealing schemes to improve the exploration of multimodal 
parameter spaces. 

The method can be employed in any sparse multiple-response regression settings, 
and is particularly suited to large molecular quantitative trait locus (QTL) 
problems, in which hotspot genetic variants, controlling many molecular levels 
at once, may be responsible for decisive regulatory mechanisms and shape the 
functional architecture underlying complex traits; our approach is a tool to 
better understand this architecture and will hopefully serve towards outlining a 
*hotspot atlas* of the human genome. 

Reference: Hele Ruffieux, Anthony C. Davison, JöroHager, Jamie Inshaw, 
Benjamin P. Fairfax, Sylvia Richardson, Leonardo Bottolo, A global-local 
approach for detecting hotspots in multiple-response regression, 
arXiv:1811.03334, 2018.

## Warning

**This is a development branch**, it is not guaranteed to be stable at any given time
and features are subject to change. Please use the [stable version](https://github.com/hruffieux/atlasqtl),
unless you want to test and report issues.

## Installation

**Important note:** the R package depends on `libgsl` which needs to be manually 
installed. For example on Ubuntu,

``` bash
$ sudo apt-get install libgsl-dev
```

To install the package in R, run the following command:

``` r
# after having installed devtools (install.packages("devtools"))
devtools::install_github("hruffieux/atlasqtl", ref = "devel")
```

## License and authors

This software uses the GPL v3 license, see [LICENSE](LICENSE).
Authors and copyright are provided in [DESCRIPTION](DESCRIPTION).

## Issues

To report an issue, please use the [atlasqtl issue tracker](https://github.com/hruffieux/atlasqtl/issues) at github.com.
