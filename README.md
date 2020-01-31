## atlasqtl - Fast global-local hotspot QTL detection

[![Travis-CI Build Status](https://travis-ci.org/hruffieux/atlasqtl.svg?branch=master)](https://travis-ci.org/hruffieux/atlasqtl)

## Overview

**atlasqtl** is an R package implementing a scalable hierarchical modelling 
framework for variable selection in regression problems with many predictors and 
many responses. The method is tailored to the detection of hotspots, i.e., 
predictors associated with several responses. The *hotspot propensity* of each 
candidate predictor is modelled using a horseshoe distribution (Carvalho et al. 
2009), whose local scale flexibly models the large hotspot effects and whose 
global scale adapts to the overall signal sparsity. Inference is performed using 
efficient batch variational inference updates which are coupled with simulated 
annealing schemes to improve the exploration of multimodal parameter spaces. 

The method can be employed in any sparse multiple-response regression setting, 
and is particularly suited to large molecular quantitative trait locus (QTL) 
problems, in which hotspot genetic variants, controlling many molecular levels 
at once, may be responsible for decisive regulatory mechanisms. Hence, our 
approach is a tool that can help towards understanding the functional 
architecture underlying complex traits and outlining a *hotspot atlas* of the 
human genome. 

Reference: Helene Ruffieux, Anthony C. Davison, Jorg Hager, Jamie Inshaw, 
Benjamin P. Fairfax, Sylvia Richardson, Leonardo Bottolo, A global-local 
approach for detecting hotspots in multiple-response regression, 
arXiv:1811.03334, 2018.

## Installation

**Important note:** the R package depends on the GSL library which needs to be manually 
installed. For example on Ubuntu,

``` bash
$ sudo apt-get install libgsl-dev
```

or on mac,

``` bash
$ brew install gsl
```

after having installed Homebrew.

Then, to install the package in R, run the following command:

``` r
# after having installed devtools (install.packages("devtools"))
devtools::install_github("hruffieux/atlasqtl")
```

## License and authors

This software uses the GPL v3 license, see [LICENSE](LICENSE).
Authors and copyright are provided in [DESCRIPTION](DESCRIPTION).

## Issues

To report an issue, please use the [atlasqtl issue tracker](https://github.com/hruffieux/atlasqtl/issues) at github.com.
