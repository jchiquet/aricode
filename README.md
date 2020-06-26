
# aricode

<!-- badges: start -->
![R-CMD-check](https://github.com/jchiquet/aricode/workflows/R-CMD-check/badge.svg?branch=master)
[![CRAN
Status](https://www.r-pkg.org/badges/version/aricode)](https://CRAN.R-project.org/package=aricode)
[![Coverage
status](https://codecov.io/gh/jchiquet/aricode/branch/master/graph/badge.svg)](https://codecov.io/gh/jchiquet/aricode)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-blue.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![](https://img.shields.io/github/last-commit/jchiquet/aricode.svg)](https://github.com/jchiquet/aricode/commits/master)
<!-- badges: end -->

A package for efficient computations of standard clustering comparison
measures

## Installation

Stable version on the
[CRAN](https://cran.rstudio.com/web/packages/aricode/).

``` r
install.packages("aricode")
```

The development version is available via:

``` r
devtools::install_github("jchiquet/aricode")
```

## Description

Computation of measures for clustering comparison (ARI, AMI, NID and
even the \(\chi^2\) distance) are usually based on the contingency
table. Traditional implementations (e.g., function `adjustedRandIndex`
of package **mclust**) are in \(\Omega(n + u v)\) where

  - \(n\) is the size of the vectors the classifications of which are to
    be compared,
  - \(u\) and \(v\) are the respective number of classes in each
    vectors.

In **aricode** we propose an implementation, based on radix sort, that
is in \(\Theta(n)\) in time and space.  
Importantly, the complexity does not depends on \(u\) and \(v\). Our
implementation of the ARI for instance is one or two order of magnitude
faster than some standard implementation in `R`.

## Available measures and functions

The functions included in aricode are:

  - `ARI`: computes the adjusted rand index
  - `Chi2`: computes the Chi-square statistics
  - `MARI/MARIraw`: computes the modified adjusted rand index (Sundqvist
    et al, in preparation)
  - `NVI`: computes the the normalized variation information
  - `NID`: computes the normalized information distance
  - `NMI`: computes the normalized mutual information
  - `AMI`: computes the adjusted mutual information
  - `expected_MI`: computes the expected mutual information
  - `entropy`: computes the conditional and joint entropies
  - `clustComp`: computes all clustering comparison measures at once

## Timings

Here are some timings to compare the cost of computing the adjusted Rand
Index with **aricode** or with the commonly used function
`adjustedRandIndex` of the *mclust* package: the cost of the latter can
be prohibitive for large vectors:

![](man/figures/timings_plot-1.png)<!-- -->
