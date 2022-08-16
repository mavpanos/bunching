
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bunching

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/bunching)](https://cran.r-project.org/package=bunching)
[![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/bunching)](https://cranlogs.r-pkg.org/badges/bunching)

<img src='man/figures/bunching_logo.png' style='width: 20%; height: 20%' />

## Overview

`bunching` allows the user to conduct bunching analysis in a kink or
notch setting and returns a rich set of results. Important features of
the package include functionality to control for (different levels of)
round-number bunching or other bunching masses within the estimation
bandwidth, options to split bins by placing the bunching point as the
minimum, median or maximum in its bin (for robustness analysis), and
produces estimates of both parametric and reduced-form versions of
elasticities associated with the bunching mass. It also returns
(residual-based) bootstrapped vectors of estimates of all the main
estimable parameters (bunching mass, elasticity, marginal buncher,
dominated region, fraction in dominated region in notch setting, etc.)
which can be used for further analysis, e.g. incorporation into
structural models that rely on bunching moments. The package also
provides an exploratory viusalization function to speed up pre-analysis,
and produces plots in the [Chetty et
al. (2011)](https://doi.org/10.1093/qje/qjr013) style with lots of
options on editing the plot appearance.

## Installation

You can install the released version of bunching from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("bunching")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mavpanos/bunching")
```

If you find a bug, please file a minimal reproducible example in the
[issues](https://github.com/mavpanos/bunching/issues).
