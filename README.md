
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bunching <img src='man/figures/bunching_logo.png' align="right" height="160" />

<!-- badges: start -->

<!-- badges: end -->

`bunching` allows the user to conduct bunching analysis in a kink or
notch setting and returns a rich set of results. Important features of
the package include functionality to control for (different levels of)
round-number bunching or other bunching masses within the estimation
bandwidth, options to split bins by placing the bunching point as the
minimum, median or maximum in its bin (for robustness analysis), and
estimates of both parametric and reduced-form versions of elasticities
associated with the bunching mass. It also provides an exploratory
viusalization function to speed up pre-analysis, and produces plots in
the Chetty et al. (2011) style with lots of options on editing the plot
appearance. Further, it returns bootstrapped vectors of estimates of all
the main estimable parameters (bunching mass, elasticity, marginal
buncher, dominated region, fraction in dominated region in notch
setting, etc.) which can be used for further analysis,
e.g. incorporation into structural models that rely on bunching
moments.

## Installation

<!--You can install the released version of bunching from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("bunching")
```
-->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mavpanos/bunching")
```

## Some Examples of bunching at a kink

This is a basic example which shows you how to solve a common problem:

``` r
library(bunching)
#> Loading required package: BB
#> Loading required package: dplyr
#> Warning: package 'dplyr' was built under R version 3.5.2
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> Loading required package: ggplot2
#> Warning: package 'ggplot2' was built under R version 3.5.2
#> Loading required package: tidyr
data("bunching_data")
kink1 <- bunchit(z_vector = bunching_data$kink, zstar = 10000, binwidth = 50,
                 bins_l = 20, bins_r = 20, poly = 4, t0 = 0, t1 = .2, p_b = TRUE)
kink1$plot
```

<img src="man/figures/README-example-1.png" width="50%" style="display: block; margin: auto;" />
