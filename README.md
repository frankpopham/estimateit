
<!-- README.md is generated from README.Rmd. Please edit that file -->

# estimateit

The goal of `estimateit` is to estimate marginal and conditional (for
the ATE only at present) effects using inverse probability weights
generated by [`weightit`](https://github.com/ngreifer/WeightIt). Effects
are generated using `svyglm` and `svycontrast` from the`survey package.`
They include the effect for the treated, the effect for the untreated or
control, their difference, their relative difference, and their odds
ratio. As the previous sentence implies `estimateit` only works with
binary exposures and binary outcomes at present.

## Installation

You can install the development version of `estimateit` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("frankpopham/estimateit")
```

## Examples

``` r
## A simple example
library(WeightIt)
library(tibble)
library(tidyr)
dfvi <- tibble(
 C = rep(0:1, each = 4),
 X = rep(0:1, times = 4),
 Y = rep(0:1, times = 2, each = 2),
 N = c(96, 36, 64, 54, 120, 120, 30, 480)
 ) %>%
 uncount(N)
 W1 <- weightit(X ~ C, data = dfvi,
               method = "ps", estimand = "ATE")
summary(W1)
E1 <- estimateit(weightitobj=W1, outcome=Y, data=dfvi)
E1
```
