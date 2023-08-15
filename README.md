
<!-- README.md is generated from README.Rmd. Please edit that file -->



# SeffCovar

<!-- badges: start -->
<!-- badges: end -->

The goal of SeffCovar is to construct covariates for slide-effect-adjustment for epigenome-wide association studies.

## Installation

You can install SeffCovar using:

``` {.r}
library(devtools)
install_github("julianhecker/SeffCovar")
```

## Example: Assessing S\_high and S\_80\_100

``` {.r}
library(SeffCovar)
length(S_high)
#> [1] 1578
length(S_80_100)
#> [1] 1685
```

## Example: Construct slide effects covariates

This function returns 10 covariates. The number of covariates that need to be included in an association analysis depends and should evaluated in the respective study.

``` {.r}
methylation_matrix=matrix(rnorm(length(S_high)*100), nrow=length(S_high), ncol=100)
rownames(methylation_matrix)=S_high
slide_covars=get_slide_effects_covariates(methylation_matrix=methylation_matrix, input_set=S_high)
```
