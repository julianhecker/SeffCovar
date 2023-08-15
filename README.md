
<!-- README.md is generated from README.Rmd. Please edit that file -->



# SeffCovar

<!-- badges: start -->
<!-- badges: end -->

The goal of SeffCovar is to construct covariates for slide-effect-adjustment for epigenome-wide association studies. **!Warning! This is a test version and not functional yet. This package will be available soon.**

## Installation

You can install SeffCovar using:

``` {.r}
library(devtools)
install_github("julianhecker/SeffCovar")
```

## Example

``` {.r}
library(SeffCovar)
# assess S_80_100 and S_high
length(S_high)
#> [1] 1578
length(S_80_100)
#> [1] 1685

# construct slide effects covariates
methylation_matrix=matrix(rnorm(length(S_high)*100), nrow=length(S_high), ncol=100)
rownames(methylation_matrix)=S_high
slide_covars=get_slide_effects_covariates(methylation_matrix=methylation_matrix, input_set=S_high)
```
