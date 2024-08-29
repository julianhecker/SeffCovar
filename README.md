
<!-- README.md is generated from README.Rmd. Please edit that file -->



# SeffCovar

<!-- badges: start -->
<!-- badges: end -->

The SeffCovar R package constructs covariates for slide-effect-adjustment in epigenome-wide association studies or integrative methylation analyses.

## Installation

You can install SeffCovar using:

``` {.r}
library(devtools)
install_github("julianhecker/SeffCovar")
```

## Sets S\_0\_20, S\_20\_40, S\_40\_60, S\_60\_80, S\_80\_100, and S\_high, as described in Hecker et al. 2023 [1].

``` {.r}
library(SeffCovar)
length(S_high)
#> [1] 1578
head(S_80_100)
#> [1] "cg06352483" "cg17517854" "cg25060657" "cg16441688" "cg06790069"
#> [6] "cg24877558"
```

## Approach 1: get\_slide\_effects\_covariates

This function extracts the DNA methylation values for CpGs in S\_high and performs a singular value decomposition, after scaling and centering. The resulting vectors are the PCs that can be used as covariates for slide effect adjustment.

``` {.r}
methylation_matrix=matrix(rnorm(length(S_high)*100), nrow=length(S_high), ncol=100)
rownames(methylation_matrix)=S_high
slide_covars=get_slide_effects_covariates(methylation_matrix=methylation_matrix, input_set=S_high)
head(slide_covars)
#>              sPC1         sPC2        sPC3        sPC4         sPC5        sPC6
#> [1,] -0.007415158  0.002678847 -0.05190384  0.01950551 -0.052571020  0.08927897
#> [2,]  0.126140851  0.011427023  0.02129938 -0.04320385  0.098750835  0.06591822
#> [3,] -0.121005282  0.053914684 -0.11519284 -0.09855958  0.003000299 -0.11682964
#> [4,]  0.001750477 -0.051667027 -0.05386999 -0.01307616 -0.204258733 -0.00522347
#> [5,]  0.091967829 -0.057450799  0.01694561  0.05728267 -0.050705562 -0.05139665
#> [6,]  0.017375320  0.004592727  0.11386451  0.01161123  0.099079753 -0.19039631
#>             sPC7        sPC8        sPC9        sPC10
#> [1,] -0.12634653  0.06724397 -0.04795647  0.009386547
#> [2,] -0.03374496  0.01367125 -0.18363086 -0.001124245
#> [3,]  0.08911991 -0.01656946  0.06720038 -0.026228370
#> [4,]  0.01373036 -0.05334789 -0.06819347 -0.045742609
#> [5,] -0.03431053 -0.22544090 -0.12627388  0.212177918
#> [6,] -0.12389500  0.03492531 -0.02910057 -0.092999069
```

## Approach 2: get\_ComBat\_slide\_effect\_covariates

This function extracts the DNA methylation values for CpGs in S\_high and performs part of the ComBat adjustment [2]. The function uses the estimated slide effect gamma's and performs a singular value decomposition on the resulting matrix. Again, this leads to PCs that can be used as covariates for slide effect adjustment. The functions is based on a modification of the ComBat function in the sva package [3]. The original code for ComBat from the sva package that can be found at <https://bioconductor.org/packages/release/bioc/html/sva.html>.

``` {.r}
methylation_matrix=matrix(rnorm(length(S_high)*100), nrow=length(S_high), ncol=100)
rownames(methylation_matrix)=S_high
slide=rep(paste0("slide",1:10),10)
slide_covars=get_ComBat_slide_effect_covariates(methylation_matrix=methylation_matrix, slide=slide, mod=NULL, input_set=S_high)
#> Found 10 slides
#> Adjusting for 0 covariate(s) or covariate level(s)
#> Fitting L/S model and finding priors
head(slide_covars)
#>             sPC1        sPC2        sPC3        sPC4        sPC5        sPC6
#> [1,] -0.05079121 -0.10007859  0.11796097  0.05423742  0.10134942  0.12336540
#> [2,]  0.09167848  0.03005084  0.13499129  0.03155485 -0.19938378 -0.05358116
#> [3,]  0.11295137  0.06668619 -0.08854503  0.03278284  0.09686079 -0.15014488
#> [4,] -0.08000584  0.03800528  0.02684123 -0.11874507  0.05332456  0.07682950
#> [5,]  0.08958204 -0.19565616 -0.12935874 -0.08339406 -0.10283713  0.06512774
#> [6,]  0.15183267  0.00307854  0.01860855  0.09597305  0.12879061  0.05293499
#>             sPC7        sPC8         sPC9       sPC10
#> [1,]  0.11383237  0.13022664  0.071697485 -0.09001980
#> [2,] -0.07544029  0.11447111  0.007344753 -0.02050830
#> [3,]  0.12799230  0.08848327 -0.087995395 -0.15647280
#> [4,] -0.10458864  0.04148537 -0.214307538 -0.26000654
#> [5,]  0.04659651 -0.04847285 -0.025836514  0.11356391
#> [6,] -0.15297080 -0.11267143  0.043545431 -0.01472454
```

# References

[1] A consistent pattern of slide effects in Illumina DNA methylation BeadChip array data. Hecker J, Lee S, Kachroo P, Prokopenko D, Maaser-Hecker A, Lutz SM, Hahn G, Irizarry R, Weiss ST, DeMeo DL, Lange C. Epigenetics. 2023 Dec;18(1):2257437. doi: 10.1080/15592294.2023.2257437. Epub 2023 Sep 20. PMID: 37731367

[2] Adjusting batch effects in microarray expression data using empirical Bayes methods. Johnson WE, Li C, Rabinovic A. Biostatistics. 2007 Jan;8(1):118-27. doi: 10.1093/biostatistics/kxj037. Epub 2006 Apr 21. PMID: 16632515

[3] The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD. Bioinformatics. 2012 Mar 15;28(6):882-3. doi: 10.1093/bioinformatics/bts034. Epub 2012 Jan 17. PMID: 22257669
