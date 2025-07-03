
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

## Sets S\_0\_20, S\_20\_40, S\_40\_60, S\_60\_80, S\_80\_100, and S\_high, as described in Hecker et al. 2023 [1,4].

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
#>             sPC1        sPC2        sPC3         sPC4        sPC5        sPC6
#> [1,] -0.03074495  0.02835304 -0.02303973  0.013008619 -0.02853287  0.01652519
#> [2,] -0.10045291 -0.04211759 -0.04996575  0.009284407  0.10414371 -0.01921734
#> [3,]  0.13452460  0.03581853  0.18358021  0.032288163 -0.14632155  0.03930918
#> [4,]  0.12582984  0.12690915  0.10732992  0.056033982 -0.11022969 -0.10760644
#> [5,] -0.21021691  0.08895420 -0.01431507  0.007530976 -0.06140415 -0.03089571
#> [6,]  0.13277350  0.02403049  0.10650773 -0.047333661 -0.15021738  0.14615809
#>             sPC7        sPC8        sPC9       sPC10
#> [1,] -0.11605651  0.06438109 -0.07995194  0.04794994
#> [2,]  0.04680788 -0.04979674 -0.04032454 -0.05890259
#> [3,] -0.02620345 -0.04184952 -0.06689604 -0.04091612
#> [4,]  0.02205297  0.10020362 -0.10742739 -0.13281197
#> [5,]  0.17344505  0.17249088  0.06990677  0.05819948
#> [6,] -0.07215218 -0.21001297 -0.11011887  0.10782760
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
#>             sPC1         sPC2         sPC3        sPC4        sPC5         sPC6
#> [1,] -0.09741275  0.119548690 -0.041750713  0.01558273 -0.06987547 -0.162089402
#> [2,]  0.05155618 -0.075864919 -0.024054478 -0.05160704 -0.18340225 -0.024822485
#> [3,]  0.17288032 -0.004374413 -0.003597050 -0.16205271  0.11203908  0.038604980
#> [4,] -0.13406621  0.017420922  0.006288485 -0.15228386  0.05294741  0.039407769
#> [5,]  0.13569768  0.161417350 -0.092539848  0.13218261  0.04738872  0.003145204
#> [6,] -0.02331851 -0.023092657  0.229797496  0.13150189  0.06103413  0.049204429
#>             sPC7        sPC8        sPC9       sPC10
#> [1,]  0.16421750  0.02700494 -0.07335335  0.60377438
#> [2,] -0.14628099  0.09698308 -0.11527438 -0.06997455
#> [3,]  0.10627597 -0.01203376 -0.09130578 -0.13897353
#> [4,] -0.03196648  0.14312004  0.15043683 -0.06856058
#> [5,] -0.08396788  0.05928110  0.08169624 -0.09136723
#> [6,]  0.01627797  0.08461162 -0.07247758 -0.05883859
```

# References

[1] A consistent pattern of slide effects in Illumina DNA methylation BeadChip array data. Hecker J, Lee S, Kachroo P, Prokopenko D, Maaser-Hecker A, Lutz SM, Hahn G, Irizarry R, Weiss ST, DeMeo DL, Lange C. Epigenetics. 2023 Dec;18(1):2257437. doi: 10.1080/15592294.2023.2257437. Epub 2023 Sep 20. PMID: 37731367

[2] Adjusting batch effects in microarray expression data using empirical Bayes methods. Johnson WE, Li C, Rabinovic A. Biostatistics. 2007 Jan;8(1):118-27. doi: 10.1093/biostatistics/kxj037. Epub 2006 Apr 21. PMID: 16632515

[3] The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD. Bioinformatics. 2012 Mar 15;28(6):882-3. doi: 10.1093/bioinformatics/bts034. Epub 2012 Jan 17. PMID: 22257669

[4] Letter to the editor: critical evaluation of the reliability of DNA methylation probes on the illumina MethylationEPIC v1.0 BeadChip microarrays. Hecker J, Weiss ST, Lasky-Su JA, DeMeo DL, Lange C. Epigenetics. 2024 Dec;19(1):2411470. doi: 10.1080/15592294.2024.2411470. Epub 2024 Oct 4. PMID: 39365898; PMCID: PMC11457593.
