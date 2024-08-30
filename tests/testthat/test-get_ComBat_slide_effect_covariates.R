library(SeffCovar)
test_that("check_approach_2", {
  methylation_matrix=matrix(rnorm(length(S_high)*100), nrow=length(S_high), ncol=100)
  rownames(methylation_matrix)=S_high
  slide=rep(paste0("slide",1:10),10)
  slide_covars=get_ComBat_slide_effect_covariates(methylation_matrix=methylation_matrix, slide=slide, mod=NULL, input_set=S_high)
  expect_true(colnames(slide_covars)[1]=="sPC1")
  expect_false(any(is.na(slide_covars)))
})
