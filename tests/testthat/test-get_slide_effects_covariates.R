library(SeffCovar)

test_that("check_approach_1", {
  methylation_matrix=matrix(rnorm(length(S_high)*100), nrow=length(S_high), ncol=100)
  rownames(methylation_matrix)=S_high
  slide_covars=get_slide_effects_covariates(methylation_matrix=methylation_matrix, input_set=S_high)
  expect_true(colnames(slide_covars)[1]=="sPC1")
  expect_false(any(is.na(slide_covars)))

})
