#' @title Constructing slide-effect-adjustment covariates
#'
#' @description This function constructs slide-effect-adjustment covariates from a beta-value or M-value DNA methylation matrix.
#' @param methylation_matrix. A matrix with DNA methylation measurements (beta-value or M-scale). Rows correspond to CpG sites and columns to samples. Rownames need to be in 'cgXXXXXX' format.
#' @param input_set. List of slide-effect-susceptible CpG sites from which the covariates will be extracted. Default is S_high.
#' @export
#' @examples
#' methylation_matrix=matrix(rnorm(length(S_high)*100), nrow=length(S_high), ncol=100)
#' rownames(methylation_matrix)=S_high
#' slide_covars=get_slide_effects_covariates(methylation_matrix=methylation_matrix, input_set=S_high)

get_slide_effects_covariates<-function(methylation_matrix, input_set=S_high)
{
		intersect_set=intersect(rownames(methylation_matrix), input_set)
		if(length(intersect_set)<10)
		{
		  stop("no sufficient overlap between provided CpG sites in methylation matrix and input set.")
		}
		methylation_matrix=methylation_matrix[rownames(methylation_matrix) %in% input_set,]
		mat=scale(t(methylation_matrix), center=TRUE, scale=TRUE) # scaling of CpG sites separately
		svd_mat=svd(mat)
		pcs=svd_mat$u[,1:10] # return first ten components
		
		return(pcs)
}