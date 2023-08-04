#' @title Constructing slide-effect-adjustment covariates
#'
#' @description This function constructs slide-effect-adjustment covariates from a beta-value or M-value DNA methylation matrix.
#' @param methylation_matrix. A matrix with DNA methylation measurements (beta-value or M-scale). Rows correspond to CpG sites and columns to samples.
#' @param input_set. List of slide-effect-susceptible CpG sites from which the covariates will be extracted
#' @export
#' @examples
#' tobedone

get_slide_effects_covariates<-function(methylation_matrix, input_set=S_high)
{
		intersect_set=intersect(rownames(methylation_matrix), input_set)
		if(length(intersect_set)<10)
		{
		  stop("no sufficient overlap between provided CpG sites and S_80_100 or S_high")
		}
		methylation_matrix=methylation_matrix[rownames(methylation_matrix) %in% input_set,]
		mat=scale(t(methylation_matrix), center=TRUE, scale=TRUE)
		svd_mat=svd(mat)
		pcs=svd_mat$u[,1:10]
		
		return(pcs)
}