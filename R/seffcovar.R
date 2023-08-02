#' @title Constructing slide-effect-adjustment covariates
#'
#' @description This function constructs slide-effect-adjustment covariates from a beta-value or M-value DNA methylation matrix.
#' @param methylation_matrix. A matrix with DNA methylation measurements (beta-value or M-scale). Rows correspond to CpG sites and columns to samples.
#' @param scale_parameter. Scale parameter for the prcomp function, defaults to FALSE.
#' @export
#' @examples
#' tobedone

get_slide_effects_covariates<-function(methylation_matrix, scale_parameter=FALSE)
{
		
		intersect_high=intersect(rownames(methylation_matrix), S_high)
		intersect_80_100=intersect(rownames(methylation_matrix), S_80_100)
		if(length(intersect_high)<10 | length(intersect_80_100)<10)
		{
		  stop("no sufficient overlap between provided CpG sites and S_80_100 or S_high")
		}
		
		methylation_matrix_high=methylation_matrix[rownames(methylation_matrix) %in% S_high,]
		methylation_matrix_80_100=methylation_matrix[rownames(methylation_matrix) %in% S_80_100,]
		
		pca_high=prcomp(methylation_matrix_high, center=TRUE, scale=scale_parameter)
		pca_80_100=prcomp(methylation_matrix_80_100, center=TRUE, scale=scale_parameter)
		
		pcs_high=pca_high$rotation[,1:10]
		pcs_80_100=pca_80_100$rotation[,1:10]
		
		vars_high=pca_high$sdev[1:10]**2/sum(pca_high$sdev**2)
		vars_80_100=pca_80_100$sdev[1:10]**2/sum(pca_80_100$sdev**2)
		
		return(list(pcs_high=pcs_high, pcs_80_100=pcs_80_100, vars_high=vars_high,
		vars_80_100=vars_80_100))

}