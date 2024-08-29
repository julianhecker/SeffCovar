utils::globalVariables(c("S_high", "S_80_100"))

#' @title Constructing slide-effect-adjustment covariates using SVD/Eigenvalue decomposition
#'
#' @description This function constructs slide-effect-adjustment covariates from a beta-value or M-value DNA methylation matrix.
#' @param methylation_matrix A matrix with DNA methylation measurements (beta-value or M-scale). Rows correspond to CpG sites and columns to samples. Rownames need to be in 'cgXXXXXX' format.
#' @param input_set List of slide-effect-susceptible CpG sites from which the covariates will be extracted. Default is S_high.
#' 
#' @return pcs Top ten slide effect PCs (
#' @export
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
		colnames(pcs)=paste0("sPC",1:10)
		return(pcs)
}



#' @title Constructing slide-effect-adjustment covariates using SVD/Eigenvalue decomposition on ComBat-estimated slide effects
#'
#' @description This function uses a modification of the ComBat implementation in the sva R package. The estimated gamma.hat (least squares estimate) are used to perform a
#' singular value decomposition.
#' The original code for ComBat from the sva package that can be found at
#' https://bioconductor.org/packages/release/bioc/html/sva.html 
#' The original code as well as this code here is distributed under the Artistic License 2.0.
#'
#' @param methylation_matrix A matrix with DNA methylation measurements (beta-value or M-scale). Rows correspond to CpG sites and columns to samples. Rownames need to be in 'cgXXXXXX' format.
#' @param slide Slide information for each sample in the dataset
#' @param mod Covariate information to adjust, default is NULL.    
#' @param input_set List of slide-effect-susceptible CpG sites from which the covariates will be extracted. Default is S_high.
#'
#' @return res List of objects corresponding to the singular value decomposition, including the top ten PCs
#'
#' @importFrom stats model.matrix var
#'
#'
#' @export
get_ComBat_slide_effect_covariates <- function(methylation_matrix, slide, mod = NULL,  input_set=S_high) {

    if(length(dim(slide))>1){
      stop("Slide information needs to be 1-dimensional")
    }   
	intersect_set=intersect(rownames(methylation_matrix), input_set)
	if(length(intersect_set)<10)
	{
		stop("no sufficient overlap between provided CpG sites in methylation matrix and input set.")
	}
	methylation_matrix <- as.matrix(methylation_matrix[rownames(methylation_matrix) %in% S_high,])
	if(any(is.na(methylation_matrix))){
	 stop("NAs found, please exclude any CpGs with missing values.")
	}
	
    
    ## find zero variance CpGs
    slide <- as.factor(slide)
    zero.rows.lst <- lapply(levels(slide), function(slide_id){
      if(sum(slide==slide_id)>1){
        return(which(apply(methylation_matrix[, slide==slide_id], 1, function(x){var(x)==0})))
      }else{
        return(which(rep(1,3)==2))
      }
    })
    zero.rows <- Reduce(union, zero.rows.lst)
    keep.rows <- setdiff(1:nrow(methylation_matrix), zero.rows)
    
    if (length(zero.rows) > 0) {
      methylation_matrix.orig <- methylation_matrix
      methylation_matrix <- methylation_matrix[keep.rows, ]
    }
  

    if(any(table(slide)==1)){stop("please provide data with at least 2 samples per slide.")}
    
    
    slidemod <- model.matrix(~-1+slide)  
    
    cat("Found",nlevels(slide),"slides\n")
	
	## A few other characteristics on the batches
    n.slide <- nlevels(slide)
    slides <- list()
    for (i in 1:n.slide) {
        slides[[i]] <- which(slide == levels(slide)[i])
    } # list of samples in each batch  
    n.slides <- sapply(slides, length)
    
    n.array <- sum(n.slides)
	
    ## combine slide variable and covariates
    design <- cbind(slidemod, mod)
  
    ## check for intercept in covariates, and drop if present
    check <- apply(design, 2, function(x) all(x == 1))
    design <- as.matrix(design[,!check])
  
    ## Number of covariates or covariate levels
    cat("Adjusting for", ncol(design)-ncol(slidemod), 'covariate(s) or covariate level(s)\n')
  
    ## Check if the design is confounded
    if(qr(design)$rank < ncol(design)) {
        ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
        if(ncol(design)==(n.slides+1)) {
            stop("The covariate is confounded with slide information! Remove the covariate and rerun.")
        }
        if(ncol(design)>(n.slides+1)) {
            if((qr(design[,-c(1:n.slides)])$rank<ncol(design[,-c(1:n.slides)]))){
                stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded.')
            } else {
                stop("At least one covariate is confounded with slide information! Please remove confounded covariates and rerun.")
            }
        }
    }
  
   
    
    
    B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(methylation_matrix)))
    
    grand.mean <- crossprod(n.slides/n.array, B.hat[1:n.slide,])
    var.pooled <- ((methylation_matrix-t(design %*% B.hat))^2) %*% rep(1/n.array,n.array) # FIXME
    stand.mean <- t(grand.mean) %*% t(rep(1,n.array)) # FIXME
    if(!is.null(design)){
        tmp <- design
        tmp[,c(1:n.slide)] <- 0
        stand.mean <- stand.mean+t(tmp %*% B.hat) #FIXME
    }  
    methylation_matrix <- (methylation_matrix-stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array))) # FIXME
  
    ##Get regression batch effect parameters
    cat("Fitting L/S model and finding priors\n")
    slide.design <- design[, 1:n.slide]
    
    gamma.hat <- solve(crossprod(slide.design), tcrossprod(t(slide.design),
                                                             as.matrix(methylation_matrix)))
    
	#############################################################
    random_effects_mat= t(gamma.hat) %*% t(slide.design)
    mat=scale(t(random_effects_mat), center=FALSE, scale=FALSE) # scaling of CpG sites separately
    svd_mat=svd(mat)
    pcs=svd_mat$u[,1:10] # return first ten components
    colnames(pcs)=paste0("sPC",1:10)
  
    return(res=list('random_effects_mat'=random_effects_mat, 'pcs'=pcs, 'svs'=svd_mat$d[1:50], 'v'=svd_mat$v[,1]))
  
}