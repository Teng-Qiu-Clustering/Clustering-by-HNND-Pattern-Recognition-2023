#' @title density estimation
#'
#' @description simlified knn density estimation (v1)
#'
#'
#' @param knn a list (knn$nn.dist: knn distance matrix)
#'
#' @return F: density vector
#'
simlified_knn_density_estimation_v1=function(knn){

  Radius = apply(knn$nn.dist,1,max)
  F = 1/Radius

  idx = which(is.infinite(F))
  if (length(idx)!=0){
    F[idx] = max(F[-idx])
  }


  return(F=F)
}
