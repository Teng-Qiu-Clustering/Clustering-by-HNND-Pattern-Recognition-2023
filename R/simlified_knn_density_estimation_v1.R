#' @title density estimation
#'
#' @description simlified knn density estimation (v1)
#'
#'
#' @param knn
#' @param k
#'
#' @export P_Amp
#' @examples
#'
simlified_knn_density_estimation_v1=function(knn,k){
  # note that the motivation for this samplification are that 1) since we only care the relative density values of the two nodes,
  # 2) the distances of each node to the neighbors are ranked in ascending order in knn$nn.dist

  Radius = apply(knn$nn.dist,1,max)
  P_Amp = 1/Radius

  # Radius = knn$nn.dist[,k]
  # P_Amp = 1/Radius

  idx = which(is.infinite(P_Amp))
  if (length(idx)!=0){
    P_Amp[idx] = max(P_Amp[-idx])
  }


  return(P_Amp=P_Amp)
}
