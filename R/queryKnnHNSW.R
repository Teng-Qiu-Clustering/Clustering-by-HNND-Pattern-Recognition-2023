#' @title Find k-nearest-neighbors (by HNSW)
#'
#' @description get the k-nearest-neighbors of the query points
#'
#' @param ann
#' @param x
#' @param K
#' @param disName
#'
#' @export nn.index
#' @export nn.dist
#'
#' @examples
#'
queryKnnHNSW = function(ann,x,K,disName) {
  # You can get distances with:
  # res <- p$getAllNNsList(data, k = 1, include_distances = TRUE)
  # res$dist contains the distance matrix, res$item stores the indexes

  knn_res <- ann$getAllNNsList(x, K, include_distances = TRUE)

  nn.dist = knn_res$dist[,-1]

  if (disName == "euclidean"){
    ## Note: the distance output of "ann$getAllNNsList"
    ## is the square of that of function "hnsw_search(x,ann,K+1)" in "hnsw" package
    nn.dist = sqrt(nn.dist)
  }

  nn.index = knn_res$item[,-1]

  if (!is.matrix(nn.index)){
    nn.index = as.matrix(nn.index,length(nn.index),1)
    nn.dist = as.matrix(nn.dist,length(nn.dist),1)
  }
  return(list(nn.index=nn.index,
              nn.dist =nn.dist)
  )
}
