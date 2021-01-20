#' @title Nearest Neighbor Descent
#'
#' @description Nearest Neighbor Descent (equivalent to Nearest Neighbor Ascent)
#'
#'
#' @param F density vector
#' @param knn a list
#'
#' @return W:   edge weight vector;
#' @return I:   parent node vector;
#' @return roots:   indexes of the root nodes;
#' @return K_optimal:   which neighbor each parent node is;
#'
#' @author Teng Qiu
#'
NND = function(F,knn){
# NND_no_density_estimation_V1_memory_saving
  N = length(F)

  W = vector(length = N)
  I = vector(length = N)
  K_optimal = vector(length = N)
  for (i in 1:N){
    j = which(F[i]<F[knn$nn.index[i,]])
    if (length(j)==1){
      W[i] = knn$nn.dist[i,j]
      I[i] = knn$nn.index[i,j]
      K_optimal[i] = j
    } else if (length(j)>1){
      id = which.min(knn$nn.dist[i,j])
      j = j[id]
      W[i] = knn$nn.dist[i,j]
      I[i] = knn$nn.index[i,j]
      K_optimal[i] = j
    } else{
      W[i] = 0
      I[i] = i
      K_optimal[i] = N
    }
  }

  root_id = which(K_optimal==N)

  return(list(W=W,I=I,roots=root_id,K_optimal = K_optimal))
}
