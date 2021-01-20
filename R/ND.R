#' @title Nearest Descent
#'
#' @description Nearest Descent (equivalent to Nearest Ascent)
#'
#'
#' @param F density vector
#' @param x test dataset
#' @param disName distanc measurement
#'
#' @return W:   edge weight vector;
#' @return I:   parent node vector;
#' @return F:   density vector;
#' @return root:   index of the root node (there is only one root node);
#'
#' @author Teng Qiu
#'
ND = function(F,x,disName = 'euclidean'){
  # ND_less_space_consumption
  N = length(F)
  W=vector("numeric",N)
  I=vector("numeric",N)
  Totoltime = vector("numeric",N)
  # print(N)
  for (i in seq(N)){
    # time_start<-Sys.time()
    # print(N)
    # print(paste(i,'/',N))
    Di=vector("numeric",N)
    Ci = F[i]-F
    Di[Ci>0]=Inf
    Di[(Ci==0)&(i>=seq(N))]=Inf
    candidate_id = which(!is.infinite(Di))
    for (j in 1:length(candidate_id)){
      if (disName == 'euclidean'){
        Di[candidate_id[j]] = sqrt(sum((x[i,]-x[candidate_id[j],])^2))
      } else if (disName == 'cosine'){
        Di[candidate_id[j]] = 1-sum(x[i,]*x[candidate_id[j],])/sqrt(sum(x[i,]^2)*sum(x[candidate_id[j],]^2))
      } else if (disName == 'l2'){
        Di[candidate_id[j]] = sum((x[i,]-x[candidate_id[j],])^2)
      } else {
        stop("the specified distance function is currently not supported")
      }
    }
    idx = which.min(Di)
    I[i] = idx
    W[i] = Di[idx]
    if (is.infinite(W[i])){
      W_Max_idx = i
      I[i]=i
      W[i]=0
    }
    # time_end<-Sys.time()

    # print(Totoltime[i]<-time_end-time_start)
    # if (i%%10 == 0){
    #   print(c("average time:",mean(Totoltime[1:10])))
    # }
  }
  return(list(W=W,I=I,root =W_Max_idx ))
}
