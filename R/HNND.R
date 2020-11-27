#' @title Hierachical Nearest Neighbor Descent
#'
#' @description Hierachical Nearest Neighbor Descent (equivalent to Hierachical Nearest Neighbor Ascent)
#'
#' @param x
#' @param K
#' @param h
#' @param disName
#' @param reuse_tree_contruction
#'
#' @export I,W,Label,Totoltime_DG
#' @export W
#' @export Label
#' @export Totoltime_DG
#'
#' @examples
#'
HNND = function(
  x,
  K = 20,
  h = 50,
  disName = "euclidean",
  LogPlot = ''){

  if (!any(disName == c('euclidean','cosine','l2'))) stop('distance metric not supported')

  if (ncol(x)<=10){
    knn_method = 'kd_tree'
  } else {
    knn_method = "hnsw"  # note: all the test datasets in HNND-PR have higher dimensions than 10
  }
  # HNND ---------------------------------------------------------------

  N = nrow(x)
  Dim = ncol(x)
  print(paste0('fast knn (',knn_method,')...'))

  switch(knn_method,
         "kd_tree" = {
           Kd_tree_construction <- WKNNF(x)
           knn <- Kd_tree_construction$query(x, k=K+1, eps=0, radius=0)
           names(knn)[1] <- "nn.index"
           names(knn)[2] <- "nn.dist"
           knn$nn.index = knn$nn.index[,-1]
           knn$nn.dist = knn$nn.dist[,-1]
         },
         "hnsw" = {
           ann = RcppHNSW::hnsw_build(x,disName)
           knn = queryKnnHNSW(ann,x,K+1,disName)
         }
  )

  print("density estimation...")
  Density_initial = simlified_knn_density_estimation_v1(knn,K)

  print("1st round NND:")
  nnd = NND(Density_initial,knn)

  I=nnd$I
  W=nnd$W

  root_id = nnd$W_Max_idx
  nnd_root_initial = nnd$W_Max_idx
  layer = 1

  print("other round NND:")

  while(1){

    if (length(root_id) <= h) break

    switch(knn_method,
           "kd_tree" = {
             Kd_tree_construction <- nabor::WKNNF(x[root_id,])
             knn <- Kd_tree_construction$query(x[root_id,], k=K+1, eps=0, radius=0)
             names(knn)[1] <- "nn.index"
             names(knn)[2] <- "nn.dist"
             knn$nn.index = knn$nn.index[,-1]
             knn$nn.dist = knn$nn.dist[,-1]
           },
           "hnsw" = {
             ann = RcppHNSW::hnsw_build(x[root_id,],disName)
             knn = RcppHNSW::hnsw_search(x[root_id,],ann,K+1)

             names(knn)[1] <- "nn.index"
             names(knn)[2] <- "nn.dist"
             knn$nn.index = knn$nn.index[,-1]
             knn$nn.dist = knn$nn.dist[,-1]
           }
    )

    if (min(Density_initial[root_id])!=max(Density_initial[root_id])){
      # nnd = NND_no_density_estimation_subset_v2(Density_initial,knn,K,root_id)
      nnd = NND(Density_initial[root_id],knn)
      I[root_id]=root_id[nnd$I]
      W[root_id]=nnd$W
      root_id = root_id[nnd$W_Max_idx]
    }else{
      I[root_id]=min(root_id)
      W[root_id]= apply(x[root_id,],1,function(z,y=x[min(root_id),]) (sum((z-y)^2))^0.5)
      root_id = min(root_id)
    }
    layer = layer + 1
  }

  Hnnd_root_id = root_id

  if (length(root_id) > 1){
    nd=ND(Density_initial[root_id],x[root_id,],disName)
    I[root_id]=root_id[nd$I]
    W[root_id]=nd$W
    Hnnd_root_id =which(I==c(1:N))
  }

  # Cut edges ---------------------------------------------------------------

  time_DG<-system.time({
    temp = list(W_forInTree = W,
                I_forInTree = I,
                Density_initial = Density_initial)
    result = DecisionGraph(temp,Log = LogPlot)
    peaks = result$peaks
  })[3]

  # search root -----

  I[peaks]=peaks
  print("dynamic update the parent node:")
  I_old=I
  I=I_old[I_old]
  while(sum(as.numeric(I-I_old))){
    I_old=I
    I=I_old[I_old]
  }
  Label = I

  return(list(I = I, W = W, Label=Label,time_DG = time_DG))
}

