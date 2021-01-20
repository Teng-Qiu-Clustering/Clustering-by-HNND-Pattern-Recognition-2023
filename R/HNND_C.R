##' @title Clustering by Hierachical Nearest Neighbor Descent
##'
##' @description Clustering by Hierachical Nearest Neighbor Descent (equivalent to: Clustering by Hierachical Nearest Neighbor Ascent)
##'
##' @param x dataset (a matrix with rows denoting samples and columns denoting features
##' @param K number of nearest neighbors (emperical value: the nearest integer to log2(N), where N is the size of the dataset)
##' @param h a threshold for deciding when to start the ND layer (usually h = k)
##' @param LogPlot whether to show the decision graph in the log scale
##' @param disName distance measurement ('euclidean','cosine', and'l2' are supported)
##'
##' @return I:    parent node vector;
##' @return W:    edge weigth vector;
##' @return Label:    cluster label vector;
##' @return Totoltime_DG:     time spent on the interactive operation on the Decision Graph (which is user-dependent and is thus suggested to be subtracted from the total runtime);
##'
##'
##' @examples
##' data('cytof.benchmark.h1')
##' x = cytof.benchmark.h1$x
##' K = ceiling(log2(nrow(x)))
##' result <- HNND_C(x = x,K = K,h = K,disName = "euclidean",LogPlot = '')
##'
##' @author Teng Qiu
HNND_C = function(
  x,
  K = 20,
  h = 20,
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
  Density_initial = simlified_knn_density_estimation_v1(knn)

  print("1st round NND:")
  nnd = NND(Density_initial,knn)

  I=nnd$I
  W=nnd$W

  root_id = nnd$roots
  nnd_root_initial = nnd$roots
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
      root_id = root_id[nnd$roots]
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
    peaks = DecisionGraph(W,I,Density_initial,Log = LogPlot)
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

