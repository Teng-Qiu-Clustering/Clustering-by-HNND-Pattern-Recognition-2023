rm(list = ls())

# libraray ----------------------------------------------------------------

library("densitycut")
library("igraph")
library("RcppHNSW")
library("cytofkit")
library(gridExtra)
# sources -----------------------------------------------------------------
disName_array = c('euclidean')
data_names_ori = c('cytof_h1','cytof_h2','cytof_one','Samusik_01','Samusik_all')
compared_methods = c('phenograph','DensityCut',"NNA")
# initialization----

Totoltime = vector();NMI = vector();Cnum = vector();
Fscore = vector();MeanF = vector();ARI = vector();
meanPurity= vector();meanPurity = vector();result_disName = vector();
result_data = vector();result_k = vector();
result_method = vector();result_exponent = vector();
result_knn_method = vector();data_dim = vector()
result_k_all = data.frame()
result_memory = vector()
tt = 1
for (disName in disName_array){
  for (method_id in 1:length(compared_methods)){
    method_name = compared_methods[method_id]
    for (file_id in 1:length(data_names_ori)){
      Totoltime_per = 0

      data_name = data_names_ori[file_id]
      TestData=load_data(data_name)
      annotation_data = TestData$annotation_data
      x = TestData$x

      Dim = ncol(x)
      N = nrow(x)
      Dim = ncol(x)

      K_array = ceiling(log2(nrow(x)))
      ## density estimation method: simple_knn,st_knn, RW_knn
      densityMethod = 'simple_knn'

      if (Dim<=10){
        knn_method = 'kd_tree'
      } else {
        knn_method = "hnsw"
      }

      NumT = 50
      alpha = 0.9
      ClustN = length(unique(annotation_data))

      for (K_id in 1:length(K_array)){

        K = K_array[K_id]
        print(paste0("data:",data_name,", K:",K,", method:",method_name,", fastKNN:",knn_method))

        result_memory_per = pryr::mem_change(
          Totoltime_per<-system.time(

            if (any(method_name == c("NA+","NA++","NA+++ (ne.mem.)","NA+++ (ne.mem.) with cluster number","NA","NNA","NA+++ (with cluster number)"))){
              result <- fast_HNND_based_on_unit_free_index_v1_function(
                x = x,
                K = K_array[K_id],
                InTree_construction_method_id =  compared_methods[method_id],
                intercluster_edges_method_id = "D_ratio",
                densityMethod = densityMethod,
                Threshold_RootNum =  1,
                knn_method = knn_method,
                disName = disName,
                ClustN = ClustN,
                h = 50,
                thresh_k = 5,
                reuse_tree_contruction = "off")
              Label = result$Label
            } else if (method_name == "km (km++)"){
              result <- KMeans_rcpp(x, ClustN,initializer = "kmeans++")
              Label = result$clusters
            } else if (method_name == "km (random)"){
              result <- KMeans_rcpp(x, ClustN,initializer = "random")
              Label = result$clusters
            } else if (method_name == "single"){
              hc <- fastcluster::hclust.vector(x, "single")
              memb <- cutree(hc, k = ClustN)
              Label = memb
            } else if (method_name == "ward"){
              hc <- fastcluster::hclust.vector(x, "ward")
              memb <- cutree(hc, k = ClustN)
              Label = memb
            } else if (method_name == "DensityCut"){
              if (knn_method == 'kd_tree'){
                denC <- DensityCut(X=x, K=K, alpha=alpha, nu=seq(0.0, 1.0, by=0.05),
                                   debug=FALSE, show.plot=FALSE,
                                   maxit=5000)
                Label = denC$cluster
              } else if (knn_method == 'hnsw'){
                denC <- DensityCut_with_hnsw(X=x, K=K, alpha=alpha,disName =disName, nu=seq(0.0, 1.0, by=0.05),
                                             debug=FALSE, show.plot=FALSE,
                                             maxit=5000)
                Label = denC$cluster
              } else {
                stop('the fastKNN method is not supported')
              }
            } else if (method_name == "Hdbscan"){
              hdb <- hdbscan(x, minPts = K)
              Label = hdb$cluster
            } else if (method_name == "phenograph"){
              if (knn_method == 'kd_tree'){
                # cluster_res <- cytofkit::cytof_cluster(xdata = x, method = "Rphenograph")
                Rphenograph_out <- Rphenograph::Rphenograph(x,k=K)
                Label=as.numeric(membership(Rphenograph_out[[2]]))
              } else if (knn_method == 'hnsw'){
                cluster_res <- as.numeric(membership(Rphenograph_with_hnsw(data = x, k=K,disName =disName)))
                Label=cluster_res
              } else {
                stop('the fastKNN method is not supported')
              }
            } else if (method_name == "adpcluster"){
              ans <- adpclust(x, centroids = "auto", draw = FALSE,nclust = ClustN)
              Label = ans$clusters
            } else if (method_name == "meanshift (blur)"){
              clustering <- bmsClustering(t(x))
              Label = clustering$labels
            } else if (method_name == "meanshift (st)"){
              clustering <- msClustering(t(x))
              Label = clustering$labels
            } else if (method_name == "k-AP"){
              apres <- apcluster::apclusterK(apcluster::negDistMat(r=2), x, K=ClustN)
              Label = as.numeric(apres@idx)
            } else if (method_name == "DBSCAN"){
              db <- dbscan(x, eps = quantile(dist(x), 0.01))
            } else if (method_name == "PAM"){
              pm <- cluster::pam(dist(x),ClustN)
              Label = pm$clustering
            } else if (method_name == "GMM"){
              gmm = ClusterR::GMM(x, ClustN, "eucl_dist", "random_subset", 10, 10)
              pr = ClusterR::predict_GMM(x, gmm$centroids, gmm$covariance_matrices, gmm$weights)
              Label = pr$cluster_labels+1
            } else if (method_name == "clara"){
              clm = ClusterR::Clara_Medoids(x, clusters = ClustN, samples = 100,sample_size = 0.2, swap_phase = TRUE)
              Label = clm$clusters
            } else if (method_name == "SC"){
              sc <- kernlab::specc(x, centers= ClustN)
              Label = sc@.Data
            } else if (method_name == 'FlowSom'){
              set.seed(1000)
              out <- FlowSOM::ReadInput(flowCore::flowFrame(x), transform = FALSE, scale = FALSE)
              out <- FlowSOM::BuildSOM(out, colsToUse = 1:ncol(x))
              out <- FlowSOM::BuildMST(out)
              labels_pre <- out$map$mapping[, 1]

              out <- FlowSOM::metaClustering_consensus(out$map$codes, k = ClustN, seed = 1000)
              Label <- out[labels_pre]

            } else if (method_name == 'flowMeans'){
              set.seed(1000)
              out=flowMeans::flowMeans(flowCore::flowFrame(x),Standardize=FALSE)
              Label = out@Label
            } else if (method_name == 'flowPeak'){
              set.seed(1000)
              out <- flowPeaks::flowPeaks(x)
              Label = out[["peaks.cluster"]]
            }
          )[3]
        )
        print(paste("Totoltime_per: ",Totoltime_per,sep = ""))
        # print(paste("Totolmemory_per: ",total(result_memory_per),sep = ""))


        if (any(method_name == c('NA+++ (ne.mem.)','NA+++'))){
          k_all_per = data.frame("k_all" = result$k_all,"data_dim" = rep(Dim,N),"K_id" = rep(K,N),'method' = rep(method_name,N))
          result_k_all = rbind(result_k_all,k_all_per)
        }
        # evaluation --------------------------------------------------------------

        if (any(data_name == c("Levine_32dim","Levine_13dim","Samusik_01","Samusik_all","Nilsson_rare","Mosmann_rare"))){
          unassigned <- is.na(annotation_data)
          NMI[tt]=densitycut::ComputeNMI(annotation_data[!unassigned],Label[!unassigned])
          # res= helper_match_evaluate_multiple(Label[!unassigned], annotation_data[!unassigned])
          # MeanF[tt] = res$mean_F1
          # Fscore[tt]=FlowSOM::FMeasure(annotation_data[!unassigned],Label[!unassigned], silent = FALSE)
          ARI[tt] = mclust::adjustedRandIndex(annotation_data[!unassigned],Label[!unassigned])
          meanPurity[tt] =FlowSOM::Purity(annotation_data[!unassigned],Label[!unassigned], weighted = TRUE)[1]
        } else if (any(data_name == c("FlowCAP_ND","FlowCAP_WNV"))) {
          unassigned <- is.na(annotation_data)
          annotation_data = annotation_data[!unassigned]
          Label = Label[!unassigned]
          sam = TestData$individual_label[!unassigned]

          NMI_per=vector();MeanF_per = vector(); Fscore_per = vector(); ARI_per = vector(); meanPurity_per = vector()
          for (sam_u in 1:length(unique(sam))){
            IDX = which(sam == sam_u)
            NMI_per[sam_u]=densitycut::ComputeNMI(annotation_data[IDX],Label[IDX])
            # res= helper_match_evaluate_multiple(Label[IDX],annotation_data[IDX])
            # MeanF_per[sam_u] = res$mean_F1
            # Fscore_per[sam_u]=FlowSOM::FMeasure(annotation_data[IDX],Label[IDX], silent = FALSE)
            ARI_per[sam_u] = mclust::adjustedRandIndex(annotation_data[IDX],Label[IDX])
            meanPurity_per[sam_u] =FlowSOM::Purity(annotation_data[IDX],Label[IDX], weighted = TRUE)[1]
          }
          NMI[tt] = mean(NMI_per)
          MeanF[tt] = mean(MeanF_per)
          # Fscore[tt]= mean(Fscore_per)
          ARI[tt]=mean(ARI_per)
          meanPurity[tt] = mean(meanPurity_per)

        } else{
          NMI[tt]=densitycut::ComputeNMI(annotation_data,Label)
          # res= helper_match_evaluate_multiple(Label,annotation_data)
          # MeanF[tt] = res$mean_F1
          # Fscore[tt]=FlowSOM::FMeasure(annotation_data,Label, silent = FALSE)
          ARI[tt] = mclust::adjustedRandIndex(annotation_data,Label)
          meanPurity[tt] =FlowSOM::Purity(annotation_data,Label, weighted = TRUE)[1]
        }
        # for cluster number, we compute the result before deleting the unm-manual gated cells
        Cnum[tt] = length(unique(Label))
        Totoltime[tt] = Totoltime_per
        # result_memory[tt] = total(result_memory_per)
        result_data[tt] = data_name
        result_k[tt] = K
        data_dim[tt] = Dim
        result_method[tt] = method_name
        result_knn_method[tt] = knn_method
        result_disName[tt] = disName
        print(paste0("NMI score: ",round(NMI[tt],digits = 3)))
         tt = tt + 1

      }
    }
  }
}

# results all -----------------------------------------------------------------

print(data_names_ori)
print(compared_methods)
print(round(NMI,digits = 3))
print(round(ARI,digits = 3))
print(round(Totoltime,digits = 3))
print(Cnum)

# method 1: ggplot based displaying method
result_data_frame = data.frame(
  result_method = result_method,
  result_data = result_data,
  result_k = result_k,
  NMI = round(NMI,digits = 3),
  result_disName= result_disName,
  ARI=round(ARI,digits = 3),
  Cnum = Cnum,
  meanPurity = meanPurity,
  Totoltime = round(Totoltime,digits = 3),
  data_dim = paste("Sample dim.: ",data_dim,sep = ""),
  result_knn_method = result_knn_method)
print(result_data_frame)
write.csv(result_data_frame,file = paste0(method_name,"_result_data_frame_new.csv"), row.names=F)

save(file=paste("F:/Download_code/TPE/tpe/R/Results_Labels/",data_name,method_name,".Rdata",sep = ""),Label = Label)


