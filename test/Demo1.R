data_name = 'cytof.benchmark.h1';
data('cytof.benchmark.h1')
x = cytof.benchmark.h1$x

N = nrow(x)
K = ceiling(log2(N))

Totoltime<-system.time(
  result <- HNND_C(
    x = x,
    K = K,
    h = K,
    disName = "euclidean",
    LogPlot = '')
)[3]

# evaluate
annotation_data = as.numeric(factor(cytof.benchmark.h1$metadata$cell))

CScores = Evaluate_clustering_result(data_name = data_name,result$Label,annotation_data,EvaIndex = c('NMI','ARI'))

Sta_result = data.frame(
  N = N,
  Dim = ncol(x),
  NMI = round(CScores$NMI,digits = 3),
  ARI = round(CScores$ARI,digits = 3),
  RunTime = round(Totoltime-result$time_DG,digits = 3),
  ClusterNumber = length(unique(result$Label)))

print(Sta_result)
