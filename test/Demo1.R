rm(list = ls())
library('HnndPR')
# datasets ----------------------------------------------------------------
 
data_names_ori = 'cytof_h1'


result_all = data.frame()
for (data_name in data_names_ori){
  #  load data ----
  TestData=load_data(data_name)
  N = nrow(TestData$x)
  K_empircal = ceiling(log2(N))

  if (any(data_name==c('bipolar','single.cell.mrna.pollen','PalmData25_uni'))) {
    disName = 'cosine'
    LogPlot = 'x'
  } else {
    disName = 'euclidean'
    LogPlot = ''
  }

  Totoltime<-system.time(
    result <- HNND(
      x = TestData$x,
      K = K_empircal,
      h = K_empircal,
      disName = disName,
      LogPlot = LogPlot)
  )[3]

  # evaluation --------------------------------------------------------------
  CScores = Evaluate_clustering_result(data_name,result$Label,TestData$annotation_data,EvaIndex = c('NMI','ARI'))

  result_all = rbind(result_all,data.frame(
    DataName = data_name,
    N = N,
    Dim = ncol(TestData$x),
    NMI = round(CScores$NMI,digits = 3),
    ARI = round(CScores$ARI,digits = 3),
    Totoltime = round(Totoltime-result$time_DG,digits = 3),
    Cnum = length(unique(result$Label))))
}
print(result_all)
