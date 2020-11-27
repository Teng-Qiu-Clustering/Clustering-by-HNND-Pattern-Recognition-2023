#' @title Load data
#' @description load data
#'
#'
#' @param data_name
#'
#' @keywords load_data
#' @export
#' @examples
#'
load_data=function(data_name){
  if (any(data_name==c("cytof_h1","cytof_benchmark_h1"))) {
    data(cytof.benchmark.h1)
    x = cytof.benchmark.h1$x
    annotation_data = as.numeric(factor(cytof.benchmark.h1$metadata$cell))
  } else if (any(data_name==c("cytof_h2","cytof_benchmark_h2"))){
    data(cytof.benchmark.h2)
    x = cytof.benchmark.h2$x
    annotation_data = as.numeric(factor(cytof.benchmark.h2$metadata$cell))
  } else if (any(data_name==c("cytof_one","cytof_benchmark_one"))){
    data(cytof.benchmark.one)
    x = cytof.benchmark.one$x
    annotation_data = as.numeric(factor(cytof.benchmark.one$metadata$cell))
  } else if (any(data_name==c("single.cell.mrna.pollen","single_cell_mrna_pollen"))){
    data('single.cell.mrna.pollen')
    x = single.cell.mrna.pollen$x
    annotation_data = as.numeric(factor(single.cell.mrna.pollen$metadata$cell))
  } else if (data_name == 'bipolar'){
    # Macosko, E. Z. et al. Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets. Cell 161, 1202–1214 (2015).
    # pre-processed by https://bitbucket.org/jerry00/scvis-dev/src/master/data/
    load("data/bipolar_pca100.Rdata")
  } else if (data_name == 'PalmData25_uni'){
    load("data/PalmData25_uni.Rdata")
  } else if (any(data_name == c("Samusik_01","Samusik_all"))){
    print('this dataset has been normalized by asinh with a factor a 5')

    files <- list(
      Samusik_01   = "data/Samusik_01.fcs",
      Samusik_all  = "data/Samusik_all.fcs"
    )
    marker_cols <- list(
      Samusik_01   = 9:47,
      Samusik_all  = 9:47
    )

    x <- flowCore::exprs(flowCore::read.FCS(files[[data_name]], transformation = FALSE, truncate_max_range = FALSE))
    annotation_data<- x[, "label"]
    x <- x[, marker_cols[[data_name]]]
  }

  if (class(x) != "matrix"){
    x = as.matrix(x)
  }

  if (!is.null(annotation_data)){
    if (class(annotation_data) != 'matrix'){
      annotation_data = as.matrix(annotation_data)
    }
  }

  message("\n",data_name,": ", nrow(x)," rows (samples), ", ncol(x), " columns (features), ", length(unique(annotation_data))," clusters")

  return(list(x=x,annotation_data=annotation_data))

}
