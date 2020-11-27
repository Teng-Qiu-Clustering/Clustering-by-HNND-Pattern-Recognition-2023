#' @title Evaluate clustering results
#' @description Evaluate clustering results
#' my_evaluate_clustering_result
#'
#'
#' @param data_name
#' @param Label
#' @param annotation_data
#' @param EvaIndex
#'
#' @param Label,annotation_data,EvaIndex
#'
#' @examples
#'
Evaluate_clustering_result = function(data_name,Label,annotation_data,EvaIndex = c('NMI','ARI','MeanFscore')){
  # Note 1 (for NMI_max_version): NMI_max_version = aricode::NMI(A,B,'max') (Note that aricode::AMI(A,B) is problematic, since it outputs the same result as aricode::NMI(A,B))
  # this equals to the NMI output by matlab code "ANMI_analytical_11" written by  Nguyen Xuan Vinh 2008-2009
  # This version is however much faster than than those in "ANMI_analytical_11" and "aricode"
  #References:
  #  'Information Theoretic Measures for Clusterings Comparison: Variants, Properties, Normalization and Correction for Chance', N.X. Vinh, Epps, J. and Bailey, J., Journal of Machine Learning Research, 2010

  # Note 2 (for densitycut::ComputeNMI):
  # densitycut::ComputeNMI = SNFtool::calNMI = aricode::NMI(A,B,'sqrt') ~= mclustcomp::mclustcomp(A,B,types = c('nmi1'))
  # the methods in note 2 refer to 'nmi1'	Normalized Mutual Information by Strehl and Ghosh.

  # Note 3: aricode contains 5 variants of NMI: NMI(c1, c2, variant = c("max", "min", "sqrt", "sum", "joint"))
  # mclustcomp contains 3 variants of NMI

  # Note 4: aricode::NMI(A,B,'sum') = mclustcomp::mclustcomp(A,B,types = 'nmi3')
  # Note 5: aricode::NVI(A,B) != mclustcomp::mclustcomp(A,B,types = 'nvi'), which is strange and problematic
  # Note 6: mclustcomp::mclustcomp(A,B,types = 'adjrand') = aricode::ARI(A,B) = mclust::adjustedRandIndex(A,B)

  NMI= c();NMI_max_version = c();ARI = c();MeanF = c()
  Label = as.numeric(Label)
  if (!is.numeric(annotation_data)){
    annotation_data = factor(annotation_data,labels =  1:length(unique(annotation_data)))
    annotation_data = as.numeric(annotation_data)
  }
  annotation_data = as.numeric(annotation_data)
  if (any(data_name == c("Samusik_01","Samusik_all"))){
    unassigned <- is.na(annotation_data)
    annotation_data = annotation_data[!unassigned]
    Label = Label[!unassigned]
    if (any('NMI' == EvaIndex)){
      NMI=densitycut::ComputeNMI(annotation_data,Label)
    }
    if (any('NMI_max_version' == EvaIndex)){
      NMI_max_version = NMI_max_version(annotation_data,Label)
    }
    if (any('MeanFscore' == EvaIndex)){
      res= helper_match_evaluate_multiple(annotation_data, Label)
      MeanF = res$mean_F1
      # Fscore=FlowSOM::FMeasure(annotation_data,Label, silent = FALSE)
    }
    if (any('ARI' == EvaIndex)){
      ARI = mclust::adjustedRandIndex(annotation_data,Label)
    }

    # meanPurity =FlowSOM::Purity(annotation_data,annotation_data, weighted = TRUE)[1]
  } else {
    if (any('NMI' == EvaIndex)){
      NMI=densitycut::ComputeNMI(annotation_data,Label) # MI(A,B)/sqrt(H(A)*H(B)) reference: 'nmi1'	Normalized Mutual Information by Strehl and Ghosh.
    }
    if (any('NMI_max_version' == EvaIndex)){
      NMI_max_version = NMI_max_version(annotation_data,Label)
    }
    if (any('MeanFscore' == EvaIndex)){
      res= helper_match_evaluate_multiple(Label,annotation_data)
      MeanF = res$mean_F1
      # Fscore=FlowSOM::FMeasure(annotation_data,Label, silent = FALSE)
    }
    if (any('ARI' == EvaIndex)) {ARI = mclust::adjustedRandIndex(annotation_data,Label)}

    # meanPurity =FlowSOM::Purity(annotation_data,Label, weighted = TRUE)[1]
  }

  ## Note that one can use the following code to show NMI and ARI
  # E2 = mclustcomp::mclustcomp(A,B,types = c('nmi1','adjrand'))
  # mclustcomp package collects 25 evaluation index

  return(list(NMI= NMI,NMI_max_version = NMI_max_version,ARI = ARI,MeanF = MeanF))

}
