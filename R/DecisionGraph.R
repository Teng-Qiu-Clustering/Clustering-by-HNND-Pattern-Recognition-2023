#' @title Decision Graph
#' @description Decision Graph for determing the inter-cluster edge nodes
#'
#' @param result
#'
#' @export result
#'
#' @examples
#'
DecisionGraph=function(result,Log = "xy"){

  result$peaks <- NA
  I = result$I_forInTree
  A = (result$Density_initial+result$Density_initial[I])/2

  B = result$W_forInTree


  if (Log == "y"){
    plot(A,B, xlab = "Mean Density",ylab = 'Edge Length (logscale)',
         cex = 0.7,log = 'y')
  } else if (Log == "x"){
    plot(A,B, xlab = "Mean Density (logscale)",ylab = 'Edge Length',
         cex = 0.7,log = 'x')
  } else if (Log == "xy"){
    plot(A,B, xlab = "Mean Density (logscale)",ylab = 'Edge Length (logscale)',
         cex = 0.7,log = 'xy')
  } else {
    plot(A,B, xlab = "Mean Density",ylab = 'Edge Length',
         cex = 0.7)
  }

  print('Click on plot to select several thresholds;push "finish" button at the top-right of the plot to end \n')
  threshold <- locator()
  print('get the coordinates of the mouse')
  peaks = vector()
  for (i in 1:length(threshold$x)){
    rho <- threshold$x[i]
    delta <- threshold$y[i]
    peaks <- c(peaks,which(A > rho &  B > delta))
  }
  result$peaks = unique(peaks)

  points(A[result$peaks],B[result$peaks], col = 2:(1 + length(result$peaks)),
         pch = 19,cex = 0.7)

  return(result)
}


