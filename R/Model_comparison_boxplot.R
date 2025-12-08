#' Makes a boxplot comparing the R2 of different Pcf fitted model
#'@param list_results A named list where each element is a \code{data.frame} or
#'                     a similar structure containing the fitting results for a specific model.
#'                     Each element **must** contain a column named \code{R2} storing the
#'                     coefficient of determination for each fit performed.
#' @return NULL This function does not return a value; it displays a plot.
#' @importFrom RColorBrewer brewer.pal
#' @export

Model_comparison_boxplot = function(list_results=list(Model_Sigmoid = Results_sigmoid,Model_exponential = Results_exponential)) {
  names = names(list_results)
  par(las=1,bty='l')
  boxplot(lapply(list_results, function(x) x$R2), ylim=c(0.4,1),xaxs='i',yaxs='i',names=names,
          xlab="Model",ylab="R2",cex.lab=1.3,col=RColorBrewer::brewer.pal(length(list_results),"Spectral"),cex.axis=0.7)
  return(invisible())
}
