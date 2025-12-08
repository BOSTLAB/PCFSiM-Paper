#' CCC_computation
#'
#' Compute the Concordance Correlation Coefficient (CCC) between two numeric vectors.
#'
#' @param x A numeric vector.
#' @param y A numeric vector of the same length as x.
#'
#' @return A numeric value giving the CCC between x and y.
#' @export


CCC_computation = function(x,y) {
  return(2*sd(x,na.rm = T)*sd(y,na.rm = T)/(var(x,na.rm = T)+var(y,na.rm = T)+(mean(x)-mean(y))^2))
}
