#' Fits the log laplace to Pcf(r)
#'
#' @param x A numeric vector.
#' @param y A numeric vector of the same length as x.
#' @param title Title for the plot.
#' @param show_plot Logical; whether to show the fitted curve.
#' @return A vector containing fitting parameters.
#' @export
#' 
Fit_log_Laplace= function(x,y,show_plot=FALSE,title=k) {
  y[0] = 1
  First_reach = min(which(y[-1]>=1))
  y[1:First_reach]=1
  y = y[x!=0]
  x = x[x!=0]
  
  data_temp = data.frame(x = x,y=y)
  
  curve.nlslrc = nlsLM(y ~ 1 + 1/(2*x*sigma)*exp(-log(abs(x-mu))/sigma)*constant,start = list(mu=1,sigma=1,constant=max(y)),
                       data = data_temp,control = nls.lm.control(maxiter=1000))
  mu = coef(curve.nlslrc)[1]
  sigma = coef(curve.nlslrc)[2]
  constant  = coef(curve.nlslrc)[3]
  
  if (show_plot) {
    par(las=1,bty="l")
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1),main = title)
    curve( 1 + 1/(2*x*sigma)*exp(-log(abs(x-mu))/sigma)*constant,add=T,lty=2,lwd=2,col="red")
    abline(h=1,lwd=1,lty=2,col='grey')
    legend("topright",legend = c(paste("mu=",round(mu,digits = 2)),
                                 paste("sigma=",round(sigma,digits = 2))),bty='n')
  }
  
  Normalizing_term = 1
  Normalised_constant = constant / Normalizing_term # The constant need to be normalized relatively to the area under the curve of the distribution...
  
  r_squared = 1- mean(residuals(curve.nlslrc)^2)/var(y)
  Output_vector = c(mu,sigma,constant,Normalised_constant,r_squared)
  
  names(Output_vector) = c("mu","sigma","C","C_normalised","R2")
  return(Output_vector)
}
