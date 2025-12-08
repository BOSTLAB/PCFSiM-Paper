#' Fit Mixed Interaction Model to PCF Data
#'
#'
#' @param x A numeric vector representing the distance values (r) for the PCF.
#' @param y A numeric vector representing the PCF values (\code{g(r)}) corresponding to \code{x}.
#' @param show_plot A logical value. If \code{TRUE}, a plot will be generated. Default is \code{FALSE}.
#'
#' @return A named numeric vector containing the fitted parameters and statistics.
#' @importFrom minpack.lm nlsLM nls.lm.control
#' @importFrom stats coef residuals var
#'
#' @export

Fit_mixed_interaction = function(x,y,show_plot=FALSE) {
  y[0] = 1
  First_reach = min(which(y[-1]>=1))
  y[1:First_reach]=1
  data_temp = data.frame(x = x,y=y)
  data_temp = data_temp[-1,]
  
  curve.nlslrc = minpack.lm::nlsLM(y ~ (x^alpha)*exp(-(x*tau)^p)*constant+1,start = list(alpha=-0.5,tau=0.1,p=1,constant=max(y)),
                       data = data_temp,control = minpack.lm::nls.lm.control(maxiter=1000),lower = c(-0.9999,0,0,0))
  
  
  
  alpha = coef(curve.nlslrc)[1]
  tau = 1/coef(curve.nlslrc)[2]
  p = coef(curve.nlslrc)[3]
  constant  = coef(curve.nlslrc)[4]
  
  if (show_plot) {
    par(las=1,bty="l")
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1))
    curve( (x^alpha)*exp(-(x/tau)^p)*constant+1,add=T,lty=2,lwd=2,col="red")
    abline(h=1,lwd=1,lty=2,col='grey')
    legend("topright",legend = c(paste("Alpha=",round(alpha,digits = 2)),
                                 paste("tau=",round(tau,digits =2 )),
                                 paste("p=",round(p,digits = 2))),bty='n')
    characteristic_size = tau*gamma((alpha+1)/p)/gamma((alpha)/p)
    abline(v=characteristic_size,lty=2,col='green')
  }
  
  Normalizing_term = p/(tau^(alpha+1))/gamma((alpha+1)/p)
  Normalised_constant = constant / Normalizing_term # The constant need to be normalized relatively to the area under the curve of the distribution...
  
  r_squared = 1- mean(residuals(curve.nlslrc)^2)/var(y)
  Output_vector = c(alpha,tau,p,constant,Normalised_constant,r_squared)
  
  #tau
  
  names(Output_vector) = c("alpha","tau","p","C","C_normalised","R2")
  return(Output_vector)
}

