#' Fits the exponential model to Pcf(r)
#'
#' @param x A numeric vector.
#' @param y A numeric vector of the same length as x.
#' @param title Title for the plot.
#' @param show_plot Logical; whether to show the fitted curve.
#' @return A vector containing fitting parameters.
#' @export
#' 
Fit_exponential = function(x,y,show_plot=FALSE,title=k) {
  y[0] = 1
  First_reach = min(which(y[-1]>=1))
  y[1:First_reach]=1
  data_temp = data.frame(x = x,y=y)
  data_temp = data_temp[-1,]
  data_temp$y = data_temp$y - 1
  m = lm(log(data_temp$y)~1+data_temp$x)
  
  tau_estimate = coef(m)[2]
  C_estimate = exp(coef(m)[1])
  
  curve.nlslrc = nlsLM(y ~ exp(x*tau)*constant,start = list(tau=tau_estimate,constant=C_estimate),
                       data = data_temp,control = nls.lm.control(maxiter=1000))
  
  tau = -1/coef(curve.nlslrc)[1]
  constant  = coef(curve.nlslrc)[2]
  
  if (show_plot) {
    par(las=1,bty="l")
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1),main = title)
    curve(  exp(-x/tau)*constant+1,add=T,lty=2,lwd=2,col="red")
    abline(h=1,lwd=1,lty=2,col='grey')
    legend("topright",legend = c(paste("tau=",round(tau,digits = 2))),bty='n')
    characteristic_size = 1/tau
    abline(v=characteristic_size,lty=2,col='green')
  }
  
  Normalizing_term = tau
  Normalised_constant = constant / Normalizing_term # The constant need to be normalized relatively to the area under the curve of the distribution...
  
  r_squared = 1- mean(residuals(curve.nlslrc)^2)/var(y)
  Output_vector = c(tau,constant,Normalised_constant,r_squared)
  
  #tau
  
  names(Output_vector) = c("tau","C","C_normalised","R2")
  return(Output_vector)
}

