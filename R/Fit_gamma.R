#' Fits the gamma model to Pcf(r)
#'
#' @param x A numeric vector.
#' @param y A numeric vector of the same length as x.
#' @param show_plot Logical; whether to show the fitted curve.
#' @param title Title for the plot.
#' @return A vector containing fitting parameters.
#' @export
#' 
Fit_gamma = function(x,y,show_plot=FALSE,title=k) {
  y[0] = 1
  First_reach = min(which(y[-1]>=1))
  y[1:First_reach]=1
  
  
  data_temp = data.frame(x = x,y=y)
  data_temp = data_temp[-1,]
  
  curve.nlslrc = nlsLM(y ~ (x^alpha)*exp(-(x*tau))*constant+1,start = list(alpha=2,tau=0.01,constant=100),
                      data = data_temp,control = nls.lm.control(maxiter=1000),lower = c(-0.999,0,0))

  alpha = coef(curve.nlslrc)[1]
  tau = 1/coef(curve.nlslrc)[2]
  constant  = coef(curve.nlslrc)[3]
  
  if (show_plot) {
    par(las=1,bty="l")
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1),main = title)
    curve( (x^alpha)*exp(-(x/tau))*constant+1,add=T,lty=2,lwd=2,col="red")
    abline(h=1,lwd=1,lty=2,col='grey')
    legend("topright",legend = c(paste("Alpha=",round(alpha,digits = 2)),
                                 paste("tau=",round(tau,digits =2 ))),bty='n')
    characteristic_size =  (alpha+1)*tau
    abline(v=characteristic_size,lty=2,col="green")
  }
  
  Normalizing_term = 1/(tau^(alpha+1))/gamma((alpha+1))
  Normalised_constant = constant / Normalizing_term # The constant need to be normalized relatively to the area under the curve of the distribution...
  
  r_squared = 1- mean(residuals(curve.nlslrc)^2)/var(y)
  Output_vector = c(alpha,tau,constant,Normalised_constant,r_squared)
  
  names(Output_vector) = c("alpha","tau","C","C_normalised","R2")
  return(Output_vector)
}

