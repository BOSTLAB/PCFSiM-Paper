#' Fits the Dagum model to Pcf(r)
#'
#' @param x A numeric vector.
#' @param y A numeric vector of the same length as x.
#' @param title Title for the plot.
#' @param show_plot Logical; whether to show fitted curve.
#' @return A vector containing fitting parameters. 
#' @export
#' 
Fit_Dagum = function(x,y,show_plot=FALSE,title=k) {
  y[0] = 1
  First_reach = min(which(y[-1]>=1))
  y[1:First_reach]=1
  
  
  data_temp = data.frame(x = x,y=y)
  
  curve.nlslrc = nlsLM(y ~ 1 + constant*(1+exp(x*tau)^-alpha)^-beta,start = list(tau=0.1,alpha=1,beta=1,constant=100),
                       data = data_temp,control = nls.lm.control(maxiter=1000),lower = c(0,0,0,0))
  tau = 1/coef(curve.nlslrc)[1]
  alpha = coef(curve.nlslrc)[2]
  beta = coef(curve.nlslrc)[3]
  constant  = coef(curve.nlslrc)[4]
  
  if (show_plot) {
    par(las=1,bty="l")
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1),main = title)
    curve( 1 + constant*(1+exp(x/tau)^-alpha)^-beta,add=T,lty=2,lwd=2,col="red")
    abline(h=1,lwd=1,lty=2,col='grey')
    legend("topright",legend = c(paste("Tau=",round(tau,digits = 2)),
                                 paste("Alpha=",round(alpha,digits = 2)),
                                 paste("Beta=",round(beta,digits =2 ))),bty='n')
  }
  
  Normalizing_term = 1
  Normalised_constant = constant / Normalizing_term # The constant need to be normalized relatively to the area under the curve of the distribution...
  
  r_squared = 1- mean(residuals(curve.nlslrc)^2)/var(y)
  Output_vector = c(tau,alpha,beta,constant,Normalised_constant,r_squared)
  
  names(Output_vector) = c("tau","alpha","beta","C","C_normalised","R2")
  return(Output_vector)
}

