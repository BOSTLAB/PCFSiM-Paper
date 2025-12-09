library(spatstat)
library(balagan)
library(minpack.lm)
library(nls.multstart)


#II)Defining functions


color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("grey","yellow","orange","red","darkred"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.99,na.rm = TRUE)
  }
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=x_prime/max_scale
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}

Compute_pcf = function(sce,N_min_points=50,  r_vector = seq(0,500,length.out = 100),computation_method = "direct",verbose=FALSE) {
 
  List_clusters = unique(colLabels(sce))
  N_ROI = length(unique(sce$ImageNumber))
  List_pcf_functions = list()
  PP_cell_type_annotation = c()
  PP_FoV_annotation = c()
  
  if (!computation_method%in%c("direct","derivative")){
    stop("Please provide a correct pcf computation method: has to be direct of derivative")
  }
  
  i = 1 
  for (j in List_clusters) {
    if (verbose) {
      print(j)
    }

    for (k in unique(sce$ImageNumber)) {
      sce_temp = sce[,sce$ImageNumber==k & colLabels(sce) %in% j]
      #print(k)
      if (ncol(sce_temp)>=N_min_points) {
        
        ppp_temp = ppp(sce_temp$Location_Center_X,sce_temp$Location_Center_Y,window = owin(xrange = range(sce_temp$Location_Center_X),yrange = range(sce_temp$Location_Center_Y)))
        
        if (computation_method=='derivative') {
          Kest_temp = Kest.fft(ppp_temp,r = r_vector,rmax = max(r_vector),sigma = 1)


          if (sum(!is.na(Kest_temp$border))<=4){
            Kest_temp$border[is.na(Kest_temp$border)] = 1
          }
          
          pcf_temp = pcf.fv(Kest_temp,method = "b") #Pcf is computed by taking the derivative of the K function and by setting pcf(0)=0 as two cells cannot be at the same place 
        }
        
        if (computation_method=='direct') {
          pcf_temp = pcf.ppp(ppp_temp,r = r_vector,rmax = max(r_vector),kernel = "rectangular",correction = "translate",divisor = "r")
        }

        List_pcf_functions[[i]] = pcf_temp
        
        PP_cell_type_annotation = c(PP_cell_type_annotation,j)
        PP_FoV_annotation = c(PP_FoV_annotation,k)
        i = i + 1 
        
      }
      
    }
    
  }
  
  if (computation_method=='derivative') {
    List_pcf_values = lapply(List_pcf_functions,FUN = function(x) {x$pcf})
  }
  if (computation_method=='direct') {
    List_pcf_values = lapply(List_pcf_functions,FUN = function(x) {x$trans})
  }
  
  
  
  List_r_values = lapply(List_pcf_functions,FUN = function(x) {x$r})
  return(list(List_pcf =List_pcf_values,List_r = List_r_values,Annotation = data.frame(ROI =PP_FoV_annotation,Cluster =PP_cell_type_annotation )))
}

Fit_power_law = function(x,y,show_plot=FALSE) {
  y[0] = 1
  First_reach = min(which(y[-1]>=1))
  y[1:First_reach]=1
  data_temp = data.frame(x = x,y=y)
  data_temp = data_temp[-1,]

  curve.nlslrc = nlsLM(y ~1 + (x/tau)^-alpha,start = list(tau=100,alpha=2),
                       data = data_temp,control = nls.lm.control(maxiter=1000),lower = c(0,0))
  tau = coef(curve.nlslrc)[1]
  alpha = coef(curve.nlslrc)[2]
  
  if (show_plot) {
    par(las=1,bty="l")
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1))
    curve( (x/tau)^-alpha+1,add=T,lty=2,lwd=2,col="red")
    abline(h=1,lwd=1,lty=2,col='grey')
    legend("topright",legend = c(paste("Alpha=",round(alpha,digits = 2)),
                                 paste("tau=",round(tau,digits =2 ))),bty='n')
    abline(h=1,lty=2,col="grey")
  }
  r_squared = 1- mean(residuals(curve.nlslrc)^2)/var(y)
  Output_vector = c(alpha,tau,r_squared)
  
  names(Output_vector) = c("alpha","tau","R2")
  return(Output_vector)
  
}

Fit_generalized_gamma = function(x,y,show_plot=FALSE) {
  y[0] = 1
  First_reach = min(which(y[-1]>=1))
  y[1:First_reach]=1
  data_temp = data.frame(x = x,y=y)
  data_temp = data_temp[-1,]
  
  curve.nlslrc = nlsLM(y ~ (x^alpha)*exp(-(x*tau)^p)*constant+1,start = list(alpha=-0.5,tau=0.1,p=1,constant=max(y)),
                       data = data_temp,control = nls.lm.control(maxiter=1000),lower = c(-0.9999,0,0,0))
  


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

Fit_exponential = function(x,y,show_plot=FALSE) {
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
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1))
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


Fit_gamma = function(x,y,show_plot=FALSE) {
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
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1))
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

Fit_sigmoid = function(x,y,show_plot=FALSE,initial_tau = 100,initial_p=2,initial_C=NULL,remove_n_points = NULL) {
  y[0] = 1
  First_reach = min(which(y[-1]>=1))
  y[1:First_reach]=1

  if (!is.null(remove_n_points)) {
    x = x[-remove_n_points]
    y = y[-remove_n_points]
    
  }
  
  data_temp = data.frame(x = x,y=y)
  data_temp = data_temp[-1,]
  
  if (is.null(initial_C)) {
    initial_C = max(y)
  }
  

  curve.nlslrc = nlsLM(y ~ constant/(1+(x*tau)^p)+1,start = list(tau=1/initial_tau,p=initial_p,constant=initial_C),
                       data = data_temp,control = nls.lm.control(maxiter=1000),na.action=na.exclude,lower = c(0,0,0))

  tau = 1/coef(curve.nlslrc)[1]
  p = coef(curve.nlslrc)[2]
  constant  = coef(curve.nlslrc)[3]
  
  if (show_plot) {
    par(las=1,bty="l")
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1))
    curve(1 + constant*(1/(1+(x/tau)^p)),add=T,lty=2,lwd=2,col="red")
    abline(h=1,lwd=1,lty=2,col='grey')
    legend("topright",legend = c(paste("Tau=",round(tau,digits = 2)),
                                 paste("p=",round(p,digits =2 ))),bty='n')
  }
  
  Normalizing_term = p/(tau)/gamma(1/p)
  Normalised_constant = constant / Normalizing_term # The constant need to be normalized relatively to the area under the curve of the distribution...
  
  r_squared = 1- mean(residuals(curve.nlslrc)^2)/var(y)
  Output_vector = c(tau,p,constant,Normalised_constant,r_squared)
  
  names(Output_vector) = c("tau","p","C","C_normalised","R2")
  return(Output_vector)
}

Fit_sigmoid_robust = function(x,y,show_plot=FALSE,initial_tau = 100) {
  y[0] = 1
  First_reach = min(which(y[-1]>=1))
  y[1:First_reach]=1
  
  data_temp = data.frame(x = x,y=y)
  data_temp = data_temp[-1,]
  
  curve.nlslrc = nlsLM(y ~ constant/(1+(x*tau)^p)+1,start = list(tau=1/initial_tau,p=2,constant=quantile(y,0.9,na.rm=TRUE)),
                       data = data_temp,control = nls.lm.control(maxiter=1000),na.action=na.exclude,lower = c(0,0,0))
  tau = 1/coef(curve.nlslrc)[1]
  p = coef(curve.nlslrc)[2]
  constant  = coef(curve.nlslrc)[3]
  
  if (show_plot) {
    par(las=1,bty="l")
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1))
    curve(1 + constant*(1/(1+(x/tau)^p)),add=T,lty=2,lwd=2,col="red")
    abline(h=1,lwd=1,lty=2,col='grey')
    legend("topright",legend = c(paste("Tau=",round(tau,digits = 2)),
                                 paste("p=",round(p,digits =2 ))),bty='n')
  }
  
  Normalizing_term = p/(tau)/gamma(1/p)
  Normalised_constant = constant / Normalizing_term # The constant need to be normalized relatively to the area under the curve of the distribution...
  
  r_squared = 1- mean(residuals(curve.nlslrc)^2)/var(y)
  Output_vector = c(tau,p,constant,Normalised_constant,r_squared)
  
  names(Output_vector) = c("tau","p","C","C_normalised","R2")
  return(Output_vector)
}


Fit_beta_prime = function(x,y,show_plot=FALSE) {
  y[0] = 1
  First_reach = min(which(y[-1]>=1))
  y[1:First_reach]=1
  
  
  data_temp = data.frame(x = x,y=y)
  data_temp = data_temp[-1,]
  
  curve.nlslrc = nlsLM(y ~ x^(alpha-1)*(x+1)^(-alpha-beta)*constant+1,start = list(alpha=15,beta=2,constant=100),
                       data = data_temp,control = nls.lm.control(maxiter=1000))
  alpha = coef(curve.nlslrc)[1]
  beta = coef(curve.nlslrc)[2]
  constant  = coef(curve.nlslrc)[3]
  
  if (show_plot) {
    par(las=1,bty="l")
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1))
    curve( x^(alpha-1)*(x+1)^(-alpha-beta)*constant+1,add=T,lty=2,lwd=2,col="red")
    abline(h=1,lwd=1,lty=2,col='grey')
    legend("topright",legend = c(paste("Alpha=",round(alpha,digits = 2)),
                                 paste("Beta=",round(beta,digits =2 ))),bty='n')
  }
  
  Normalizing_term = gamma(alpha)*gamma(beta)/gamma(alpha+beta)
  Normalised_constant = constant / Normalizing_term # The constant need to be normalized relatively to the area under the curve of the distribution...
  
  r_squared = 1- mean(residuals(curve.nlslrc)^2)/var(y)
  Output_vector = c(alpha,beta,constant,Normalised_constant,r_squared)
  
  names(Output_vector) = c("alpha","beta","C","C_normalised","R2")
  return(Output_vector)
}

Fit_generalized_beta_prime = function(x,y,show_plot=FALSE) {
  y[0] = 1
  First_reach = min(which(y[-1]>=1))
  y[1:First_reach]=1
  
  y = y[x!=0]
  x = x[x!=0]
  
  data_temp = data.frame(x = x,y=y)
  data_temp = data_temp[-1,]
  
  curve.nlslrc = nlsLM(y ~ 1+constant*(x/tau)^(alpha*p-1)*(1+(x/tau)^p)^(-alpha-beta),
                       start = list(alpha=15,beta=2,tau=1,p=1,constant=max(y)),
                       data = data_temp,control = nls.lm.control(maxiter=1000),lower = c(0,0,0,0,0))
  alpha = coef(curve.nlslrc)[1]
  beta = coef(curve.nlslrc)[2]
  tau = coef(curve.nlslrc)[3]
  p = coef(curve.nlslrc)[4]
  constant  = coef(curve.nlslrc)[5]
  

  
  if (show_plot) {
    par(las=1,bty="l")
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1))
    curve( 1+constant*(x/tau)^(alpha*p-1)*(1+(x/tau)^p)^(-alpha-beta),add=T,lty=2,lwd=2,col="red")
    abline(h=1,lwd=1,lty=2,col='grey')
    legend("topright",legend = c(paste("Alpha=",round(alpha,digits = 2)),
                                 paste("Beta=",round(beta,digits =2 )),
                                 paste("Tau=",round(tau,digits =2 )),
                                 paste("p=",round(p,digits =2 ))),bty='n')
  }
  
  Normalizing_term = p/tau*gamma(alpha)*gamma(beta)/gamma(alpha+beta)
  Normalised_constant = constant / Normalizing_term # The constant need to be normalized relatively to the area under the curve of the distribution...
  
  r_squared = 1- mean(residuals(curve.nlslrc)^2)/var(y)
  Output_vector = c(alpha,beta,tau,p,constant,Normalised_constant,r_squared)
  
  names(Output_vector) = c("alpha","beta","tau","p","C","C_normalised","R2")
  return(Output_vector)
}

Fit_log_cauchy = function(x,y,show_plot=FALSE) {
  y[0] = 1
  First_reach = min(which(y[-1]>=1))
  y[1:First_reach]=1
  
  y = y[x!=0]
  x = x[x!=0]
  
  data_temp = data.frame(x = x,y=y)
  
  
  curve.nlslrc = nlsLM(y ~ sigma/(pi*x)/((log(x)-mu)^2+sigma^2)*constant+1,start = list(mu=1,sigma=5,constant=100),
                       data = data_temp,control = nls.lm.control(maxiter=1000))
  mu = coef(curve.nlslrc)[1]
  sigma = coef(curve.nlslrc)[2]
  constant  = coef(curve.nlslrc)[3]
  
  if (show_plot) {
    par(las=1,bty="l")
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1))
    curve( sigma/(pi*x)/((log(x)-mu)^2+sigma^2)*constant+1,add=T,lty=2,lwd=2,col="red")
    abline(h=1,lwd=1,lty=2,col='grey')
    legend("topright",legend = c(paste("Mu=",round(mu,digits = 2)),
                                 paste("Sigma=",round(sigma,digits =2 ))),bty='n')
  }
  
  Normalizing_term = 1
  Normalised_constant = constant / Normalizing_term # The constant need to be normalized relatively to the area under the curve of the distribution...
  
  r_squared = 1- mean(residuals(curve.nlslrc)^2)/var(y)
  Output_vector = c(mu,sigma,constant,Normalised_constant,r_squared)
  
  names(Output_vector) = c("mu","sigma","C","C_normalised","R2")
  return(Output_vector)
}

Fit_Dagum = function(x,y,show_plot=FALSE) {
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
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1))
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

Fit_log_Laplace= function(x,y,show_plot=FALSE) {
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
    plot(x,y,xaxs='i',yaxs='i',xlab="r",ylab="Pcf(r)",cex.lab=1.3,ylim=c(0,max(y)*1.1))
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

Fit_parametric_pcf_model = function(List_pcf,model="Generalized_gamma") {
  
  N_pcf = length(List_pcf$List_pcf)
  Table_fitted_parameters = c()
  
  if (model=="Power_law") {
    for (k in 1:N_pcf) {
      print(k)
      temp_fit = try(Fit_power_law(x = List_pcf$List_r[[k]],y = List_pcf$List_pcf[[k]],show_plot = FALSE))
      if (class(temp_fit)=="try-error") {
        temp_fit = rep(NA,3)
      }
      Table_fitted_parameters = rbind(Table_fitted_parameters,temp_fit)
    }
    rownames(Table_fitted_parameters) = 1:nrow(Table_fitted_parameters)
  }
  
  
  if (model=="Generalized_gamma") {
    for (k in 1:N_pcf) {
      print(k)
      temp_fit = try(Fit_generalized_gamma(x = List_pcf$List_r[[k]],y = List_pcf$List_pcf[[k]],show_plot = FALSE))
      if (class(temp_fit)=="try-error") {
        temp_fit = rep(NA,6)
      }
      Table_fitted_parameters = rbind(Table_fitted_parameters,temp_fit)
    }
    rownames(Table_fitted_parameters) = 1:nrow(Table_fitted_parameters)
  }
  
  if (model=="Exponential") {
    for (k in 1:N_pcf) {
      print(k)
      temp_fit = try(Fit_exponential(x = List_pcf$List_r[[k]],y = List_pcf$List_pcf[[k]],show_plot = FALSE))
      if (class(temp_fit)=="try-error") {
        temp_fit = rep(NA,4)
      }
      Table_fitted_parameters = rbind(Table_fitted_parameters,temp_fit)
    }
    rownames(Table_fitted_parameters) = 1:nrow(Table_fitted_parameters)
  }
  
  
  
  if (model=="Beta_prime") {
    for (k in 1:N_pcf) {
      temp_fit = try(Fit_beta_prime(x = List_pcf$List_r[[k]],y = List_pcf$List_pcf[[k]],show_plot = FALSE))
      if (class(temp_fit)=="try-error") {
        temp_fit = rep(NA,5)
      }
      Table_fitted_parameters = rbind(Table_fitted_parameters,temp_fit)
    }
    rownames(Table_fitted_parameters) = 1:nrow(Table_fitted_parameters)
  }
  
  if (model=="Log_Cauchy") {
    for (k in 1:N_pcf) {
      temp_fit = Fit_log_cauchy(x = List_pcf$List_r[[k]],y = List_pcf$List_pcf[[k]],show_plot = FALSE)
      Table_fitted_parameters = rbind(Table_fitted_parameters,temp_fit)
    }
    rownames(Table_fitted_parameters) = 1:nrow(Table_fitted_parameters)
  }
  
  
  if (model=="Generalized_beta_prime") {
    for (k in 1:N_pcf) {
      print(k)
      temp_fit = try(Fit_generalized_beta_prime(x = List_pcf$List_r[[k]],y = List_pcf$List_pcf[[k]],show_plot = FALSE))
      if (class(temp_fit)=="try-error") {
        temp_fit = rep(NA,7)
      }
      Table_fitted_parameters = rbind(Table_fitted_parameters,temp_fit)
    }
    rownames(Table_fitted_parameters) = 1:nrow(Table_fitted_parameters)
  }
  
  if (model=="Dagum") {
    for (k in 1:N_pcf) {
      temp_fit = try(Fit_Dagum(x = List_pcf$List_r[[k]],y = List_pcf$List_pcf[[k]],show_plot = FALSE))
      if (class(temp_fit)=="try-error") {
        temp_fit = rep(NA,5)
      }
      Table_fitted_parameters = rbind(Table_fitted_parameters,temp_fit)
    }
    rownames(Table_fitted_parameters) = 1:nrow(Table_fitted_parameters)
  }
  
  
  if (model=="Gamma") {
    for (k in 1:N_pcf) {
      temp_fit = try(Fit_gamma(x = List_pcf$List_r[[k]],y = List_pcf$List_pcf[[k]],show_plot = FALSE))
      if (class(temp_fit)=="try-error") {
        temp_fit = rep(NA,6)
      }
      Table_fitted_parameters = rbind(Table_fitted_parameters,temp_fit)
    }
    rownames(Table_fitted_parameters) = 1:nrow(Table_fitted_parameters)
  }
  
  if (model=="Sigmoid") {
    for (k in 1:N_pcf) {
      print(k)
      temp_fit = try(Fit_sigmoid(x = List_pcf$List_r[[k]],y = List_pcf$List_pcf[[k]],show_plot = FALSE))
      if (class(temp_fit)=="try-error") {
        temp_fit = rep(NA,5)
      }
      Table_fitted_parameters = rbind(Table_fitted_parameters,temp_fit)
    }
    rownames(Table_fitted_parameters) = 1:nrow(Table_fitted_parameters)
  }
  
  if (model=="Sigmoid_robust") {
    for (k in 1:N_pcf) {
      print(k)
      temp_fit = try(Fit_sigmoid(x = List_pcf$List_r[[k]],y = List_pcf$List_pcf[[k]],show_plot = FALSE,remove_n_points = 1))
      if (class(temp_fit)=="try-error") {
        temp_fit = rep(NA,5)
      }
      Table_fitted_parameters = rbind(Table_fitted_parameters,temp_fit)
    }
    rownames(Table_fitted_parameters) = 1:nrow(Table_fitted_parameters)
  }
  
  
  if (model=="Log_Laplace") {
    for (k in 1:N_pcf) {
      temp_fit = Fit_log_Laplace(x = List_pcf$List_r[[k]],y = List_pcf$List_pcf[[k]],show_plot = FALSE)
      Table_fitted_parameters = rbind(Table_fitted_parameters,temp_fit)
    }
    rownames(Table_fitted_parameters) = 1:nrow(Table_fitted_parameters)
  }
  
  
  Table_fitted_parameters = as.data.frame(Table_fitted_parameters)
  return(Table_fitted_parameters)
  
}




Fit_mixed_interaction = function(x,y,show_plot=FALSE) {
  y[0] = 1
  First_reach = min(which(y[-1]>=1))
  y[1:First_reach]=1
  data_temp = data.frame(x = x,y=y)
  data_temp = data_temp[-1,]
  
  curve.nlslrc = nlsLM(y ~ (x^alpha)*exp(-(x*tau)^p)*constant+1,start = list(alpha=-0.5,tau=0.1,p=1,constant=max(y)),
                       data = data_temp,control = nls.lm.control(maxiter=1000),lower = c(-0.9999,0,0,0))
  
  
  
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




