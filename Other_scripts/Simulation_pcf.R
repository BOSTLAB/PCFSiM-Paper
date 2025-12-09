library(spatstat)
source("Desktop/Projet_analyse_pcf/List_function_pcf_analysis.R")


#I) Defining functions to simulate various point patterns


###Thomas process : random cluster center + random number of point in each cluster in a given zone

List_pcf_Thomas = c()

Generate_Thomas_pattern = function(win_size = 10,scale_parameter = 0.15 ,N_simulations = 100) {
  
  
  List_pcf = c()
  List_r = c()
  
  for (k in 1:N_simulations) {
    ppp_temp = rThomas(kappa = 1,scale = scale_parameter,mu = 50,win=square(win_size))
    u = Kest(ppp_temp,r = seq(0,3,length.out=50))
    u = pcf.fv(u,method = "b")
    List_pcf[[k]]=u$pcf
    List_r[[k]] = seq(0,3,length.out=50)
  }
  plot(ppp_temp)
  
  return(List_pcf_Thomas = list(List_r = List_r,List_pcf=List_pcf))
  
}

###
Generate_Cauchy_pattern = function(win_size = 10,scale_parameter = 0.15 ,N_simulations = 100) {
  
  
  List_pcf = c()
  List_r = c()
  
  for (k in 1:N_simulations) {
    ppp_temp = rCauchy(kappa = 1,scale = scale_parameter,mu = 50,win=square(win_size))
    u = Kest(ppp_temp,r = seq(0,3,length.out=50))
    u = pcf.fv(u,method = "b")
    List_pcf[[k]]=u$pcf
    List_r[[k]] = seq(0,3,length.out=50)
  }
  plot(ppp_temp)
  
  return(List_pcf_Thomas = list(List_r = List_r,List_pcf=List_pcf))
  
}



#Same thing but with a percentage of noise 

Generate_Thomas_pattern_noise = function(win_size = 10,scale_parameter = 0.15 ,N_simulations = 100,proportion_noise = 0.2) {
  
  
  List_pcf = c()
  List_r = c()
  
  for (k in 1:N_simulations) {
    ppp_temp = rThomas(kappa = 1,scale = scale_parameter,mu = 50,win=square(win_size))
    N_points = ppp_temp$n
    ppp_noise = rpoispp(lambda = proportion_noise*N_points/area(square(win_size)),win = square(win_size))
    ppp_temp = superimpose(ppp_temp,ppp_noise)
    u = Kest(ppp_temp,r = seq(0,3,length.out=50))
    u = pcf.fv(u,method = "b")
    List_pcf[[k]]=u$pcf
    List_r[[k]] = seq(0,3,length.out=50)
  }
  plot(ppp_temp)
  
  return(List_pcf_Thomas = list(List_r = List_r,List_pcf=List_pcf))
  
}


### Mimic of epithelial layer (gut): sinusoidal

Generate_sinus_patterns = function(N_cells = 1000,tau=2,noise = 0.25,win_size = 10,scaling_factor =2,N_simulations = 50) {
  
  List_pcf = c()
  List_r = c()
  
  for (k in 1:N_simulations) {
    
    x = runif(n = N_cells,min = 0,max = win_size)
    y = scaling_factor*sin(x*tau)+runif(n = N_cells,min = -noise*scaling_factor,max = noise*scaling_factor)
    y = y+win_size/2
    ppp_temp = ppp(x,y,window = square(win_size))
    ppp_temp = rotate(ppp_temp,angle = runif(n = 1,min = -pi,max = pi),centre = c(mean(x),mean(y)))
    x_temp = ppp_temp$x
    y_temp = ppp_temp$y
    ppp_temp = ppp(x_temp,y_temp,window = square(win_size))
    
    u = Kest(ppp_temp,r = seq(0,3,length.out=200))
    u = pcf.fv(u,method = "c")
    List_pcf[[k]]=u$pcf
    List_r[[k]] = seq(0,3,length.out=200)
  }
  plot(ppp_temp,pch=21,bg='grey',lwd=0.5)
  return(List_pcf_ring = list(List_r = List_r,List_pcf=List_pcf))
}

### Ring :

Generate_ring_patterns = function(N_cells = 1000,win_size = 10,width_factor =0.2,radius = 2,N_simulations = 50) {
  
  List_pcf= c()
  List_r= c()
  
  for (k in 1:N_simulations) {
    
    
    random_theta = runif(n = N_cells,min = -pi,max = pi)
    random_radius = radius* (1+runif(N_cells,min = -width_factor,max = width_factor))
    
    x = random_radius * cos(random_theta)
    y = random_radius * sin(random_theta)
    
    x = x + win_size/2
    y = y + win_size/2
    
    ppp_temp = ppp(x,y,window = square(win_size))
    u = Kest(ppp_temp,r = seq(0,3,length.out=200))
    u = pcf.fv(u,method = "c")
    List_pcf[[k]]=u$pcf
    List_r[[k]]=u$r
  }
  plot(ppp_temp,pch=21,bg="grey",lwd=0.5)
  return(List_pcf_ring = list(List_r = List_r,List_pcf=List_pcf))
}

Plot_average_pcf = function(List_pcf) {
  
  Matrix_pcf = matrix(unlist(List_pcf$List_pcf),nrow = length(List_pcf$List_pcf[[1]]),byrow = FALSE)
  Mean_pcf_value = rowMeans(Matrix_pcf)
  Quantile_top = apply(Matrix_pcf,MARGIN = 1,FUN = function(x) {quantile(x,0.99)})
  Quantile_bottom= apply(Matrix_pcf,MARGIN = 1,FUN = function(x) {quantile(x,0.01)})
  
  par(las=1,bty='l')
  plot(NULL,xlim=c(0,max(List_pcf$List_r[[1]])),ylim=c(0,max(Quantile_top)*1.2),xlab="r",ylab="pcf(r)",xaxs='i',yaxs='i',cex.lab=1.3)
  abline(h=1,lwd=2,lty=2,col="grey")
  points(List_pcf$List_r[[1]],Mean_pcf_value,type="l",lwd=2,col="black")
  points(List_pcf$List_r[[1]],Quantile_top,type="l",lwd=1,lty=2,col="grey")
  points(List_pcf$List_r[[1]],Quantile_bottom,type="l",lwd=1,lty=2,col="grey")
  
}

#II)Performing simulations to check the impact of the parameters 

#A)Start with Thomas Process with various scale values and fit with an exponential function

List_scale_parameter = seq(0.05,0.25,length.out=20) 
Table_mean_values_Thomas_fit = c()
Table_sd_values_Thomas_fit = c()

for (k in 1:length(List_scale_parameter)) {
  Thomas_temp = Generate_Thomas_pattern(win_size = 10,N_simulations = 50,scale_parameter = List_scale_parameter[k])
  Fit_temp = Fit_parametric_pcf_model(Thomas_temp,model = "Exponential")
  Table_mean_values_Thomas_fit = rbind(Table_mean_values_Thomas_fit,apply(Fit_temp,MARGIN = 2,FUN = mean))
  Table_sd_values_Thomas_fit = rbind(Table_sd_values_Thomas_fit,apply(Fit_temp,MARGIN = 2,FUN = sd))
  
}

par(las=1,bty='l')
plot(List_scale_parameter,Table_mean_values_Thomas_fit[,1],xlab="Scale parameter",ylab="Estimated tau parameter",
     cex=1.5,pch=21,bg="red3")
segments(x0 = List_scale_parameter,x1 = List_scale_parameter,
         y0 = Table_mean_values_Thomas_fit[,1]+Table_sd_values_Thomas_fit[,1],
         y1 =Table_mean_values_Thomas_fit[,1]-Table_sd_values_Thomas_fit[,1],)

m = lm(Table_mean_values_Thomas_fit[,1]~List_scale_parameter)
abline(coef(m),lty=2,col='grey')
legend("topleft",bty = "n",legend = paste("R2=",round(summary(m)$r.squared,digits = 3)))

#B)Same process but with a sigmoid curve

List_scale_parameter = seq(0.05,0.25,length.out=20) 
Table_mean_values_Thomas_fit = c()
Table_sd_values_Thomas_fit = c()

for (k in 1:length(List_scale_parameter)) {
  Thomas_temp = Generate_Thomas_pattern(win_size = 10,N_simulations = 50,scale_parameter = List_scale_parameter[k])
  Fit_temp = Fit_parametric_pcf_model(Thomas_temp,model = "Sigmoid")
  Table_mean_values_Thomas_fit = rbind(Table_mean_values_Thomas_fit,apply(Fit_temp,MARGIN = 2,FUN = mean))
  Table_sd_values_Thomas_fit = rbind(Table_sd_values_Thomas_fit,apply(Fit_temp,MARGIN = 2,FUN = sd))
  
}

par(las=1,bty='l')
plot(List_scale_parameter,Table_mean_values_Thomas_fit[,1],xlab="Scale parameter",ylab="Estimated tau parameter",
     cex=1.5,pch=21,bg="red3")
segments(x0 = List_scale_parameter,x1 = List_scale_parameter,
         y0 = Table_mean_values_Thomas_fit[,1]+Table_sd_values_Thomas_fit[,1],
         y1 =Table_mean_values_Thomas_fit[,1]-Table_sd_values_Thomas_fit[,1],)
m = lm(Table_mean_values_Thomas_fit[,1]~List_scale_parameter)
abline(coef(m),lty=2,col='grey')
legend("topleft",bty = "n",legend = paste("R2=",round(summary(m)$r.squared,digits = 2)))

### Both methods are efficient for simple structures. What about more complex shapes ??

#C) Sigmoid fit on a ring process

List_scale_parameter = seq(0.1,0.3,length.out=20) 
Table_mean_values_ring_fit = c()
Table_sd_values_ring_fit = c()

for (k in 1:length(List_scale_parameter)) {
  Ring_temp = Generate_ring_patterns(win_size = 10,N_simulations = 50,width_factor  = List_scale_parameter[k])
  Fit_temp = Fit_parametric_pcf_model(Ring_temp,model = "Sigmoid")
  Table_mean_values_ring_fit = rbind(Table_mean_values_ring_fit,apply(Fit_temp,MARGIN = 2,FUN = mean))
  Table_sd_values_ring_fit = rbind(Table_sd_values_ring_fit,apply(Fit_temp,MARGIN = 2,FUN = sd))
}

par(las=1,bty='l')
plot(List_scale_parameter,Table_mean_values_ring_fit[,1],xlab="Scale parameter",ylab="Estimated tau parameter",
     cex=1.5,pch=21,bg="red3",cex.lab=1.3)
segments(x0 = List_scale_parameter,x1 = List_scale_parameter,
         y0 = Table_mean_values_ring_fit[,1]+Table_sd_values_ring_fit[,1],
         y1 =Table_mean_values_ring_fit[,1]-Table_sd_values_ring_fit[,1],)

m = lm(Table_mean_values_ring_fit[,1]~List_scale_parameter)
abline(coef(m),lty=2,col='grey')
legend("topleft",bty = "n",legend = paste("R2=",round(summary(m)$r.squared,digits = 2)))
#D) Thomas process with an increasing fraction of noise, i.e. Homogenous Poisson Point Process


List_noise_parameter = seq(0,2,length.out=20) 
Table_mean_values_Thomas_fit = c()
Table_sd_values_Thomas_fit = c()

for (k in 1:length(List_noise_parameter)) {
  Thomas_temp = Generate_Thomas_pattern_noise(win_size = 10,N_simulations = 50,scale_parameter = 0.2,proportion_noise = List_noise_parameter[k])
  Fit_temp = Fit_parametric_pcf_model(Thomas_temp,model = "Exponential")
  Table_mean_values_Thomas_fit = rbind(Table_mean_values_Thomas_fit,apply(Fit_temp,MARGIN = 2,FUN = mean))
  Table_sd_values_Thomas_fit = rbind(Table_sd_values_Thomas_fit,apply(Fit_temp,MARGIN = 2,FUN = sd))
  
}

par(las=1,bty='l')
plot(List_noise_parameter,Table_mean_values_Thomas_fit[,3],xlab="Noise parameter",ylab="Estimated C normalised parameter",
     cex=1.5,pch=21,bg="red3",cex.lab=1.3,log="y")
segments(x0 = List_noise_parameter,x1 = List_noise_parameter,
         y0 = Table_mean_values_Thomas_fit[,3]+Table_sd_values_Thomas_fit[,3],
         y1 =Table_mean_values_Thomas_fit[,3]-Table_sd_values_Thomas_fit[,3],)
m = lm(log10(Table_mean_values_Thomas_fit[,3])~List_noise_parameter)
abline(coef(m),lty=2,col='grey')
legend("topright",bty = "n",legend = paste("R2=",round(summary(m)$r.squared,digits = 2)))


#E) Same thing with a sigmoid curve

List_noise_parameter = seq(0,2,length.out=20) 
Table_mean_values_Thomas_fit = c()
Table_sd_values_Thomas_fit = c()

for (k in 1:length(List_noise_parameter)) {
  Thomas_temp = Generate_Thomas_pattern_noise(win_size = 10,N_simulations = 50,scale_parameter = 0.2,proportion_noise = List_noise_parameter[k])
  Fit_temp = Fit_parametric_pcf_model(Thomas_temp,model = "Sigmoid")
  Table_mean_values_Thomas_fit = rbind(Table_mean_values_Thomas_fit,apply(Fit_temp,MARGIN = 2,FUN = mean))
  Table_sd_values_Thomas_fit = rbind(Table_sd_values_Thomas_fit,apply(Fit_temp,MARGIN = 2,FUN = sd))
  
}

par(las=1,bty='l')
plot(List_noise_parameter,Table_mean_values_Thomas_fit[,4],xlab="Scale parameter",ylab="Estimated C normalised parameter",
     cex=1.5,pch=21,bg="red3",cex.lab=1.3,log="y",xlim=c(0,1.05))
segments(x0 = List_noise_parameter,x1 = List_noise_parameter,
         y0 = Table_mean_values_Thomas_fit[,4]+Table_sd_values_Thomas_fit[,4],
         y1 =Table_mean_values_Thomas_fit[,4]-Table_sd_values_Thomas_fit[,4],)
m = lm(log10(Table_mean_values_Thomas_fit[,4])~List_noise_parameter)
abline(coef(m),lty=2,col='grey')
legend("topright",bty = "n",legend = paste("R2=",round(summary(m)$r.squared,digits = 2)))
dev.off()

plot(List_noise_parameter,Table_mean_values_Thomas_fit[,1])
m_tau = lm((Table_mean_values_Thomas_fit[,1])~List_noise_parameter)

#F) Link with the number of cells 



List_window_size =seq(5,30,length.out=20)
Table_mean_values_Thomas_fit = c()
Table_sd_values_Thomas_fit = c()

for (k in 1:length(List_window_size)) {
  Thomas_temp = Generate_Thomas_pattern(win_size = List_window_size[k],N_simulations = 50,scale_parameter = 0.2)
  Fit_temp = Fit_parametric_pcf_model(Thomas_temp,model = "Sigmoid")
  Table_mean_values_Thomas_fit = rbind(Table_mean_values_Thomas_fit,apply(Fit_temp,MARGIN = 2,FUN = mean))
  Table_sd_values_Thomas_fit = rbind(Table_sd_values_Thomas_fit,apply(Fit_temp,MARGIN = 2,FUN = sd))
  
}

    useDingbats = FALSE,width = 6.5,height = 6.5)
plot(List_window_size^2,Table_mean_values_Thomas_fit[,4],xlab="FoV area",ylab="Estimated C normalised parameter",
     cex=1.5,pch=21,bg="red3",cex.lab=1.3,log="y",xlim=c(0,1000),xaxs='i',ylim=c(0.45,0.65))
segments(x0 = List_window_size^2,x1 = List_window_size^2,
         y0 = Table_mean_values_Thomas_fit[,4]+Table_sd_values_Thomas_fit[,4],
         y1 =Table_mean_values_Thomas_fit[,4]-Table_sd_values_Thomas_fit[,4],)
y = List_window_size^2
m = lm(log10(Table_mean_values_Thomas_fit[,4])~y)
abline(coef(m),lty=2,col='grey')
legend("topright",bty = "n",legend = paste("R2=",round(summary(m)$r.squared,digits = 2)))

plot(List_window_size^2,Table_mean_values_Thomas_fit[,1],xlab="FoV area",ylab="Estimated tau parameter",
     cex=1.5,pch=21,bg="red3",cex.lab=1.3,log="y",xlim=c(0,1000),xaxs='i',ylim=c(0.3,0.4))
segments(x0 = List_window_size^2,x1 = List_window_size^2,
         y0 = Table_mean_values_Thomas_fit[,1]+Table_sd_values_Thomas_fit[,1],
         y1 =Table_mean_values_Thomas_fit[,1]-Table_sd_values_Thomas_fit[,1],)
y = List_window_size^2
m = lm(log10(Table_mean_values_Thomas_fit[,1])~y)
abline(coef(m),lty=2,col='grey')
legend("topright",bty = "n",legend = paste("R2=",round(summary(m)$r.squared,digits = 2)))


plot(List_window_size^2,Table_mean_values_Thomas_fit[,2],xlab="FoV area",ylab="Estimated p parameter",
     cex=1.5,pch=21,bg="red3",cex.lab=1.3,log="y",xlim=c(0,1000),xaxs='i',ylim=c(2.5,4.5))
segments(x0 = List_window_size^2,x1 = List_window_size^2,
         y0 = Table_mean_values_Thomas_fit[,2]+Table_sd_values_Thomas_fit[,2],
         y1 =Table_mean_values_Thomas_fit[,2]-Table_sd_values_Thomas_fit[,2],)
y = List_window_size^2
m = lm(log10(Table_mean_values_Thomas_fit[,2])~y)
abline(coef(m),lty=2,col='grey')
legend("topright",bty = "n",legend = paste("R2=",round(summary(m)$r.squared,digits = 2)))


### Plot for visualisation

ppp_temp = rThomas(kappa = 1,scale = 0.05,mu = 50,win=square(10))
plot(ppp_temp,pch=21,bg="grey",lwd=0.5,main="Scale = 0.05")

ppp_temp = rThomas(kappa = 1,scale = 0.15,mu = 50,win=square(10))
plot(ppp_temp,pch=21,bg="grey",lwd=0.5,main="Scale = 0.15")

ppp_temp = rThomas(kappa = 1,scale = 0.25,mu = 50,win=square(10))
plot(ppp_temp,pch=21,bg="grey",lwd=0.5,main="Scale = 0.25")



#II) Looking at the p parameter 

#A) Comparing Neymann-Scott point process (Cauchy kernel) with Thomas process

List_scale_parameter = seq(0.05,0.25,length.out=20) 
Table_mean_values_Thomas_fit = c()
Table_sd_values_Thomas_fit = c()

for (k in 1:length(List_scale_parameter)) {
  Thomas_temp = Generate_Thomas_pattern(win_size = 10,N_simulations = 50,scale_parameter = List_scale_parameter[k])
  Fit_temp = Fit_parametric_pcf_model(Thomas_temp,model = "Sigmoid")
  Table_mean_values_Thomas_fit = rbind(Table_mean_values_Thomas_fit,apply(Fit_temp,MARGIN = 2,FUN = mean))
  Table_sd_values_Thomas_fit = rbind(Table_sd_values_Thomas_fit,apply(Fit_temp,MARGIN = 2,FUN = sd))
}

Table_mean_values_Cauchy_fit = c()
Table_sd_values_Cauchy_fit = c()

for (k in 1:length(List_scale_parameter)) {
  Cauchy_temp = Generate_Cauchy_pattern(win_size = 10,N_simulations = 50,scale_parameter = List_scale_parameter[k])
  Fit_temp = Fit_parametric_pcf_model(Cauchy_temp,model = "Sigmoid")
  Table_mean_values_Cauchy_fit = rbind(Table_mean_values_Cauchy_fit,apply(Fit_temp,MARGIN = 2,FUN = mean))
  Table_sd_values_Cauchy_fit = rbind(Table_sd_values_Cauchy_fit,apply(Fit_temp,MARGIN = 2,FUN = sd))
}

plot(Table_mean_values_Cauchy_fit[,2],Table_mean_values_Thomas_fit[,2])

par(las=1,bty='l')
test_temp = t.test(Table_mean_values_Thomas_fit[,2],Table_mean_values_Cauchy_fit[,2],paired = TRUE)
boxplot(Table_mean_values_Thomas_fit[,2],Table_mean_values_Cauchy_fit[,2],outline=F,
        names=c("Thomas","Cauchy"),ylab="p parameter value",main=paste("p = ",round(test_temp$p.value,100)))


#III)Plotting function



ppp_temp = rThomas(kappa = 1,scale = 0.05,mu = 50,win=square(10))
plot(ppp_temp,pch=21,bg="grey",lwd=0.5,main="Scale = 0.05")

ppp_temp = rThomas(kappa = 1,scale = 0.15,mu = 50,win=square(10))
plot(ppp_temp,pch=21,bg="grey",lwd=0.5,main="Scale = 0.15")

ppp_temp = rThomas(kappa = 1,scale = 0.25,mu = 50,win=square(10))
plot(ppp_temp,pch=21,bg="grey",lwd=0.5,main="Scale = 0.25")


Generate_ring_patterns(win_size = 10,N_simulations = 50,width_factor  = 0.1)
Generate_ring_patterns(win_size = 10,N_simulations = 50,width_factor  = 0.3)




