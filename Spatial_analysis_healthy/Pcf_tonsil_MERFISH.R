
#I)Loading function/object
library(SingleCellExperiment)
library(RColorBrewer)
library(balagan)

CCC_computation = function(x,y) {
  return(2*sd(x,na.rm = T)*sd(y,na.rm = T)/(var(x,na.rm = T)+var(y,na.rm = T)+(mean(x)-mean(y))^2))
}

#A)List functions

source("Desktop/Projet_analyse_pcf/List_function_pcf_analysis.R")

Raw_metadata = read.delim("Annotation_human_tonsil_MERFISH.csv",sep=",")


sce = SingleCellExperiment(assays = list(Raw_intensity = matrix(0,nrow = 2,ncol =nrow(Raw_metadata) )), 
                           metadata = list(dimension = "2D", 
                                           N_core = 8, Is_nuc_cyt = F))
colLabels(sce) = as.numeric(Raw_metadata$Clustering)
sce$Location_Center_X = Raw_metadata$cell_centroid_x
sce$Location_Center_Y = Raw_metadata$cell_centroid_y
sce$ImageNumber =  1



#II)Computation itself 

#A)Computing pcf
List_pcf = Compute_pcf(sce,r_vector = seq(0,2000,length.out=50),computation_method = "derivative",verbose = TRUE)


#B)Checking the functions on an example

k = 19
x = List_pcf$List_r[[k]]
y = List_pcf$List_pcf[[k]]

m_GG = try(Fit_generalized_gamma(x,y,show_plot = TRUE),silent = TRUE)
m_G = try(Fit_gamma(x,y,show_plot = TRUE))
m_expo = try(Fit_exponential(x,y,show_plot = TRUE))
m_sigmoid = try(Fit_sigmoid(x,y,show_plot = TRUE))
m_power = try(Fit_power_law(x,y,show_plot = TRUE))
m_bp = try(Fit_beta_prime(x,y,show_plot = TRUE))

par(las=1,bty='l')
plot(x,y,pch=21,bg="grey70",xlim=c(0,2050),ylim=c(0,12.5),xaxs='i',yaxs='i',xlab="r (µm)",ylab="Pcf(r)",cex.lab=1.3,cex=1.2)
abline(h=1,lty=2,col='black')
curve(1+ (x/m_power[2])^-m_power[1],from=0,to=2020,add=T,col=brewer.pal(6, "Spectral")[1])
curve(1+m_expo[2]*exp(-(x/m_expo[1])),from=0,to=2020,add=T,col=brewer.pal(6, "Spectral")[2])
curve(1+m_G[3]*exp(-(x/m_G[2]))*x^(m_G[1]),from=0,to=2020,add=T,col=brewer.pal(6, "Spectral")[3])
curve(1+m_GG[4]*exp(-(x/m_GG[2])^m_GG[3])*x^(m_GG[1]),from=0,to=2020,add=T,col=brewer.pal(6, "Spectral")[4])
curve(1+m_sigmoid[3]/(1+(x/m_sigmoid[1])^m_sigmoid[2]),from=0,to=2020,add=T,col=brewer.pal(6, "Spectral")[5])
curve( x^(m_bp[1]-1)*(x+1)^(-m_bp[1]-m_bp[2])*m_bp[3]+1,from=0,to=2020,add=T,col=brewer.pal(6, "Spectral")[6])
#C)Checking the functions on an example

Results_generalised_gamma = Fit_parametric_pcf_model(List_pcf,model = "Generalized_gamma")
Results_gamma= Fit_parametric_pcf_model(List_pcf,model = "Gamma")
Results_exponential = Fit_parametric_pcf_model(List_pcf,model = "Exponential")
Results_power_law = Fit_parametric_pcf_model(List_pcf,model = "Power_law")
Results_sigmoid= Fit_parametric_pcf_model(List_pcf,model = "Sigmoid")
Results_beta_prime= Fit_parametric_pcf_model(List_pcf,model = "Beta_prime")


par(las=1,bty='l')
boxplot(Results_power_law$R2,Results_exponential$R2,Results_gamma$R2,
        Results_generalised_gamma$R2,Results_sigmoid$R2,Results_beta_prime$R2,ylim=c(0.4,1),
        xaxs='i',yaxs='i',names = c("Power law","Exponential","Gamma","GG","Sigmoid","Beta prime"),
        xlab="Model",ylab="R2",cex.lab=1.3,col=brewer.pal(6, "Spectral"),cex.axis=0.7)



selected_clustering = 5
Plot_selected_cluster = function(selected_cluster,title_show=NULL) {
  ppp_temp = ppp(x = Raw_metadata$cell_centroid_x[Raw_metadata$Clustering==selected_cluster],
                 y = Raw_metadata$cell_centroid_y[Raw_metadata$Clustering==selected_cluster],
                 window = owin(xrange = range(Raw_metadata$cell_centroid_x),yrange = range(Raw_metadata$cell_centroid_y)))
  plot(ppp_temp,cex=0.2,main=title_show)
  
}

Plot_selected_cluster(selected_cluster = 22,title_show = "Cluster 22, p = 2.71")
Plot_selected_cluster(selected_cluster = 1,title_show = "Cluster 1, p = 1.03")
Plot_selected_cluster(selected_cluster = 5,title_show = "Cluster 5, p = 0.67")


Plot_selected_cluster(selected_cluster = 4,title_show = "Cluster 4, tau = 355µm")
Plot_selected_cluster(selected_cluster = 8,title_show = "Cluster 8, tau = 157µm")
Plot_selected_cluster(selected_cluster = 16,title_show = "Cluster 16, tau = 73µm")

Plot_selected_cluster(selected_cluster = 12,title_show = "Cluster 12, C = 257")
Plot_selected_cluster(selected_cluster = 23,title_show = "Cluster 23, C = 1164")
Plot_selected_cluster(selected_cluster = 18,title_show = "Cluster 18, C = 8633")

#III) Sampling robustness

#A)Creating a function to simulate 'small FoV'

Perform_spatial_sampling_fitting = function(sce,N_min_cells = 1000,FoV_width=500,N_sampling = 500) {
  List_sampled_cells = c()
  cat("Perfoming sampling....")
  for (k in 1:N_sampling) {
    print(k)
    x = Random_spatial_sampling(sce,width_FOV = FoV_width,height_FOV = FoV_width,N_samplings = 1,plot_result = FALSE)
    List_sampled_cells[[k]]= x$List_sampled_cells
  }
  cat("done ! \n")
  
  N_cells_sampling = unlist(lapply(List_sampled_cells, length))
  List_sampled_cells_filtered = List_sampled_cells[N_cells_sampling>=N_min_cells] # Some empty regions -> removing them
  
  cat("Perfoming fitting...")
  
  
  Table_sampling_tau = c()
  Table_sampling_p = c()
  Table_sampling_C_normalised = c()
  Table_sampling_R2 = c()
  
  for (k in 1:length(List_sampled_cells_filtered)) {
    print(k)
    sce_temp = sce[,List_sampled_cells_filtered[[k]]]
    List_pcf_temp = Compute_pcf(sce_temp,r_vector = seq(0,FoV_width/sqrt(2),length.out=50),computation_method = "derivative")
    Results_sigmoid_temp= Fit_parametric_pcf_model(List_pcf_temp,model = "Sigmoid")
    rownames(Results_sigmoid_temp) = List_pcf_temp$Annotation$Cluster
    Results_sigmoid_temp = Results_sigmoid_temp[as.character(List_pcf$Annotation$Cluster),]
    Table_sampling_tau = rbind(Table_sampling_tau,Results_sigmoid_temp$tau)
    Table_sampling_p = rbind(Table_sampling_p,Results_sigmoid_temp$p)
    Table_sampling_C_normalised = rbind(Table_sampling_C_normalised,Results_sigmoid_temp$C_normalised)
    Table_sampling_R2 = rbind(Table_sampling_R2,Results_sigmoid_temp$R2)
  }
  cat("done ! \n")
  return(list(Tau = Table_sampling_tau,P = Table_sampling_p,C_normalised = Table_sampling_C_normalised,R2 = Table_sampling_R2))
  
}

#Typical size of TMA: 600µm radius core ~ 500µm width square,

#B)Performing the simulations
Sampling_500µm = Perform_spatial_sampling_fitting(sce,FoV_width = 500,N_sampling = 500)
Sampling_1000µm = Perform_spatial_sampling_fitting(sce,FoV_width = 1000,N_sampling = 500)
Sampling_1500µm = Perform_spatial_sampling_fitting(sce,FoV_width = 1500,N_sampling = 500)
Sampling_2000µm = Perform_spatial_sampling_fitting(sce,FoV_width = 2000,N_sampling = 500)


par(las=1,bty='l',mfrow=c(4,1))
Selected_cluster = 16
hist(Sampling_500µm[[1]][,Selected_cluster],50,xlim=c(0,500),main="500µm FoV",xlab="Estimated tau parameter")
abline(v=Results_sigmoid$tau[Selected_cluster],lty=2)
abline(v=mean(Sampling_500µm[[1]][,Selected_cluster],na.rm=T),lty=2,col='red')

hist(Sampling_1000µm[[1]][,Selected_cluster],50,xlim=c(0,500),main="1000µm FoV",xlab="Estimated tau parameter")
abline(v=Results_sigmoid$tau[Selected_cluster],lty=2)
abline(v=mean(Sampling_1000µm[[1]][,Selected_cluster],na.rm=T),lty=2,col='red',xlab="Estimated tau parameter")

hist(Sampling_1500µm[[1]][,Selected_cluster],50,xlim=c(0,500),main="1500µm FoV",xlab="Estimated tau parameter")
abline(v=Results_sigmoid$tau[Selected_cluster],lty=2)
abline(v=mean(Sampling_1500µm[[1]][,Selected_cluster],na.rm=T),lty=2,col='red')

hist(Sampling_2000µm[[1]][,Selected_cluster],50,xlim=c(0,500),main="2000µm FoV",xlab="Estimated tau parameter")
abline(v=Results_sigmoid$tau[Selected_cluster],lty=2)
abline(v=mean(Sampling_2000µm[[1]][,Selected_cluster],na.rm=T),lty=2,col='red')




#C)What is the quality of the inference ?

par(las=1,bty='l')
boxplot(colMeans(Sampling_500µm[[4]],na.rm = T),
        colMeans(Sampling_1000µm[[4]],na.rm = T),
        colMeans(Sampling_1500µm[[4]],na.rm = T),
        colMeans(Sampling_2000µm[[4]],na.rm = T),
        Results_sigmoid$R2,outline=F,
        names=c("500µm","1000µm","1500µm","2000µm","Initial data"),ylim=c(0,1),yaxs='i',ylab="Fitting R2",cex.lab=1.3)


boxplot(colMeans(is.na(Sampling_500µm[[1]]),na.rm = T)*100,
        colMeans(is.na(Sampling_1000µm[[1]]),na.rm = T)*100,
        colMeans(is.na(Sampling_1500µm[[1]]),na.rm = T)*100,
        colMeans(is.na(Sampling_2000µm[[1]]),na.rm = T)*100,
        names=c("500µm","1000µm","1500µm","2000µm"),ylim=c(0,100),yaxs='i',ylab="Proportion of absence of fit (%)",cex.lab=1.3,outline=F)


Compute_QC_sampling = function(x) {
  Sum_NA = colSums(is.na(x))
  Sum_low_quality = colSums(x<0.6,na.rm = T)
  Sum_pass = colSums(x>=0.6 & !is.na(x),na.rm = T)
  X = rbind(Sum_pass,Sum_low_quality,Sum_NA)
  X = t(t(X)/colSums(X))
  return(X)
}
barplot(Compute_QC_sampling(Sampling_1500µm[[4]])*100,col=c("darkseagreen","lightsalmon","firebrick"),xlab="Clusters",
        ylab="Proportion of fittings (%)",cex.lab=1.3)

plot(Results_sigmoid$tau,colMedians(Sampling_500µm[[4]],na.rm = T))
plot(Results_sigmoid$tau,colMeans(is.na(Sampling_500µm[[4]]),na.rm = T))

#D)Comparing the average estimations and the real values

#i) For tau parameter

List_coef = c()
List_coef_sd = c()

par(las=1,bty='l')
plot(Results_sigmoid$tau,colMeans(Sampling_500µm[[1]],na.rm = T),xlim=c(0,400),ylim=c(0,400),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="Tau parameter values (µm)",ylab="Mean estimated tau values (µm)",main="500µm FoV")
abline(0,1,lwd=2,lty=2,col="grey")
m = lm(colMeans(Sampling_500µm[[1]],na.rm = T)~0+Results_sigmoid$tau)
abline(0,coef(m),lwd=2,lty=2,col="black")
List_coef = c(List_coef,coef(m))
List_coef_sd = c(List_coef_sd,summary(m)$coefficients[2])

plot(Results_sigmoid$tau,colMeans(Sampling_1000µm[[1]],na.rm = T),xlim=c(0,400),ylim=c(0,400),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="Tau parameter values (µm)",ylab="Mean estimated tau values (µm)",main="1000µm FoV")
abline(0,1,lwd=2,lty=2,col="grey")
m = lm(colMeans(Sampling_1000µm[[1]],na.rm = T)~0+Results_sigmoid$tau)
abline(0,coef(m),lwd=2,lty=2,col="black")
List_coef = c(List_coef,coef(m))
List_coef_sd = c(List_coef_sd,summary(m)$coefficients[2])

plot(Results_sigmoid$tau,colMeans(Sampling_1500µm[[1]],na.rm = T),xlim=c(0,400),ylim=c(0,400),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="Tau parameter values (µm)",ylab="Mean estimated tau values (µm)",main="1500µm FoV")
abline(0,1,lwd=2,lty=2,col="grey")
m = lm(colMeans(Sampling_1500µm[[1]],na.rm = T)~0+Results_sigmoid$tau)
abline(0,coef(m),lwd=2,lty=2,col="black")
List_coef = c(List_coef,coef(m))
List_coef_sd = c(List_coef_sd,summary(m)$coefficients[2])

plot(Results_sigmoid$tau,colMeans(Sampling_2000µm[[1]],na.rm = T),xlim=c(0,400),ylim=c(0,400),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="Tau parameter values (µm)",ylab="Mean estimated tau values (µm)",main="2000µm FoV")
abline(0,1,lwd=2,lty=2,col="grey")
m = lm(colMeans(Sampling_2000µm[[1]],na.rm = T)~0+Results_sigmoid$tau)
abline(0,coef(m),lwd=2,lty=2,col="black")
List_coef = c(List_coef,coef(m))
List_coef_sd = c(List_coef_sd,summary(m)$coefficients[2])

par(las=1,bty='l')
x =barplot(List_coef,ylim=c(0,1),names.arg = c("500µm","1000µm","1500µm","2000µm"),xlab="Fov size",ylab="Estimated/real tau value slope",cex.lab=1.3)
segments(x0 = x,x1 =x,y0 = List_coef-List_coef_sd,y1 = List_coef+List_coef_sd)



par(las=1,bty='l')
barplot(c(CCC_computation(Results_sigmoid$tau,colMeans(Sampling_500µm[[1]],na.rm = T)),
        CCC_computation(Results_sigmoid$tau,colMeans(Sampling_1000µm[[1]],na.rm = T)),
        CCC_computation(Results_sigmoid$tau,colMeans(Sampling_1500µm[[1]],na.rm = T)),
        CCC_computation(Results_sigmoid$tau,colMeans(Sampling_2000µm[[1]],na.rm = T))),
        ylim=c(0,1),ylab="CCC coefficient",cex.lab=1.3,names.arg = c("500µm","1000µm","1500µm","2000µm"))

barplot(c(cor(Results_sigmoid$tau,colMeans(Sampling_500µm[[1]],na.rm = T)),
          cor(Results_sigmoid$tau,colMeans(Sampling_1000µm[[1]],na.rm = T)),
          cor(Results_sigmoid$tau,colMeans(Sampling_1500µm[[1]],na.rm = T)),
          cor(Results_sigmoid$tau,colMeans(Sampling_2000µm[[1]],na.rm = T))),
        ylim=c(0,1),ylab="Pearson's R",cex.lab=1.3,names.arg = c("500µm","1000µm","1500µm","2000µm"))

barplot(List_coef,ylim=c(0,1),ylab="Slope coefficient",cex.lab=1.3,names.arg = c("500µm","1000µm","1500µm","2000µm"))

#ii) For p parameter

List_coef = c()
List_coef_sd = c()

par(las=1,bty='l')

plot(Results_sigmoid$p,colMeans(Sampling_500µm[[2]],na.rm = T),xlim=c(0,10),ylim=c(0,10),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="p parameter values (µm)",ylab="Mean estimated p values (µm)",main="500µm FoV")
abline(0,1,lwd=2,lty=2,col="grey")
m = lm(colMeans(Sampling_500µm[[2]],na.rm = T)~0+Results_sigmoid$p)
abline(0,coef(m),lwd=2,lty=2,col="black")
List_coef = c(List_coef,coef(m))
List_coef_sd = c(List_coef_sd,summary(m)$coefficients[2])

plot(Results_sigmoid$p,colMeans(Sampling_1000µm[[2]],na.rm = T),xlim=c(0,10),ylim=c(0,10),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="p parameter values (µm)",ylab="Mean estimated p values (µm)",main="1000µm FoV")
abline(0,1,lwd=2,lty=2,col='grey')
m = lm(colMeans(Sampling_1000µm[[2]],na.rm = T)~0+Results_sigmoid$p)
abline(0,coef(m),lwd=2,lty=2,col="black")
List_coef = c(List_coef,coef(m))
List_coef_sd = c(List_coef_sd,summary(m)$coefficients[2])

plot(Results_sigmoid$p,colMeans(Sampling_1500µm[[2]],na.rm = T),xlim=c(0,10),ylim=c(0,10),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="p parameter values (µm)",ylab="Mean estimated p values (µm)",main="1500µm FoV")
abline(0,1,lwd=2,lty=2,col='grey')
m = lm(colMeans(Sampling_1500µm[[2]],na.rm = T)~0+Results_sigmoid$p)
abline(0,coef(m),lwd=2,lty=2,col="black")
List_coef = c(List_coef,coef(m))
List_coef_sd = c(List_coef_sd,summary(m)$coefficients[2])

plot(Results_sigmoid$p,colMeans(Sampling_2000µm[[2]],na.rm = T),xlim=c(0,10),ylim=c(0,10),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="p parameter values (µm)",ylab="Mean estimated p values (µm)",main="2000µm FoV")
abline(0,1,lwd=2,lty=2,col='grey')
m = lm(colMeans(Sampling_2000µm[[2]],na.rm = T)~0+Results_sigmoid$p)
abline(0,coef(m),lwd=2,lty=2,col="black")
List_coef = c(List_coef,coef(m))
List_coef_sd = c(List_coef_sd,summary(m)$coefficients[2])


par(las=1,bty='l')
x =barplot(List_coef,ylim=c(0,10),names.arg = c("500µm","1000µm","1500µm","2000µm"),xlab="Fov size",ylab="Estimated/real p value slope",cex.lab=1.3)
segments(x0 = x,x1 =x,y0 = List_coef-List_coef_sd,y1 = List_coef+List_coef_sd)


par(las=1,bty='l')
barplot(c(CCC_computation(Results_sigmoid$p,colMeans(Sampling_500µm[[2]],na.rm = T)),
          CCC_computation(Results_sigmoid$p,colMeans(Sampling_1000µm[[2]],na.rm = T)),
          CCC_computation(Results_sigmoid$p,colMeans(Sampling_1500µm[[2]],na.rm = T)),
          CCC_computation(Results_sigmoid$p,colMeans(Sampling_2000µm[[2]],na.rm = T))),
        ylim=c(0,1),ylab="CCC coefficient",cex.lab=1.3,names.arg = c("500µm","1000µm","1500µm","2000µm"))


barplot(c(cor(Results_sigmoid$p,colMeans(Sampling_500µm[[2]],na.rm = T)),
          cor(Results_sigmoid$p,colMeans(Sampling_1000µm[[2]],na.rm = T)),
          cor(Results_sigmoid$p,colMeans(Sampling_1500µm[[2]],na.rm = T)),
          cor(Results_sigmoid$p,colMeans(Sampling_2000µm[[2]],na.rm = T))),
        ylim=c(0,1),ylab="Pearson's R",cex.lab=1.3,names.arg = c("500µm","1000µm","1500µm","2000µm"))



#E) Looking at the convergence of the estimator (even though biased...)

#i) For tau
Sd_estimator_matrix= c()
for (k in 1:ncol(Sampling_1000µm[[1]])) {
  temp_vector_sampling = Sampling_1000µm[[1]][,k]
  temp_vector_sampling = temp_vector_sampling[!is.na(temp_vector_sampling)]
  Sd_estimating_vector = c()
  for (i in 1:20) {
    u = apply(matrix(rep(1,500)),MARGIN = 1,FUN = function(x) {mean(temp_vector_sampling[sample(x = 1:length(temp_vector_sampling),size = i)])})
    Sd_estimating_vector = c(Sd_estimating_vector,sd(u))
  }
  Sd_estimator_matrix = cbind(Sd_estimator_matrix,Sd_estimating_vector)
}

CV_estimator_matrix = t(apply(Sd_estimator_matrix,MARGIN = 1,FUN = function(x) {x/colMeans(Sampling_1000µm[[1]],na.rm = T)}))

library(RColorBrewer)
N_cluster = ncol(CV_estimator_matrix)
optimal_palette =colorRamp(brewer.pal(11, "Spectral"))
optimal_palette = optimal_palette((1:N_cluster)/N_cluster)
optimal_palette = optimal_palette / 255
optimal_palette = rgb(optimal_palette)

par(las=1,bty='l')
plot(NULL,xlim=c(0,20),ylim=c(0,1.07),xaxs='i',yaxs='i',xlab="Number of FoVs sampled",ylab="CV",cex.lab=1.3)
for (k in 1:ncol(CV_estimator_matrix)) {
  x = 1:20
  y = CV_estimator_matrix[,k]
  points(x,y,pch=21,bg=optimal_palette[k])
  x_2 = 1/sqrt(x)
  m = lm(y~0+x_2)
  curve(expr = coef(m)/sqrt(x),from=0.1,to=20,add=T,col=optimal_palette[k],lwd=1.5)
  
}
abline(h=0.25,lwd=2,lty=2)

N_minimal_sampling_tau = apply(CV_estimator_matrix,MARGIN = 2,FUN = function(x) {min(which(x<0.25))})

#ii) For p


Sd_estimator_matrix= c()
for (k in 1:ncol(Sampling_1000µm[[2]])) {
  temp_vector_sampling = Sampling_1000µm[[2]][,k]
  temp_vector_sampling = temp_vector_sampling[!is.na(temp_vector_sampling)]
  Sd_estimating_vector = c()
  for (i in 1:20) {
    u = apply(matrix(rep(1,500)),MARGIN = 1,FUN = function(x) {mean(temp_vector_sampling[sample(x = 1:length(temp_vector_sampling),size = i)])})
    Sd_estimating_vector = c(Sd_estimating_vector,sd(u))
  }
  Sd_estimator_matrix = cbind(Sd_estimator_matrix,Sd_estimating_vector)
}

CV_estimator_matrix = t(apply(Sd_estimator_matrix,MARGIN = 1,FUN = function(x) {x/colMeans(Sampling_1000µm[[2]],na.rm = T)}))

library(RColorBrewer)
N_cluster = ncol(CV_estimator_matrix)
optimal_palette =colorRamp(brewer.pal(11, "Spectral"))
optimal_palette = optimal_palette((1:N_cluster)/N_cluster)
optimal_palette = optimal_palette / 255
optimal_palette = rgb(optimal_palette)

par(las=1,bty='l')
plot(NULL,xlim=c(0,20),ylim=c(0,4),xaxs='i',yaxs='i',xlab="Number of FoVs sampled",ylab="CV",cex.lab=1.3)
for (k in 1:ncol(CV_estimator_matrix)) {
  x = 1:20
  y = CV_estimator_matrix[,k]
  points(x,y,pch=21,bg=optimal_palette[k])
  x_2 = 1/sqrt(x)
  m = lm(y~0+x_2)
  curve(expr = coef(m)/sqrt(x),from=0.1,to=20,add=T,col=optimal_palette[k],lwd=1.5)
  
}
abline(h=0.25,lwd=2,lty=2)

N_minimal_sampling_p = apply(CV_estimator_matrix,MARGIN = 2,FUN = function(x) {min(which(x<0.25))})

#IV)Benchmark against other methods

#A) Comparison with Clark Evans index

CE_index = CE_interaction_tensor(sce)
CE_index_vector = 2^diag((CE_index[,,1])@data)

par(las=1,bty='l')
plot(CE_index_vector,(Results_sigmoid$C_normalised),xlab="CE index",ylab="C normalised",pch=21,bg="red3",cex=2,log='y',cex.lab=1.3)
m = lm(log10(Results_sigmoid$C_normalised)~CE_index_vector)
abline(coef(m),lty=2,col="grey",lwd=1.5)
m = summary(m)
legend("topright",legend = c(paste("R2=",round(m$r.squared,3)),paste("p=",round(m$coefficients[2,4],6))),bty='n')

par(las=1,bty='l')
plot(CE_index_vector,(Results_sigmoid$tau),xlab="CE index",ylab="Tau",pch=21,bg="red3",cex=2,log='',cex.lab=1.3)
m = lm((Results_sigmoid$tau)~CE_index_vector)
abline(coef(m),lty=2,col="grey",lwd=1.5)
m = summary(m)
legend("topright",legend = c(paste("R2=",round(m$r.squared,3)),paste("p=",round(m$coefficients[2,4],6))),bty='n')

par(las=1,bty='l')
plot(CE_index_vector,(Results_sigmoid$p),xlab="CE index",ylab="p",pch=21,bg="red3",cex=2,log='',cex.lab=1.3)
m = lm((Results_sigmoid$p)~CE_index_vector)
abline(coef(m),lty=2,col="grey",lwd=1.5)
m = summary(m)
legend("topright",legend = c(paste("R2=",round(m$r.squared,3)),paste("p=",round(m$coefficients[2,4],6))),bty='n')

barplot(c(cor(CE_index_vector,Results_sigmoid$tau,use = 'pairwise.complete'),
          cor(CE_index_vector,Results_sigmoid$p,use = 'pairwise.complete'),
          cor(CE_index_vector,Results_sigmoid$C_normalised,use = 'pairwise.complete')),ylim=c(-1,1),
        ylab="Correlation with CE index",xlab="Parameters",names.arg = c("tau","p","C_normalised"),cex.lab=1.3)


#B) Comparison with the maximal value of the normmalised L function

List_L_function =c()

for (k in List_pcf$Annotation$Cluster) {
  sce_temp = sce[,colLabels(sce)==k]
  ppp_temp = ppp(sce_temp$Location_Center_X,sce_temp$Location_Center_Y,window = owin(xrange = range(sce_temp$Location_Center_X),yrange = range(sce_temp$Location_Center_Y)))
  r_vector = seq(0,2000,length.out=50)
  Kest_temp = Kest.fft(ppp_temp,r = r_vector,rmax = max(r_vector),sigma = 1)
  List_L_function[[k]] = sqrt(Kest_temp$border/pi)-Kest_temp$r
}

plot(List_L_function[[5]]) #For most of the clusters : no local maxima

#C) Study of the over-dispersion index

library(MASS)
Overdispersion_vector = c()
for (k in List_pcf$Annotation$Cluster) {
  sce_temp = sce[,colLabels(sce)==k]
  ppp_temp = ppp(sce_temp$Location_Center_X,sce_temp$Location_Center_Y,window = owin(xrange = range(sce_temp$Location_Center_X),yrange = range(sce_temp$Location_Center_Y)))
  temp_count = as.numeric(quadratcount(ppp_temp,nx = 20,ny = 20))
  #parameter_fit_temp = fitdistr(x = temp_count,densfun = "negative binomial")
  parameter_fit_temp = fitdistrplus::fitdist(temp_count,distr  = "nbinom",method = "mm")
  Overdispersion_vector = c(Overdispersion_vector,coef(parameter_fit_temp)[1])
}

#Cluster 18 was removed : too rare and too concentrated
Overdispersion_vector = 1/Overdispersion_vector

par(las=1,bty='l')
plot(Overdispersion_vector[-18],Results_sigmoid$C_normalised[-18],xlab="Overdispersion index",ylab="C normalised",pch=21,bg="red3",cex=2,log='',cex.lab=1.3)
m = lm(Results_sigmoid$C_normalised[-18]~Overdispersion_vector[-18])
abline(coef(m),lty=2,col="grey",lwd=1.5)
m = summary(m)
legend("topleft",legend = c(paste("R2=",round(m$r.squared,3)),paste("p=",round(m$coefficients[2,4],10))),bty='n')

par(las=1,bty='l')
plot(Overdispersion_vector[-18],(Results_sigmoid$tau[-18]),xlab="Overdispersion index",ylab="Tau",pch=21,bg="red3",cex=2,log='',cex.lab=1.3)
m = lm(Results_sigmoid$tau[-18]~Overdispersion_vector[-18])
abline(coef(m),lty=2,col="grey",lwd=1.5)
m = summary(m)
legend("topright",legend = c(paste("R2=",round(m$r.squared,3)),paste("p=",round(m$coefficients[2,4],6))),bty='n')

par(las=1,bty='l')
plot(Overdispersion_vector[-18],Results_sigmoid$p[-18],xlab="Overdispersion index",ylab="p",pch=21,bg="red3",cex=2,log='',cex.lab=1.3)
m = lm(Results_sigmoid$p[-18]~Overdispersion_vector[-18])
abline(coef(m),lty=2,col="grey",lwd=1.5)
m = summary(m)
legend("topright",legend = c(paste("R2=",round(m$r.squared,3)),paste("p=",round(m$coefficients[2,4],6))),bty='n')

barplot(c(cor(Overdispersion_vector,Results_sigmoid$tau,use = 'pairwise.complete'),
          cor(Overdispersion_vector,Results_sigmoid$p,use = 'pairwise.complete'),
          cor(Overdispersion_vector,Results_sigmoid$C_normalised,use = 'pairwise.complete')),ylim=c(-1,1),
        ylab="Correlation with overdispersion index",xlab="Parameters",names.arg = c("tau","p","C_normalised"),cex.lab=1.3)



#Conclusion : only pcf based appraoch can capture the structure size....


#V) Testing the robustness of our approch to bad quality labelling 

#A) Simulation of 20% error of the clustering 
N_simulation = 100

Proportion_swapping = 0.2
List_parameter_table = c()
for (j in unique(colLabels(sce))) {
  Table_parameter_noise = c()
  Selected_cluster =j
  for (k in 1:N_simulation) {
    print(k)
    sce_temp = sce
    x = colLabels(sce_temp)
    x[x!=Selected_cluster] = 0
    x[x!=0] = 1
    N_cells_to_swapp = round(Proportion_swapping*length(which(x==1)))
    Cell_swap_positive = base::sample(which(x==1),size = N_cells_to_swapp,replace = F )
    Cell_swap_negative = base::sample(which(x==0),size = N_cells_to_swapp,replace = F )
    x[Cell_swap_positive] = 0
    x[Cell_swap_negative] = 1
    colLabels(sce_temp) = x
    List_pcf = Compute_pcf(sce_temp,r_vector = seq(0,2000,length.out=50),computation_method = "derivative")
    m_sigmoid = Fit_parametric_pcf_model(List_pcf,model = "Sigmoid")
    Table_parameter_noise = rbind(Table_parameter_noise,m_sigmoid[which(List_pcf$Annotation$Cluster==1),])
  }
  List_parameter_table[[j]] = Table_parameter_noise
}

#The corrupted PCF table computation can take some time to be computed. Can be saved as an RDS file
#List_parameter_table = readRDS("Tonsil_MERFISH_noise_effect_20_percent.rds")
Merged_parameter_table = c()

for (k in 1:length(List_parameter_table)) {
  u = List_parameter_table[[k]]
  u$Cluster = k
  Merged_parameter_table = rbind(Merged_parameter_table,u)
}
boxplot(Merged_parameter_table$tau~Merged_parameter_table$Cluster)
points(1:nrow(Results_sigmoid),Results_sigmoid$tau,pch=21,bg="red3")

#B) Plotting of the result 

#Concordance Correlation Coefficient 

par(las=1,bty="l")
Mean_tau_estimation = aggregate(Merged_parameter_table$tau,by=list(Merged_parameter_table$Cluster),FUN = mean)
Sd_tau_estimation = aggregate(Merged_parameter_table$tau,by=list(Merged_parameter_table$Cluster),FUN = sd)

plot(Results_sigmoid$tau,Mean_tau_estimation$x,xlim=c(0,500),ylim=c(0,500),xaxs='i',yaxs='i',cex=1.5,
     xlab="Initial tau estimation (µm)",ylab="Tau estimation with noise (µm)",pch=21,bg="red3",cex.lab=1.3)
abline(0,1,lty=2)
segments(x0 =Results_sigmoid$tau,x1 = Results_sigmoid$tau,y0 = Mean_tau_estimation$x + Sd_tau_estimation$x,
         y1 = Mean_tau_estimation$x - Sd_tau_estimation$x)
legend("bottomright",legend = paste("CCC = ",round(CCC_computation(Results_sigmoid$tau,Mean_tau_estimation$x),digits = 3)),bty='n')

Mean_p_estimation = aggregate(Merged_parameter_table$p,by=list(Merged_parameter_table$Cluster),FUN = mean)
Sd_p_estimation = aggregate(Merged_parameter_table$p,by=list(Merged_parameter_table$Cluster),FUN = sd)

plot(Results_sigmoid$p,Mean_p_estimation$x,xlim=c(0,3),ylim=c(0,3),xaxs='i',yaxs='i',cex=1.5,
     xlab="Initial p estimation",ylab="P estimation with noise",pch=21,bg="red3",cex.lab=1.3)
abline(0,1,lty=2)
segments(x0 =Results_sigmoid$p,x1 = Results_sigmoid$p,y0 = Mean_p_estimation$x + Sd_p_estimation$x,
         y1 = Mean_p_estimation$x - Sd_p_estimation$x)
legend("bottomright",legend = paste("CCC = ",round(CCC_computation(Results_sigmoid$p,Mean_p_estimation$x),digits = 3)),bty='n')

Mean_C_estimation = aggregate(Merged_parameter_table$C_normalised,by=list(Merged_parameter_table$Cluster),FUN = mean)
Sd_C_estimation = aggregate(Merged_parameter_table$C_normalised,by=list(Merged_parameter_table$Cluster),FUN = sd)

plot(Results_sigmoid$C_normalised,Mean_C_estimation$x,xlim=c(50,9000),ylim=c(50,9000),xaxs='i',yaxs='i',cex=1.5,
     xlab="Cnormalised p estimation",ylab="Cnormalised estimation with noise",pch=21,bg="red3",cex.lab=1.3,log="xy")
abline(0,1,lty=2)
segments(x0 =Results_sigmoid$C_normalised,x1 = Results_sigmoid$C_normalised,y0 = Mean_C_estimation$x + Sd_C_estimation$x,
         y1 = Mean_C_estimation$x - Sd_C_estimation$x)
legend("bottomright",legend = paste("CCC = ",round(CCC_computation(Results_sigmoid$C_normalised,Mean_C_estimation$x),digits = 3)),bty='n')


