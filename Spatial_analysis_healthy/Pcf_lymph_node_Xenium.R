
#I)Loading function/object
library(RColorBrewer)
CCC_computation = function(x,y) {
  return(2*sd(x,na.rm = T)*sd(y,na.rm = T)/(var(x,na.rm = T)+var(y,na.rm = T)+(mean(x,na.rm=T)-mean(y,na.rm=T))^2))
}

#A)List functions

source("Desktop/Projet_analyse_pcf/List_function_pcf_analysis.R")

Raw_metadata = read.delim("Annotation_lymph_node_Xenium.csv",sep=",")


sce = SingleCellExperiment(assays = list(Raw_intensity = matrix(0,nrow = 2,ncol =nrow(Raw_metadata) )), 
                           metadata = list(dimension = "2D", 
                                           N_core = 8, Is_nuc_cyt = F))
colLabels(sce) = as.numeric(Raw_metadata$Clustering)
sce$Location_Center_X = Raw_metadata$cell_centroid_x
sce$Location_Center_Y = Raw_metadata$cell_centroid_y
sce$ImageNumber =  1



#II)Computation itself 

#A)Computing pcf
List_pcf = Compute_pcf(sce,r_vector = seq(0,2000,length.out=50),computation_method = "derivative")


#B)Checking the functions on an example
k = 18
x = List_pcf$List_r[[k]]
y = List_pcf$List_pcf[[k]]

m = try(Fit_generalized_gamma(x,y,show_plot = TRUE),silent = TRUE)
m = try(Fit_sigmoid(x,y,show_plot = TRUE))

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


selected_clustering = 18
ppp_temp = ppp(x = Raw_metadata$cell_centroid_x[Raw_metadata$Clustering==selected_clustering],
               y = Raw_metadata$cell_centroid_y[Raw_metadata$Clustering==selected_clustering],
               window = owin(xrange = range(Raw_metadata$cell_centroid_x),yrange = range(Raw_metadata$cell_centroid_y)))
plot(ppp_temp,cex=0.2)



#III) Comparing the Pcf approach with more regular approaches

Plot_selected_cluster = function(selected_cluster,title_show=NULL) {
  ppp_temp = ppp(x = Raw_metadata$cell_centroid_x[Raw_metadata$Clustering==selected_cluster],
                 y = Raw_metadata$cell_centroid_y[Raw_metadata$Clustering==selected_cluster],
                 window = owin(xrange = range(Raw_metadata$cell_centroid_x),yrange = range(Raw_metadata$cell_centroid_y)))
  plot(ppp_temp,cex=0.2,main=title_show)
  
}

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

Overdispersion_vector = 1/Overdispersion_vector

par(las=1,bty='l')
plot(Overdispersion_vector,Results_sigmoid$C_normalised,xlab="Overdispersion index",ylab="C normalised",pch=21,bg="red3",cex=2,log='',cex.lab=1.3)
m = lm((Results_sigmoid$C_normalised)~Overdispersion_vector)
abline(coef(m),lty=2,col="grey",lwd=1.5)
m = summary(m)
legend("topleft",legend = c(paste("R2=",round(m$r.squared,3)),paste("p=",round(m$coefficients[2,4],10))),bty='n')

par(las=1,bty='l')
plot(Overdispersion_vector,(Results_sigmoid$tau),xlab="Overdispersion index",ylab="Tau",pch=21,bg="red3",cex=2,log='',cex.lab=1.3)
m = lm(Results_sigmoid$tau~Overdispersion_vector)
abline(coef(m),lty=2,col="grey",lwd=1.5)
m = summary(m)
legend("topright",legend = c(paste("R2=",round(m$r.squared,3)),paste("p=",round(m$coefficients[2,4],6))),bty='n')

par(las=1,bty='l')
plot(Overdispersion_vector,Results_sigmoid$p,xlab="Overdispersion index",ylab="p",pch=21,bg="red3",cex=2,log='',cex.lab=1.3)
m = lm(Results_sigmoid$p~Overdispersion_vector)
abline(coef(m),lty=2,col="grey",lwd=1.5)
m = summary(m)
legend("topright",legend = c(paste("R2=",round(m$r.squared,3)),paste("p=",round(m$coefficients[2,4],6))),bty='n')


barplot(c(cor(Overdispersion_vector,Results_sigmoid$tau,use = 'pairwise.complete'),
          cor(Overdispersion_vector,Results_sigmoid$p,use = 'pairwise.complete'),
          cor(Overdispersion_vector,Results_sigmoid$C_normalised,use = 'pairwise.complete')),ylim=c(-1,1),
        ylab="Correlation with overdispersion index",xlab="Parameters",names.arg = c("tau","p","C_normalised"),cex.lab=1.3)



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
Selected_cluster = 1
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
par(las=1,bty='l')
plot(Results_sigmoid$tau,colMeans(Sampling_500µm[[1]],na.rm = T),xlim=c(0,1000),ylim=c(0,1000),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="Tau parameter values (µm)",ylab="Mean estimated tau values (µm)",main="500µm FoV")
abline(0,1,lwd=2,lty=2,col="grey")
m = lm(colMeans(Sampling_500µm[[1]],na.rm = T)~0+Results_sigmoid$tau)
abline(0,coef(m),lwd=2,lty=2,col="black")
List_coef = c(List_coef,coef(m))

plot(Results_sigmoid$tau,colMeans(Sampling_1000µm[[1]],na.rm = T),xlim=c(0,1000),ylim=c(0,1000),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="Tau parameter values (µm)",ylab="Mean estimated tau values (µm)",main="1000µm FoV")
abline(0,1,lwd=2,lty=2,col="grey")
m = lm(colMeans(Sampling_1000µm[[1]],na.rm = T)~0+Results_sigmoid$tau)
abline(0,coef(m),lwd=2,lty=2,col="black")
List_coef = c(List_coef,coef(m))

plot(Results_sigmoid$tau,colMeans(Sampling_1500µm[[1]],na.rm = T),xlim=c(0,1000),ylim=c(0,1000),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="Tau parameter values (µm)",ylab="Mean estimated tau values (µm)",main="1500µm FoV")
abline(0,1,lwd=2,lty=2,col="grey")
m = lm(colMeans(Sampling_1500µm[[1]],na.rm = T)~0+Results_sigmoid$tau)
abline(0,coef(m),lwd=2,lty=2,col="black")
List_coef = c(List_coef,coef(m))

plot(Results_sigmoid$tau,colMeans(Sampling_2000µm[[1]],na.rm = T),xlim=c(0,1000),ylim=c(0,1000),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="Tau parameter values (µm)",ylab="Mean estimated tau values (µm)",main="2000µm FoV")
abline(0,1,lwd=2,lty=2,col="grey")
m = lm(colMeans(Sampling_2000µm[[1]],na.rm = T)~0+Results_sigmoid$tau)
abline(0,coef(m),lwd=2,lty=2,col="black")
List_coef = c(List_coef,coef(m))

par(las=1,bty='l')
barplot(c(CCC_computation(Results_sigmoid$tau,colMeans(Sampling_500µm[[1]],na.rm = T)),
          CCC_computation(Results_sigmoid$tau,colMeans(Sampling_1000µm[[1]],na.rm = T)),
          CCC_computation(Results_sigmoid$tau,colMeans(Sampling_1500µm[[1]],na.rm = T)),
          CCC_computation(Results_sigmoid$tau,colMeans(Sampling_2000µm[[1]],na.rm = T))),
        ylim=c(0,1),ylab="CCC coefficient",cex.lab=1.3,names.arg = c("500µm","1000µm","1500µm","2000µm"))

barplot(c(cor(Results_sigmoid$tau,colMeans(Sampling_500µm[[1]],na.rm = T),use = 'pairwise.complete'),
          cor(Results_sigmoid$tau,colMeans(Sampling_1000µm[[1]],na.rm = T),use = 'pairwise.complete'),
          cor(Results_sigmoid$tau,colMeans(Sampling_1500µm[[1]],na.rm = T),use = 'pairwise.complete'),
          cor(Results_sigmoid$tau,colMeans(Sampling_2000µm[[1]],na.rm = T),use = 'pairwise.complete')),
        ylim=c(0,1),ylab="Pearson's R",cex.lab=1.3,names.arg = c("500µm","1000µm","1500µm","2000µm"))

barplot(List_coef,ylim=c(0,1),ylab="Slope coefficient",cex.lab=1.3,names.arg = c("500µm","1000µm","1500µm","2000µm"))

#ii) For p parameter

par(las=1,bty='l')
plot(Results_sigmoid$p,colMeans(Sampling_500µm[[2]],na.rm = T),xlim=c(0,10),ylim=c(0,10),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="p parameter values (µm)",ylab="Mean estimated p values (µm)",main="500µm FoV")
abline(0,1,lwd=2,lty=2,col="grey")
m = lm(colMeans(Sampling_500µm[[2]],na.rm = T)~0+Results_sigmoid$p)
abline(0,coef(m),lwd=2,lty=2,col="black")
plot(Results_sigmoid$p,colMeans(Sampling_1000µm[[2]],na.rm = T),xlim=c(0,10),ylim=c(0,10),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="p parameter values (µm)",ylab="Mean estimated p values (µm)",main="1000µm FoV")
abline(0,1,lwd=2,lty=2,col='grey')
m = lm(colMeans(Sampling_1000µm[[2]],na.rm = T)~0+Results_sigmoid$p)
abline(0,coef(m),lwd=2,lty=2,col="black")
plot(Results_sigmoid$p,colMeans(Sampling_1500µm[[2]],na.rm = T),xlim=c(0,10),ylim=c(0,10),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="p parameter values (µm)",ylab="Mean estimated p values (µm)",main="1500µm FoV")
abline(0,1,lwd=2,lty=2,col='grey')
m = lm(colMeans(Sampling_1500µm[[2]],na.rm = T)~0+Results_sigmoid$p)
abline(0,coef(m),lwd=2,lty=2,col="black")
plot(Results_sigmoid$p,colMeans(Sampling_2000µm[[2]],na.rm = T),xlim=c(0,10),ylim=c(0,10),xaxs='i',yaxs='i',pch=21,bg="red3",cex=2,cex.lab=1.3,xlab="p parameter values (µm)",ylab="Mean estimated p values (µm)",main="2000µm FoV")
abline(0,1,lwd=2,lty=2,col='grey')
m = lm(colMeans(Sampling_2000µm[[2]],na.rm = T)~0+Results_sigmoid$p)
abline(0,coef(m),lwd=2,lty=2,col="black")

par(las=1,bty='l')
barplot(c(CCC_computation(Results_sigmoid$p,colMeans(Sampling_500µm[[2]],na.rm = T)),
          CCC_computation(Results_sigmoid$p,colMeans(Sampling_1000µm[[2]],na.rm = T)),
          CCC_computation(Results_sigmoid$p,colMeans(Sampling_1500µm[[2]],na.rm = T)),
          CCC_computation(Results_sigmoid$p,colMeans(Sampling_2000µm[[2]],na.rm = T))),
        ylim=c(0,1),ylab="CCC coefficient",cex.lab=1.3,names.arg = c("500µm","1000µm","1500µm","2000µm"))


barplot(c(cor(Results_sigmoid$p,colMeans(Sampling_500µm[[2]],na.rm = T),use = "pairwise.complete"),
          cor(Results_sigmoid$p,colMeans(Sampling_1000µm[[2]],na.rm = T),use = "pairwise.complete"),
          cor(Results_sigmoid$p,colMeans(Sampling_1500µm[[2]],na.rm = T),use = "pairwise.complete"),
          cor(Results_sigmoid$p,colMeans(Sampling_2000µm[[2]],na.rm = T),use = "pairwise.complete")),
        ylim=c(0,1),ylab="Pearson's R",cex.lab=1.3,names.arg = c("500µm","1000µm","1500µm","2000µm"))



#E) Computing coefficient of variation

#i) For tau
CV_500µm_tau = apply(Sampling_500µm[[1]],MARGIN = 2,FUN = function(x) {sd(x,na.rm = T)/mean(x,na.rm=T)})
CV_1000µm_tau = apply(Sampling_1000µm[[1]],MARGIN = 2,FUN = function(x) {sd(x,na.rm = T)/mean(x,na.rm=T)})
CV_1500µm_tau = apply(Sampling_1500µm[[1]],MARGIN = 2,FUN = function(x) {sd(x,na.rm = T)/mean(x,na.rm=T)})
CV_2000µm_tau = apply(Sampling_2000µm[[1]],MARGIN = 2,FUN = function(x) {sd(x,na.rm = T)/mean(x,na.rm=T)})


boxplot(CV_500µm_tau,CV_1000µm_tau,CV_1500µm_tau,CV_2000µm_tau,names = c("500µm","1000µm","1500µm","2000µm"),outline=F,ylim=c(0,1),yaxs='i')

#ii) For p
CV_500µm_p = apply(Sampling_500µm[[2]],MARGIN = 2,FUN = function(x) {sd(x,na.rm = T)/mean(x,na.rm=T)})
CV_1000µm_p = apply(Sampling_1000µm[[2]],MARGIN = 2,FUN = function(x) {sd(x,na.rm = T)/mean(x,na.rm=T)})
CV_1500µm_p = apply(Sampling_1500µm[[2]],MARGIN = 2,FUN = function(x) {sd(x,na.rm = T)/mean(x,na.rm=T)})
CV_2000µm_p = apply(Sampling_2000µm[[2]],MARGIN = 2,FUN = function(x) {sd(x,na.rm = T)/mean(x,na.rm=T)})

boxplot(CV_500µm_p,CV_1000µm_p,CV_1500µm_p,CV_2000µm_p,names = c("500µm","1000µm","1500µm","2000µm"),outline=F,ylim=c(0,3),yaxs='i')


