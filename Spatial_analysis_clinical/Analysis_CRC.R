library(fifer)
library(SingleCellExperiment)
library(RColorBrewer)
library(balagan)
library(pheatmap)
color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("grey","yellow","orange","red"))
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

#I)Loading function/object

#A)List functions

source("Desktop/Projet_analyse_pcf/List_function_pcf_analysis.R")

Raw_metadata = read.delim("Annotated_CRC_data.txt",sep="\t")


Cluster_temp = Raw_metadata$Cluster
Cluster_temp[Cluster_temp%in%c(7,9)] = "Macrophages"
Cluster_temp[Cluster_temp%in%c(14)] = "PD-L1+ cells"
Cluster_temp[Cluster_temp%in%c(16)] = "CD8a+ T cells"
Cluster_temp[Cluster_temp%in%c(10)] = "CD4+ T cells"
Cluster_temp[Cluster_temp%in%c(18)] = "B cells"
Cluster_temp[Cluster_temp%in%c(13)] = "Treg"
Cluster_temp[Cluster_temp%in%c(1)] = "Endothelial cells"
Cluster_temp[Cluster_temp%in%c(6)] = "Vessel cells"
Cluster_temp[Cluster_temp%in%c(11,3)] = "Smooth muscular cells"
Cluster_temp[Cluster_temp%in%c(5,8)] = "Fibroblast"
Cluster_temp[Cluster_temp%in%c(4)] = "Cancer 1"
Cluster_temp[Cluster_temp%in%c(2)] = "Cancer 2"
Cluster_temp[Cluster_temp%in%c(15)] = "Cancer 3"
Cluster_temp[Cluster_temp%in%c(19)] = "Cancer 4"
Cluster_temp[Cluster_temp%in%c(12,17)] = "Other"

sce = SingleCellExperiment(assays = list(Raw_intensity = matrix(0,nrow = 2,ncol =nrow(Raw_metadata) )), 
                           metadata = list(dimension = "2D", 
                                           N_core = 8, Is_nuc_cyt = F))
colLabels(sce) = (Cluster_temp)
sce$Location_Center_X = Raw_metadata$Xt
sce$Location_Center_Y = Raw_metadata$Yt
sce$ImageNumber =  Raw_metadata$Sample



Patient_metadata = read.delim("Desktop/CRC_pcf_analysis/Metadata_patient.txt",row.names = 1)
#removing CRC01
Patient_metadata = Patient_metadata[-1,]
#II)Computation itself 

#A)Computing pcf
List_pcf = Compute_pcf(sce,r_vector = seq(0,5000,length.out=50),computation_method = "derivative",verbose = TRUE)


#B)Checking the functions on an example

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

Table_fits = Results_sigmoid
Table_fits$Cluster = List_pcf$Annotation$Cluster
Table_fits$Sample = List_pcf$Annotation$ROI
#C) Computing mean behavior of each cell type 

Mean_cluster_tau = aggregate(Table_fits$tau,by=list(Table_fits$Cluster),FUN =function(x) {mean(x[is.finite(x)])})

Mean_cluster_C = aggregate(Table_fits$C_normalised,by=list(Table_fits$Cluster),FUN =function(x) {mean(x[is.finite(x)])})

Mean_cluster_p = aggregate(Table_fits$p,by=list(Table_fits$Cluster),FUN =function(x) {mean(x[is.finite(x)])})

N_cell_cluster = table(colLabels(sce))
N_cell_cluster = N_cell_cluster/sum(N_cell_cluster)*100

par(las=1,bty='l')
plot(Mean_cluster_tau$x,Mean_cluster_p$x,cex=sqrt(N_cell_cluster)*1.3,
     xlab="Tau parameter (µm)",ylab="p parameter",pch=21,log='',
     bg="red3",cex.lab=1.3,xlim=c(1,1300),ylim=c(0.5,3),xaxs='i',yaxs='i')
text(Mean_cluster_tau$x,Mean_cluster_p$x,Mean_cluster_p$Group.1)



#V)Looking at the spatial structure of multipke cell types

#A)Reshaping tau, p and C normalised

List_clusters = unique(colLabels(sce))
Table_fit_reshape_tau = c()
for (k in unique(Table_fits$Sample)) {
  x = Table_fits[Table_fits$Sample==k,]
  rownames(x) = x$Cluster
  x = x[as.character(List_clusters),]
  Table_fit_reshape_tau = rbind(Table_fit_reshape_tau,x$tau)
}

colnames(Table_fit_reshape_tau) = List_clusters
Table_fit_reshape_tau = as.data.frame(Table_fit_reshape_tau)
Table_fit_reshape_tau = Table_fit_reshape_tau[,colSums(!is.na(Table_fit_reshape_tau))>=5]


List_clusters = unique(colLabels(sce))
Table_fit_reshape_p = c()
for (k in unique(Table_fits$Sample)) {
  x = Table_fits[Table_fits$Sample==k,]
  rownames(x) = x$Cluster
  x = x[as.character(List_clusters),]
  Table_fit_reshape_p = rbind(Table_fit_reshape_p,x$p)
}

colnames(Table_fit_reshape_p) = List_clusters
Table_fit_reshape_p = as.data.frame(Table_fit_reshape_p)
Table_fit_reshape_p = Table_fit_reshape_p[,colSums(!is.na(Table_fit_reshape_p))>=5]

#B)Computing the robust correlation

library(WRS2)

Robust_correlation = matrix(0,ncol = ncol(Table_fit_reshape_tau),nrow = ncol(Table_fit_reshape_tau))
Robust_correlation_pvalue = matrix(1,ncol = ncol(Table_fit_reshape_tau),nrow = ncol(Table_fit_reshape_tau))
List_interaction = matrix(NA,ncol = ncol(Table_fit_reshape_tau),nrow = ncol(Table_fit_reshape_tau))
N_shared_observation = matrix(0,ncol = ncol(Table_fit_reshape_tau),nrow = ncol(Table_fit_reshape_tau))
for (i in 1:nrow(Robust_correlation)) {
  for (j in 1:ncol(Robust_correlation)) {
    x = Table_fit_reshape_tau[,i]
    y = Table_fit_reshape_tau[,j]
    shared_obs = !is.na(x) & !is.na(y)
    if ( sum(shared_obs)>1) {
      u = pball(cbind(x[shared_obs],y[shared_obs]))
      Robust_correlation[i,j] = u$pbcorm[1,2]
      Robust_correlation_pvalue [i,j] = u$p.values[1,2]
    }
    List_interaction[i,j] = paste(colnames(Table_fit_reshape_tau)[i],colnames(Table_fit_reshape_tau)[j])
    N_shared_observation[i,j]=sum(shared_obs)
  }
}
colnames(Robust_correlation) = colnames(Table_fit_reshape_tau)
rownames(Robust_correlation) = colnames(Table_fit_reshape_tau)
pheatmap((Robust_correlation),clustering_method = 'ward.D2')


Corrected_p_value = -log10(p.adjust(as.numeric(Robust_correlation_pvalue),method = "fdr"))
Corrected_p_value[is.na(Corrected_p_value)] = 0
Corrected_p_value[is.infinite(Corrected_p_value)]=0

#C) Generating plots 

#Volcano

par(las=1,bty='l')
plot(as.numeric(Robust_correlation),Corrected_p_value,
     cex=as.numeric(N_shared_observation)/5,pch=21,
     bg=string.to.colors(Corrected_p_value>2,colors = c("grey","red3")),
     xlim=c(-1,1),xaxs='i',ylim=c(0,3),yaxs='i',xlab='Robust correlation coefficient',cex.lab=1.3,
     ylab="-Log10(fdr")
abline(h=2,lty=2,col="red3")
plot(as.numeric(Robust_correlation),Corrected_p_value,
     cex=as.numeric(N_shared_observation)/5,pch=21,
     bg=string.to.colors(Corrected_p_value>2,colors = c("grey","red3")),
     xlim=c(-1,1),xaxs='i',ylim=c(0,3),yaxs='i',xlab='Robust correlation coefficient',cex.lab=1.3,
     ylab="-Log10(fdr")
abline(h=2,lty=2,col="red3")
List_interaction[lower.tri(Robust_correlation)]=NA
text(as.numeric(Robust_correlation)[Corrected_p_value>2],Corrected_p_value[Corrected_p_value>2],as.character(List_interaction)[Corrected_p_value>2])



View(data.frame(as.numeric(Robust_correlation),-log10(p.adjust(as.numeric(Robust_correlation_pvalue),method = 'fdr')),as.character(List_interaction)))


par(las=1,bty='l')
plot(Table_fit_reshape_tau$`CD4+ T cells`,Table_fit_reshape_tau$Treg,
     xlim=c(0,1000),ylim=c(0,1200),xaxs='i',yaxs='i',
     xlab="tau parameter in CD4+ T cells",
     ylab="tau parameter in Treg cells",cex=2,bg="red3",pch=21,cex.lab=1.3)
m = lm(Table_fit_reshape_tau$Treg~Table_fit_reshape_tau$`CD4+ T cells`)
abline(coef(m),lwd=2,lty=2,col="black")


plot(Table_fit_reshape_tau$`Cancer 2`,Table_fit_reshape_tau$`Vessel cells`,
     xlim=c(0,4000),ylim=c(0,3000),xaxs='i',yaxs='i',
     xlab="tau parameter in Cancer 2 cells",
     ylab="tau parameter in Vessel cells",cex=2,bg="red3",pch=21,cex.lab=1.3)
m = lm(Table_fit_reshape_tau$`Vessel cells`~Table_fit_reshape_tau$`Cancer 2`)
abline(coef(m),lwd=2,lty=2,col="black")


plot(Table_fit_reshape_tau$`CD8a+ T cells`,Table_fit_reshape_tau$`Smooth muscular cells`,
     xlim=c(0,2000),ylim=c(0,4000),xaxs='i',yaxs='i',
     ylab="tau parameter in smooth muscular cells",
     xlab="tau parameter in CD8a+ T cells",cex=2,bg="red3",pch=21,cex.lab=1.3)
m = lm(Table_fit_reshape_tau$`Smooth muscular cells`~Table_fit_reshape_tau$`CD8a+ T cells`)
abline(coef(m),lwd=2,lty=2,col="black")

plot(Table_fit_reshape_tau$`CD8a+ T cells`,Table_fit_reshape_tau$`CD4+ T cells`,
     xlim=c(0,2000),ylim=c(0,1000),xaxs='i',yaxs='i',
     ylab="tau parameter in CD4+ T cells",
     xlab="tau parameter in CD8a+ T cells",cex=2,bg="red3",pch=21,cex.lab=1.3)
m = lm(Table_fit_reshape_tau$`CD4+ T cells`~Table_fit_reshape_tau$`CD8a+ T cells`)
abline(coef(m),lwd=2,lty=2,col="black")



#V) Co-organisation analysis 

#A) Computation

Robust_correlation = matrix(0,ncol = ncol(Table_fit_reshape_p),nrow = ncol(Table_fit_reshape_p))
Robust_correlation_pvalue = matrix(1,ncol = ncol(Table_fit_reshape_p),nrow = ncol(Table_fit_reshape_p))
List_interaction = matrix(NA,ncol = ncol(Table_fit_reshape_p),nrow = ncol(Table_fit_reshape_p))
N_shared_observation = matrix(0,ncol = ncol(Table_fit_reshape_p),nrow = ncol(Table_fit_reshape_p))
for (i in 1:nrow(Robust_correlation)) {
  for (j in 1:ncol(Robust_correlation)) {
    x = Table_fit_reshape_p[,i]
    y = Table_fit_reshape_p[,j]
    shared_obs = !is.na(x) & !is.na(y)
    if ( sum(shared_obs)>1) {
      u = pball(cbind(x[shared_obs],y[shared_obs]))
      Robust_correlation[i,j] = u$pbcorm[1,2]
      Robust_correlation_pvalue [i,j] = u$p.values[1,2]
    }
    List_interaction[i,j] = paste(colnames(Table_fit_reshape_p)[i],colnames(Table_fit_reshape_p)[j])
    N_shared_observation[i,j]=sum(shared_obs)
  }
}
colnames(Robust_correlation) = colnames(Table_fit_reshape_p)
rownames(Robust_correlation) = colnames(Table_fit_reshape_p)
pheatmap((Robust_correlation),clustering_method = 'ward.D2')


Corrected_p_value = -log10(p.adjust(as.numeric(Robust_correlation_pvalue),method = "fdr"))
Corrected_p_value[is.na(Corrected_p_value)] = 0
Corrected_p_value[is.infinite(Corrected_p_value)]=0

#B) Generating plots 

#Volcano

par(las=1,bty='l')
plot(as.numeric(Robust_correlation),Corrected_p_value,
     cex=as.numeric(N_shared_observation)/5,pch=21,
     bg=string.to.colors(Corrected_p_value>2,colors = c("grey","red3")),
     xlim=c(-1,1),xaxs='i',ylim=c(0,10),yaxs='i',xlab='Robust correlation coefficient',cex.lab=1.3,
     ylab="-Log10(fdr")
abline(h=2,lty=2,col="red3")
plot(as.numeric(Robust_correlation),Corrected_p_value,
     cex=as.numeric(N_shared_observation)/5,pch=21,
     bg=string.to.colors(Corrected_p_value>2,colors = c("grey","red3")),
     xlim=c(-1,1),xaxs='i',ylim=c(0,3),yaxs='i',xlab='Robust correlation coefficient',cex.lab=1.3,
     ylab="-Log10(fdr")
abline(h=2,lty=2,col="red3")
List_interaction[lower.tri(Robust_correlation)]=NA
text(as.numeric(Robust_correlation)[Corrected_p_value>2],Corrected_p_value[Corrected_p_value>2],as.character(List_interaction)[Corrected_p_value>2])


par(las=1,bty='l')
plot(Table_fit_reshape_p$`CD4a+ T cells`,Table_fit_reshape_p$`Endothelial cells`,
     xlim=c(0.5,2),ylim=c(0.5,2),xaxs='i',yaxs='i',xlab="P parameter in CD4+ T cells",
     ylab="P parameter in endothelial cells",cex=2,bg="red3",pch=21,cex.lab=1.3)
m = lm(Table_fit_reshape_p$`Endothelial cells`~Table_fit_reshape_p$`CD4a+ T cells`)
abline(coef(m),lwd=2,lty=2,col="black")
plot(Table_fit_reshape_p$`CD4a+ T cells`,Table_fit_reshape_p$`CD8a+ T cells`,
     xlim=c(0.5,2.5),ylim=c(0.5,2.5),xaxs='i',yaxs='i',xlab="P parameter in CD4+ T cells",
     ylab="P parameter in CD8a+ T cells",cex=2,bg="red3",pch=21,cex.lab=1.3)
m = lm(Table_fit_reshape_p$`CD8a+ T cells`~Table_fit_reshape_p$`CD4a+ T cells`)
abline(coef(m),lwd=2,lty=2,col="black")


#VII)Plotting examples of cell type co-scaling

Plot_two_cell_types = function(Sample="CRC02",Cell_type_1 = "CD8a+ T cells",Cell_type_2 = "CD4a+ T cells") {
  Temp_sce = sce[,sce$ImageNumber==Sample]
  range_temp = range(c(Temp_sce$Location_Center_X,Temp_sce$Location_Center_Y))
  par(las=1,bty='l')
  plot(Temp_sce$Location_Center_X,Temp_sce$Location_Center_Y,xlim=range_temp,ylim=range_temp,cex=0.06,pch=16,col="grey",
       xlab="X location (µm)",ylab="Y location (µm)",xaxs='i',yaxs='i',main=Sample)
  points(Temp_sce$Location_Center_X[colLabels(Temp_sce)==Cell_type_1],
         Temp_sce$Location_Center_Y[colLabels(Temp_sce)==Cell_type_1],cex=0.06,col="orange")
  points(Temp_sce$Location_Center_X[colLabels(Temp_sce)==Cell_type_2],
         Temp_sce$Location_Center_Y[colLabels(Temp_sce)==Cell_type_2],cex=0.06,col="darkblue")
  
}

#A) For co-scaling

#CD4 T cells / Treg
Plot_two_cell_types("CRC08",Cell_type_1 = "CD4a+ T cells",Cell_type_2 = "Treg")

#CD8+ T cells / smooth muscular cells 

Plot_two_cell_types("CRC09",Cell_type_2 = "CD8a+ T cells",Cell_type_1= "Smooth muscular cells")

#B) For co-organisation

Plot_two_cell_types("CRC08",Cell_type_2 = "CD4a+ T cells",Cell_type_1= "CD8a+ T cells")
Plot_two_cell_types("CRC03",Cell_type_2 = "CD4a+ T cells",Cell_type_1= "CD8a+ T cells")
Plot_two_cell_types("CRC10",Cell_type_2 = "CD4a+ T cells",Cell_type_1= "Endothelial cells")
Plot_two_cell_types("CRC03",Cell_type_2 = "CD4a+ T cells",Cell_type_1= "Endothelial cells")


