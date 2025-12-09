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

#I)Loading function/object and checking cell composition

#A)List functions

source("Desktop/Projet_analyse_pcf/List_function_pcf_analysis.R")

Raw_metadata = read.delim("Updated_metadata_UM.csv",sep=",")

Cluster_temp = Raw_metadata$Clustering
Cluster_temp[Cluster_temp%in%c(4)] = "Healthy hepatocytes"
Cluster_temp[Cluster_temp%in%c(15)] = "Inflammed hepatocytes"
Cluster_temp[Cluster_temp%in%c(13)] = "Neutrophils"
Cluster_temp[Cluster_temp%in%c(7)] = "Kupfer cells"
Cluster_temp[Cluster_temp%in%c(16)] = "Plasma cells"
Cluster_temp[Cluster_temp%in%c(10)] = "CD8+ T cells"
Cluster_temp[Cluster_temp%in%c(6)] = "Unannotated cells"
Cluster_temp[Cluster_temp%in%c(5)] = "Cancer cells 6"
Cluster_temp[Cluster_temp%in%c(2)] = "Cancer cells 5"
Cluster_temp[Cluster_temp%in%c(17)] = "Cholangiocytes"
Cluster_temp[Cluster_temp%in%c(18)] = "Lymphatic endothelial cells"
Cluster_temp[Cluster_temp%in%c(11)] = "Hepatic stellate cells"
Cluster_temp[Cluster_temp%in%c(9)] = "CCN+ cells"
Cluster_temp[Cluster_temp%in%c(3)] = "Fibroblasts"
Cluster_temp[Cluster_temp%in%c(8)] = "Dividing cancer cells"
Cluster_temp[Cluster_temp%in%c(14)] = "Cancer cells 1"
Cluster_temp[Cluster_temp%in%c(19)] = "Cancer cells 2"
Cluster_temp[Cluster_temp%in%c(12)] = "Cancer cells 3"
Cluster_temp[Cluster_temp%in%c(1)] = "Cancer cells 4"


sce = SingleCellExperiment(assays = list(Raw_intensity = matrix(0,nrow = 2,ncol =nrow(Raw_metadata) )), 
                           metadata = list(dimension = "2D", 
                                           N_core = 8, Is_nuc_cyt = F))
colLabels(sce) = (Cluster_temp)

sce$Location_Center_X = Raw_metadata$cell_centroid_x
sce$Location_Center_Y = Raw_metadata$cell_centroid_y
sce$ImageNumber =  Raw_metadata$sample

#B) Studyin composition

table_cell_count = table(sce$ImageNumber,colLabels(sce))
table_cell_count_nornalised = table_cell_count/rowSums(table_cell_count)
table_cell_count_nornalised = as.data.frame.matrix(table_cell_count_nornalised)
U = c()
for (k in colnames(table_cell_count_nornalised)) {
  x = table_cell_count_nornalised[,k]
  U = rbind(U,cbind(x,rep(k,nrow(table_cell_count_nornalised))))
}
U = as.data.frame(U)
colnames(U) = c("Proportion","Cell_type")
U$Proportion = as.numeric(U$Proportion)*100
Mean_proportion = sapply(table_cell_count_nornalised,FUN = median)
U$Cell_type = factor(U$Cell_type,names(Mean_proportion)[order(Mean_proportion)])

par(las=1,bty='l')
plot(NULL,xlim=c(0.5,20.5),ylim=c(0,51),log='',xaxs='i',yaxs='i',xaxt='n',ylab="Cell proportion type (%)",xlab="Cell type",cex.lab=1.3)
prism.plots(Proportion~Cell_type,U,add=T,cex.axis=0.1,pch=21,bg="red3",col='black',cex=1.5)

table_cell_count_nornalised_cell_type = t(table_cell_count)/rowSums(t(table_cell_count))
entropy_cell_type = apply(as.data.frame.matrix(table_cell_count_nornalised_cell_type),MARGIN = 1,function(x) {-sum(x*log(x),na.rm = T)})
U = data.frame(Entropy = entropy_cell_type,Type=names(entropy_cell_type),Is_cancer = ifelse(grepl(names(entropy_cell_type),pattern = "Cancer",ignore.case = T),yes = "Cancer",no = "Non cancer"))
U = U[order(U$Is_cancer),]

par(las=1,bty='l')
plot(NULL,xlim=c(0.7,2.3),ylim=c(0,2.1),xaxt='n',xaxs='i',yaxs='i',xlab="Type",ylab="Entropy",cex.lab=1.3)
prism.plots(Entropy~Is_cancer,U,add=T,cex=2,pch=21,bg=string.to.colors(U$Is_cancer,colors = c("red3","grey")),col='black')

#II)Computation itself 

#A)Computing pcf
List_pcf = Compute_pcf(sce,r_vector = seq(0,1000,length.out=50),computation_method = "derivative",verbose = TRUE)


#B)Checking the functions on an example

Results_generalised_gamma = Fit_parametric_pcf_model(List_pcf,model = "Generalized_gamma")
Results_gamma= Fit_parametric_pcf_model(List_pcf,model = "Gamma")
Results_exponential = Fit_parametric_pcf_model(List_pcf,model = "Exponential")
Results_power_law = Fit_parametric_pcf_model(List_pcf,model = "Power_law")
Results_sigmoid= Fit_parametric_pcf_model(List_pcf,model = "Sigmoid")
Results_sigmoid_robust= Fit_parametric_pcf_model(List_pcf,model = "Sigmoid_robust")
Results_beta_prime= Fit_parametric_pcf_model(List_pcf,model = "Beta_prime")

par(las=1,bty='l')
boxplot(Results_power_law$R2,Results_exponential$R2,Results_gamma$R2,
        Results_generalised_gamma$R2,Results_sigmoid$R2,Results_beta_prime$R2,ylim=c(0.4,1),
        xaxs='i',yaxs='i',names = c("Power law","Exponential","Gamma","GG","Sigmoid","Beta prime"),
        xlab="Model",ylab="R2",cex.lab=1.3,col=brewer.pal(6, "Spectral"),cex.axis=0.7)

boxplot(Results_sigmoid$R2,Results_sigmoid_robust$R2)
Table_fits = Results_sigmoid_robust
Table_fits$Cluster = List_pcf$Annotation$Cluster
Table_fits$Sample = List_pcf$Annotation$ROI

#C) Computing mean behavior of each cell type 

Mean_cluster_tau = aggregate(Table_fits$tau,by=list(Table_fits$Cluster),FUN =function(x) {median(x[is.finite(x)])})

Mean_cluster_C = aggregate(Table_fits$C_normalised,by=list(Table_fits$Cluster),FUN =function(x) {median(x[is.finite(x)])})

Mean_cluster_p = aggregate(Table_fits$p,by=list(Table_fits$Cluster),FUN =function(x) {median(x[is.finite(x)])})

N_cell_cluster = table(colLabels(sce))
N_cell_cluster = N_cell_cluster/sum(N_cell_cluster)*100

#D) Plots

par(las=1,bty='l')
plot(Mean_cluster_tau$x,Mean_cluster_p$x,cex=sqrt(N_cell_cluster),
     xlab="Tau parameter (µm)",ylab="p parameter",pch=21,log='',
     bg="red3",cex.lab=1.3,xaxs='i',yaxs='i',xlim=c(0,2500),ylim=c(0.5,3))
text(Mean_cluster_tau$x,Mean_cluster_p$x,Mean_cluster_p$Group.1)


#III)Looking at individual spatial patterns 

#A)Prioriterizing the cell types to study

Median_cluster_C = aggregate(Table_fits$C_normalised,by=list(Table_fits$Cluster),FUN =function(x) {median(x[is.finite(x)])})


Order_cluster_ranked_C = Median_cluster_C$Group.1[order(Median_cluster_C$x,decreasing = TRUE)]

U = Table_fits
U$Cluster = factor(U$Cluster,levels = Order_cluster_ranked_C)
par(las=1,bty='l')
plot(NULL,xlim=c(0.5,20.5),ylim=c(1,max(U$C_normalised,na.rm = T)*1.05),log='y',xaxs='i',yaxs='i',xaxt='n',ylab="Cnormalised",xlab="Cell type")
prism.plots(C_normalised~Cluster,U,add=T,cex.axis=0.1,pch=21,bg="red3",col='black',cex=1.5)


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


Table_fit_reshape_C = c()
for (k in unique(Table_fits$Sample)) {
  x = Table_fits[Table_fits$Sample==k,]
  rownames(x) = x$Cluster
  x = x[as.character(List_clusters),]
  Table_fit_reshape_C = rbind(Table_fit_reshape_C,x$C_normalised)
}

colnames(Table_fit_reshape_C) = List_clusters
Table_fit_reshape_C = as.data.frame(Table_fit_reshape_C)
Table_fit_reshape_C = Table_fit_reshape_C[,colSums(!is.na(Table_fit_reshape_C))>=5]


#C)Computing the robust correlation

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

#D) Generating plots 

#Volcano

par(las=1,bty='l')
plot(as.numeric(Robust_correlation),Corrected_p_value,
     cex=as.numeric(N_shared_observation)/5,pch=21,
     bg=string.to.colors(Corrected_p_value>1.3,colors = c("grey","red3")),
     xlim=c(-1,1),xaxs='i',ylim=c(0,4),yaxs='i',xlab='Robust correlation coefficient',cex.lab=1.3,
     ylab="-Log10(fdr")
 
abline(h=2,lty=2,col="red3")
plot(as.numeric(Robust_correlation),Corrected_p_value,
     cex=as.numeric(N_shared_observation)/3,pch=21,
     bg=string.to.colors(Corrected_p_value>2,colors = c("grey","red3")),
     xlim=c(-1.05,1.05),xaxs='i',ylim=c(0,4),yaxs='i',xlab='Robust correlation coefficient',cex.lab=1.3,
     ylab="-Log10(FDR)")
abline(h=2,lty=2,col="red3")
List_interaction[lower.tri(Robust_correlation)]=NA
text(as.numeric(Robust_correlation)[Corrected_p_value>2],Corrected_p_value[Corrected_p_value>2],as.character(List_interaction)[Corrected_p_value>2])


View(data.frame(as.numeric(Robust_correlation),-log10(p.adjust(as.numeric(Robust_correlation_pvalue),method = 'fdr')),as.character(List_interaction)))

par(las=1,bty='l')
plot(Table_fit_reshape_tau$`Lymphatic endothelial cells`,Table_fit_reshape_tau$`CD8+ T cells`,
     xlim=c(0,400),ylim=c(0,750),xaxs='i',yaxs='i',log='',
     ylab="tau parameter in CD8+ T cells (µm)",
     xlab="tau parameter in lymphatic endothelial cells (µm)",cex=2,bg="red3",pch=21,cex.lab=1.3)
m = lm(Table_fit_reshape_tau$`CD8+ T cells`~Table_fit_reshape_tau$`Lymphatic endothelial cells`)
abline(coef(m),lty=2,lwd=2)


#VII)Spatial plots

Plot_sample_clustering = function(ImageNumber="0028976_Region_1",cluster = "Lymphatic endothelial cells") {
  par(las=1,bty='l')
  sce_temp = sce[,sce$ImageNumber==ImageNumber]
  
  x_range = range(sce_temp$Location_Center_X)
  y_range = range(sce_temp$Location_Center_Y)
  Delta_x = x_range[2]-x_range[1]
  Delta_y = y_range[2]-y_range[1]
  Delta_max = max(Delta_x,Delta_y)

  
  plot(sce_temp$Location_Center_X,sce_temp$Location_Center_Y,
       col=string.to.colors(colLabels(sce_temp)==cluster,colors = c("grey","red3")),
       cex=ifelse(colLabels(sce_temp)==cluster,yes = 0.35,no = 0.1),xaxs='i',yaxs='i',
       xlim = c(x_range[1],x_range[1]+Delta_max),
       ylim = c(y_range[1],y_range[1]+Delta_max),
       xlab="X location (µm)",ylab="Y location (µm)",cex.lab=1.3,pch=16,main = paste(ImageNumber,"Cluster",cluster))
}

Plot_two_cell_types = function(ImageNumber="0028976_Region_1",
                               Cell_type_2 = "Lymphatic endothelial cells",
                               Cell_type_1 = "CD8+ T cells") {
  Temp_sce = sce[,sce$ImageNumber==ImageNumber]
  
  x_range = range(Temp_sce$Location_Center_X)
  y_range = range(Temp_sce$Location_Center_Y)
  Delta_x = x_range[2]-x_range[1]
  Delta_y = y_range[2]-y_range[1]
  Delta_max = max(Delta_x,Delta_y)
  
  par(las=1,bty='l')
  plot(Temp_sce$Location_Center_X,Temp_sce$Location_Center_Y,
       xlim = c(x_range[1],x_range[1]+Delta_max),
       ylim = c(y_range[1],y_range[1]+Delta_max),
       cex=0.1,pch=16,col="grey",
       xlab="X location (µm)",ylab="Y location (µm)",xaxs='i',yaxs='i',main=ImageNumber)
  points(Temp_sce$Location_Center_X[colLabels(Temp_sce)==Cell_type_1],
         Temp_sce$Location_Center_Y[colLabels(Temp_sce)==Cell_type_1],cex=0.1,col="orange")
  points(Temp_sce$Location_Center_X[colLabels(Temp_sce)==Cell_type_2],
         Temp_sce$Location_Center_Y[colLabels(Temp_sce)==Cell_type_2],cex=0.1,col="darkblue")
  
}

Plot_two_cell_types("0029448_Region_2")
Plot_two_cell_types("0029448_Region_1")

Plot_sample_clustering(cluster = "Inflammed hepatocytes",ImageNumber = "0029448_Region_3")

Plot_sample_clustering(cluster = "Lymphatic endothelial cells",ImageNumber = "0028976_Region_2")


Plot_sample_clustering(cluster = "CD8+ T cells",ImageNumber = "0029448_Region_1")


Plot_sample_clustering(cluster = "Plasma cells",ImageNumber = "0029448_Region_4")

Plot_sample_clustering(cluster = "Kupfer cells",ImageNumber = "0028976_Region_1")

Plot_two_cell_types("0029448_Region_2")

Plot_two_cell_types("0029448_Region_1")

pot(Table_fit_reshape_p$`Lymphatic endothelial cells`,Table_fit_reshape_p$`CD8+ T cells`,xlim=c(0,3),ylim=c(0,3))

