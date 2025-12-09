library(RColorBrewer)
library(pheatmap)

#I)Loading function/object

#A)List functions

source("Desktop/Projet_analyse_pcf/List_function_pcf_analysis.R")

#B)Sce object
#Reading the SCE object 
sce= readRDS("/Thyroid_sce_object.rds")

{Renamed_clustering = List_pcf$Annotation$Cluster
  Renamed_clustering[Renamed_clustering==24]="B_cell_H2AX"
  Renamed_clustering[Renamed_clustering==23]="B_cell_CD45RA"
  Renamed_clustering[Renamed_clustering==4]="B_cell_activated"
  Renamed_clustering[Renamed_clustering==19]="Tfh"
  Renamed_clustering[Renamed_clustering==15]="pDC"
  Renamed_clustering[Renamed_clustering==21]="CD8_T_cells_activated"
  Renamed_clustering[Renamed_clustering==12]="CD8_T_cells_resting"
  Renamed_clustering[Renamed_clustering==13]="Naive_CD4_T_cells"
  Renamed_clustering[Renamed_clustering==20]="FC_1"
  Renamed_clustering[Renamed_clustering==8]="FC_2"
  Renamed_clustering[Renamed_clustering==14]="FC_3"
  Renamed_clustering[Renamed_clustering==16]="Fibronecting_cells"
  Renamed_clustering[Renamed_clustering==3]="Segmentation_errors"
  Renamed_clustering[Renamed_clustering==22]="migDCs"
  Renamed_clustering[Renamed_clustering==9]="Vessels"
  Renamed_clustering[Renamed_clustering==2]="FC_4"
  Renamed_clustering[Renamed_clustering==11]="Stromal_cells"
  Renamed_clustering[Renamed_clustering%in%c(1,5)]="Unassigned cells"
  Renamed_clustering[Renamed_clustering==6]="Plasma_cells_MX1"
  Renamed_clustering[Renamed_clustering==10]="Plasma_cells"
  Renamed_clustering[Renamed_clustering==18]="Neutrophils"
  Renamed_clustering[Renamed_clustering==7]="GSN_cells"
  Renamed_clustering[Renamed_clustering==17]="Macrophages"
  List_pcf$Annotation$Cluster_renamed = Renamed_clustering
}

x = colLabels(sce)
{
  x[x==24]="B_cell_H2AX"
  x[x==23]="B_cell_CD45RA"
  x[x==4]="B_cell_activated"
  x[x==19]="Tfh"
  x[x==15]="pDC"
  x[x==21]="CD8_T_cells_activated"
  x[x==12]="CD8_T_cells_resting"
  x[x==13]="Naive_CD4_T_cells"
  x[x==20]="FC_1"
  x[x==8]="FC_2"
  x[x==14]="FC_3"
  x[x==16]="Fibronecting_cells"
  x[x==3]="Segmentation_errors"
  x[x==22]="migDCs"
  x[x==9]="Vessels"
  x[x==2]="FC_4"
  x[x==11]="Stromal_cells"
  x[x%in%c(1,5)]="Unassigned cells"
  x[x==6]="Plasma_cells_MX1"
  x[x==10]="Plasma_cells"
  x[x==18]="Neutrophils"
  x[x==7]="GSN_cells"
  x[x==17]="Macrophages"
}


#II)Computation and fitting of the pcf  

#A)Computing pcf

List_pcf = Compute_pcf(sce,r_vector = seq(0,500,length.out=50),computation_method = "direct")


#B)Checking the functions on an example
k = 26
x = List_pcf$List_r[[k]]
y = List_pcf$List_pcf[[k]]
plot(x,y)
Fit_sigmoid(x,y,show_plot = TRUE)


#C)Fitting the different models
Results_generalised_gamma = Fit_parametric_pcf_model(List_pcf,model = "Generalized_gamma")
Results_gamma= Fit_parametric_pcf_model(List_pcf,model = "Gamma")
Results_exponential = Fit_parametric_pcf_model(List_pcf,model = "Exponential")
Results_power_law = Fit_parametric_pcf_model(List_pcf,model = "Power_law")
Results_sigmoid= Fit_parametric_pcf_model(List_pcf,model = "Sigmoid")
Results_beta_prime= Fit_parametric_pcf_model(List_pcf,model = "Beta_prime")


#Checking the quaity of the fit
par(las=1,bty='l')
boxplot(Results_power_law$R2,Results_exponential$R2,Results_gamma$R2,
        Results_generalised_gamma$R2,Results_sigmoid$R2,Results_beta_prime$R2,ylim=c(0.,1),
        xaxs='i',yaxs='i',names = c("Power law","Exponential","Gamma","GG","Sigmoid","Beta prime"),
        xlab="Model",ylab="R2",cex.lab=1.3,col=brewer.pal(6, "Spectral"),cex.axis=0.7)


plot(Results_generalised_gamma$R2,Results_gamma$R2,xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i',
     xlab="Fitting quality of the Gamma model (R2)",
     ylab="Fitting quality of the Sigmoid model (R2)",pch=21,bg="red3",cex=1.5,cex.lab=1.5)
abline(0,1,lty=2,col="grey")


#III) Analysis of the average behavior of the different cell types    

#A)Filtering low quality fitting/Pcfs

Results_sigmoid$Cluster = List_pcf$Annotation$Cluster_renamed
Results_sigmoid$ROI = List_pcf$Annotation$ROI

Results_sigmoid_filtered = Results_sigmoid[Results_sigmoid$R2>0.6,]
Renamed_clustering_filtering = Renamed_clustering[Results_sigmoid$R2>0.6]

#B)Aggregating across cell types

Aggregated_tau_table = data.frame(Mean_value = aggregate(Results_sigmoid_filtered$tau,by=list(Renamed_clustering_filtering),FUN=mean),
                                  Sd_value = aggregate(Results_sigmoid_filtered$tau,by=list(Renamed_clustering_filtering),FUN=sd)[,2])

Aggregated_p_table = data.frame(Mean_value = aggregate(Results_sigmoid_filtered$p,by=list(Renamed_clustering_filtering),FUN=mean),
                                  Sd_value = aggregate(Results_sigmoid_filtered$p,by=list(Renamed_clustering_filtering),FUN=sd)[,2])

Aggregated_C_table = data.frame(Mean_value = aggregate(Results_sigmoid_filtered$C_normalised,by=list(Renamed_clustering_filtering),FUN=mean),
                                Sd_value = aggregate(Results_sigmoid_filtered$C_normalised,by=list(Renamed_clustering_filtering),FUN=sd)[,2])

#C)Resulting plots 

Ordered_cell_types = Aggregated_C_table$Mean_value.Group.1[order(Aggregated_C_table$Mean_value.x,decreasing = T)]
Ordered_cell_types = Ordered_cell_types[!Ordered_cell_types%in%c("Segmentation_errors","Unassigned")]

par(las=1,bty='l',mar=c(5,10,2,2))
boxplot(Results_sigmoid_filtered$C_normalised~factor(Renamed_clustering_filtering,levels = Ordered_cell_types),
        outline=F,horizontal=TRUE,xlab="Normalised C values",ylab="",yaxs='i',ylim=c(0,800))
dev.off()

N_cells = table(x)
N_cells = N_cells[Aggregated_tau_table$Mean_value.Group.1]
par(las=1,bty='l')
plot(Aggregated_tau_table$Mean_value.x,Aggregated_p_table$Mean_value.x,
     xlab="Mean tau value",ylab="Mean p value",cex.lab=1.3,pch=21,bg="red3",
     cex = sqrt(N_cells)/50 )
plot(Aggregated_tau_table$Mean_value.x,Aggregated_p_table$Mean_value.x,
     xlab="Mean tau value",ylab="Mean p value",cex.lab=1.3,pch=21,bg="red3",
     cex = sqrt(N_cells)/50 )
graphics::text(x = Aggregated_tau_table$Mean_value.x,y=Aggregated_p_table$Mean_value.x,labels=Aggregated_p_table$Mean_value.Group.1)


#IV)Comparison of the different conditions

#A) Reshapping of the meta-data
Meta_data = read.delim("Desktop/Thyroid_project/Final_acquisition/Sample_description.txt",row.names = 1)
Tissue_type = read.delim("Desktop/Thyroid_project/Final_acquisition/Sample_annotation.txt")
Tissue_type$Tissue_type[Tissue_type$Tissue_type=="q"] = "Test"
x = strsplit(Tissue_type$Sample,split = " ")
x = unlist(lapply(x,FUN = function(x) {x[1]}))

Tissue_type$Donor = x
Tissue_type$Status = Meta_data[Tissue_type$Donor,1]
Tissue_type$Merged_type = paste(Tissue_type$Status,Tissue_type$Tissue_type)
rownames(Tissue_type) = Tissue_type$Sample

Results_sigmoid_filtered$Merged_type = Tissue_type[Results_sigmoid_filtered$ROI,"Merged_type"]


#B) Comparing the spatial paraemter 

List_clusters = unique(Results_sigmoid_filtered$Cluster)
List_clusters = List_clusters[!is.na(List_clusters)]

List_p_values = c()

for (k in List_clusters) {
  
  x = Results_sigmoid_filtered[Results_sigmoid_filtered$Cluster==k,]
  
  p_temp = NA
  if (length(unique(x$Merged_type))>1 & nrow(x)>=20) {
    print(k)
    m = anova(aov(x$tau~x$Merged_type))
    p_temp = m$`Pr(>F)`[1]
  }
  List_p_values = c(List_p_values,p_temp)
}

names(List_p_values) = List_clusters
List_p_values = List_p_values[!is.na(List_p_values)]
List_p_values_tau = p.adjust(List_p_values,method = "BH")
View(data.frame(List_p_values_tau))


List_p_values = c()

for (k in List_clusters) {
  
  x = Results_sigmoid_filtered[Results_sigmoid_filtered$Cluster==k,]
  
  p_temp = NA
  if (length(unique(x$Merged_type))>1 & nrow(x)>=20) {
    print(k)
    m = anova(aov(x$p~x$Merged_type))
    p_temp = m$`Pr(>F)`[1]
  }
  List_p_values = c(List_p_values,p_temp)
}

names(List_p_values) = List_clusters
List_p_values = List_p_values[!is.na(List_p_values)]
List_p_values_p = p.adjust(List_p_values,method = "BH")
View(data.frame(List_p_values_p))



#C) Quick function to check the results

Plot_tau_condition = function(Cluster_temp = "FC_1") {
  x = Results_sigmoid_filtered[Results_sigmoid_filtered$Cluster==Cluster_temp,]
  x = x[!is.na(x$tau),]
  x = x[order(x$Merged_type),]
  if (length(unique(grepl(x$Merged_type,pattern = "Hashimoto")))>1) {
    plot(NULL,xlim=c(0.5,3.5),ylim=c(0,max(x$tau*1.1)),xaxs='i',yaxs='i',xaxt='n',ylab="Tau parameter value",xlab="",main=Cluster_temp)
    prism.plots(tau~Merged_type,data = x,add=TRUE,pch=ifelse(grepl(x$Merged_type,pattern = " T"),yes = 23,no = 21),
                cex=2,col="black",
                bg=string.to.colors(grepl(x$Merged_type,pattern = "Hashimoto"),colors = c("grey","orange")),cex.axis=0.7)
    
  }
  else {
    plot(NULL,xlim=c(0.5,2.5),ylim=c(0,max(x$tau*1.1)),xaxs='i',yaxs='i',xaxt='n',ylab="Tau parameter value",xlab="",main=Cluster_temp)
    prism.plots(tau~Merged_type,data = x,add=TRUE,pch=ifelse(grepl(x$Merged_type,pattern = " T"),yes = 23,no = 21),
                cex=2,col="black",
                bg="orange",cex.axis=0.7)
    
  }
}

Plot_p_condition = function(Cluster_temp = "FC_1") {
  x = Results_sigmoid_filtered[Results_sigmoid_filtered$Cluster==Cluster_temp,]
  x = x[!is.na(x$tau),]
  x = x[order(x$Merged_type),]
  plot(NULL,xlim=c(0.5,3.5),ylim=c(0,max(x$p*1.1)),xaxs='i',yaxs='i',xaxt='n',ylab="P parameter value",xlab="")
  prism.plots(p~Merged_type,x,pch=21,bg=string.to.colors(x$Merged_type),cex=2,col="black",add=TRUE)
}


l = unique(c(names(List_p_values_p),names(List_p_values_tau)))
pdf("Desktop/Papier_Bost_2/Plot_Figure_Thyroide//P_value_anova_tau_p.pdf",width = 6.5,height = 6.5,useDingbats = FALSE)
par(las=1,bty='l')
plot(-log10(List_p_values_tau[l]),-log10(List_p_values_p[l]),
     xlab="-Log10(fdr) for tau",ylab="-Log10(fdr) for p",
     xaxs='i',yaxs='i',ylim=c(0,4),xlim=c(0,4),cex=2,pch=21,bg='red3')
abline(v=2,lty=2,col="grey")
abline(h=2,lty=2,col="grey")
plot(-log10(List_p_values_tau[l]),-log10(List_p_values_p[l]),
     xlab="-Log10(fdr) for tau",ylab="-Log10(fdr) for p",
     xaxs='i',yaxs='i',ylim=c(0,4),xlim=c(0,4),cex=2,pch=21,bg='red3')
abline(v=2,lty=2,col="grey")
abline(h=2,lty=2,col="grey")
text(-log10(List_p_values_tau[l]),-log10(List_p_values_p[l]),labels = l)
dev.off()


par(las=1,bty='l')
Plot_tau_condition(Cluster_temp = "Tfh")
Plot_tau_condition(Cluster_temp = "Naive_CD4_T_cells")
Plot_tau_condition(Cluster_temp = "B_cell_activated")
Plot_tau_condition(Cluster_temp = "FC_4")

#V)Looking at the spatial structure of multipke cell types

#A) Reshaping tau, p and C normalised

Table_fit_reshape_tau = c()
for (k in unique(Results_sigmoid_filtered$ROI)) {
  x = Results_sigmoid_filtered[Results_sigmoid_filtered$ROI==k,]
  x = x[!is.na(x$tau),]
  x = x[x$Cluster!="Unassigned cells",]
  rownames(x) = x$Cluster
  x = x[unique(List_pcf$Annotation$Cluster_renamed),"tau"]
  Table_fit_reshape_tau = rbind(Table_fit_reshape_tau,x)
}
colnames(Table_fit_reshape_tau)=unique(List_pcf$Annotation$Cluster_renamed)

rownames(Table_fit_reshape_tau) = unique(Results_sigmoid_filtered$ROI)
Table_fit_reshape_tau = as.data.frame(Table_fit_reshape_tau)
Table_fit_reshape_tau = Table_fit_reshape_tau[,colSums(!is.na(Table_fit_reshape_tau))>=20]

#B) Computing the assocation score : Correlation and MI

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
     cex=as.numeric(N_shared_observation)/10,pch=21,
     bg=string.to.colors(Corrected_p_value>2,colors = c("grey","red3")),
     xlim=c(-1,1),xaxs='i',ylim=c(0,8),yaxs='i',xlab='Robust correlation coefficient',cex.lab=1.3,
     ylab="-Log10(fdr")
abline(h=2,lty=2,col="red3")
plot(as.numeric(Robust_correlation),Corrected_p_value,
     cex=as.numeric(N_shared_observation)/10,pch=21,
     bg=string.to.colors(Corrected_p_value>2,colors = c("grey","red3")),
     xlim=c(-1,1),xaxs='i',ylim=c(0,8),yaxs='i',xlab='Robust correlation coefficient',cex.lab=1.3,
     ylab="-Log10(fdr")
abline(h=2,lty=2,col="red3")

List_interaction[lower.tri(Robust_correlation)]=NA
text(as.numeric(Robust_correlation)[Corrected_p_value>2],Corrected_p_value[Corrected_p_value>2],as.character(List_interaction)[Corrected_p_value>2])

#Scatter plots 

par(las=1,bty='l')
plot(Table_fit_reshape_tau$Naive_CD4_T_cells,Table_fit_reshape_tau$B_cell_activated,
     xlim=c(0,300),ylim=c(0,250),xaxs='i',yaxs='i',
     xlab="Tau parameter for naive CD4+ T cells",
     ylab="Tau parameter for activated B cells",cex=2,cex.lab=1.3,
     bg= string.to.colors(grepl(Results_sigmoid_filtered$Merged_type,pattern = "Hashimoto"),colors = c("orange","grey")),
     pch=ifelse(grepl(Results_sigmoid_filtered$Merged_type,pattern = " T"),yes = 23,no = 21))
m = lm(Table_fit_reshape_tau$B_cell_activated~Table_fit_reshape_tau$Naive_CD4_T_cells)
abline(coef(m),lty=2,lwd=2)

plot(Table_fit_reshape_tau$Tfh,Table_fit_reshape_tau$FC_2,
     xlim=c(0,250),ylim=c(0,150),xaxs='i',yaxs='i',
     xlab="Tau parameter for Tfh cells",
     ylab="Tau parameter for FC 2",cex=2,cex.lab=1.3,
     bg= string.to.colors(grepl(Results_sigmoid_filtered$Merged_type,pattern = "Hashimoto"),colors = c("orange","grey")),
     pch=ifelse(grepl(Results_sigmoid_filtered$Merged_type,pattern = " T"),yes = 23,no = 21))
m = lm(Table_fit_reshape_tau$FC_2~Table_fit_reshape_tau$Tfh)
abline(coef(m),lty=2,lwd=2)

plot(Table_fit_reshape_tau$Naive_CD4_T_cells,Table_fit_reshape_tau$FC_2,
     xlim=c(0,200),ylim=c(0,150),xaxs='i',yaxs='i',
     xlab="Tau parameter for naive CD4+ T cells",
     ylab="Tau parameter for FC 2",cex=2,cex.lab=1.3,
     bg= string.to.colors(grepl(Results_sigmoid_filtered$Merged_type,pattern = "Hashimoto"),colors = c("orange","grey")),
     pch=ifelse(grepl(Results_sigmoid_filtered$Merged_type,pattern = " T"),yes = 23,no = 21))
m = lm(Table_fit_reshape_tau$FC_2~Table_fit_reshape_tau$Tfh)
abline(coef(m),lty=2,lwd=2)

#V) Various plots for figures

Plot_sample_clustering = function(ImageNumber="B21_22037 1",cluster = 1,width_FoV=1000) {
  par(las=1,bty='l')
  sce_temp = sce[,sce$ImageNumber==ImageNumber]
  plot(sce_temp$Location_Center_X,sce_temp$Location_Center_Y,
       col=string.to.colors(colLabels(sce_temp)==cluster,colors = c("grey","red3")),
       cex=ifelse(colLabels(sce_temp)==cluster,yes = 0.6,no = 0.35),xaxs='i',yaxs='i',
       xlim = c(0,width_FoV),ylim = c(0,width_FoV),
       xlab="X location (µm)",ylab="Y location (µm)",cex.lab=1.3,pch=16,main = paste(ImageNumber,"Cluster",cluster))
}

Plot_sample_two_clustering= function(ImageNumber="B21_22037 1",cluster_1 = 1,cluster_2 = 2,width_FoV=1000) {
  par(las=1,bty='l')
  sce_temp = sce[,sce$ImageNumber==ImageNumber]
  plot(sce_temp$Location_Center_X,sce_temp$Location_Center_Y,col="grey",
       cex = 0.35,xaxs='i',yaxs='i',
       xlim = c(0,width_FoV),ylim = c(0,width_FoV),
       xlab="X location (µm)",ylab="Y location (µm)",cex.lab=1.3,pch=16,main = paste(ImageNumber,": Clusters",cluster_1,"and",cluster_2))
  points(sce_temp$Location_Center_X[colLabels(sce_temp)==cluster_1],sce_temp$Location_Center_Y[colLabels(sce_temp)==cluster_1],pch=16,cex=0.5,col="orange")
  points(sce_temp$Location_Center_X[colLabels(sce_temp)==cluster_2],sce_temp$Location_Center_Y[colLabels(sce_temp)==cluster_2],pch=16,cex=0.5,col="darkblue")
}


Plot_sample_three_clustering= function(ImageNumber="B21_22037 1",cluster_1 = 1,cluster_2 = 2,cluster_3 = 3,width_FoV=1000) {
  par(las=1,bty='l')
  sce_temp = sce[,sce$ImageNumber==ImageNumber]
  plot(sce_temp$Location_Center_X,sce_temp$Location_Center_Y,col="grey",
       cex = 0.35,xaxs='i',yaxs='i',
       xlim = c(0,width_FoV),ylim = c(0,width_FoV),
       xlab="X location (µm)",ylab="Y location (µm)",cex.lab=1.3,pch=16,main = paste(ImageNumber,": Clusters",cluster_1,"and",cluster_2))
  points(sce_temp$Location_Center_X[colLabels(sce_temp)==cluster_1],sce_temp$Location_Center_Y[colLabels(sce_temp)==cluster_1],pch=16,cex=0.5,col="orange")
  points(sce_temp$Location_Center_X[colLabels(sce_temp)==cluster_2],sce_temp$Location_Center_Y[colLabels(sce_temp)==cluster_2],pch=16,cex=0.5,col="darkblue")
  points(sce_temp$Location_Center_X[colLabels(sce_temp)==cluster_3],sce_temp$Location_Center_Y[colLabels(sce_temp)==cluster_3],pch=16,cex=0.5,col="darkgreen")
}


# FC_2 : 8 ; Tfh : 19 ; 4 : B_cell_activated ;Naive_CD4_T_cells:  13 

#Naive T-cells and Activated B-cells
Plot_sample_two_clustering(ImageNumber = "B22_04111 2",cluster_1 = 4,cluster_2 = 13)
Plot_sample_two_clustering(ImageNumber = "B21_22037 3",cluster_1 = 4,cluster_2 = 13)

#FC 2 cells and naive T-cells

Plot_sample_three_clustering(ImageNumber = "B21_22037 4",cluster_1 = 8,cluster_2 = 13,cluster_3 = 19)
Plot_sample_three_clustering(ImageNumber = "B22_04111 1",cluster_1 = 8,cluster_2 = 13,cluster_3 = 19)





