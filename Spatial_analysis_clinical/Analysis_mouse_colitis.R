library(pheatmap)
library(fifer)
library(FactoMineR)

#I)Data and function loading

source("Desktop/Projet_analyse_pcf/List_function_pcf_analysis.R")

Cell_data = read.delim("cell_properties.csv",sep=",")

View(x)

Plot_mice_sample = function(sample=1) {
  sample_ID = unique(Cell_data$Slice_ID)[sample]
  x_temp = Cell_data[Cell_data$Slice_ID==sample_ID,]
  plot(x_temp$x,x_temp$y,col=string.to.colors(x_temp$Tier2),cex=0.1)
}

Plot_mice_sample_clustering = function(sample="062921_D0_m3a_1_slice_2",cluster = "Arterial EC",width_FoV=2500) {
  par(las=1,bty='l')
  x_temp = Cell_data[Cell_data$Slice_ID==sample,]
  x_temp$x = x_temp$x-min(x_temp$x)
  x_temp$y = x_temp$y-min(x_temp$y)
  
  plot(x_temp$x,x_temp$y,col=string.to.colors(x_temp$Tier2==cluster,colors = c("grey","red3")),
       cex=ifelse(x_temp$Tier2==cluster,yes = 0.5,no = 0.2),xaxs='i',yaxs='i',
       xlim = c(0,width_FoV),ylim = c(0,width_FoV),
       xlab="X location (µm)",ylab="Y location (µm)",cex.lab=1.3,pch=16,main = paste(sample,cluster))
}

Plot_mice_sample(25)
Plot_mice_sample_clustering(cluster = "B cell 1")

#II) Compositional analysis

#A) Simple analysis
Table_composition = table(factor(Cell_data$Slice_ID,levels = unique(Cell_data$Slice_ID)),Cell_data$Tier2)
CA_anaysis = CA(Table_composition)

Total_cell_number = rowSums(Table_composition)
l = rownames(Table_composition)
l = strsplit(l,split = "_",fixed = TRUE)
l = unlist(lapply(l,FUN = function(x) {x[2]}))

P_value_DA = c()
for (k in colnames(Table_composition)) {
  m_1 = glm(Table_composition[,k] ~ 1 + offset(log(Total_cell_number)) + factor(l),family = "poisson")
  m_0 = glm(Table_composition[,k] ~ 1 + offset(log(Total_cell_number)),family = "poisson")
  u = anova(m_0,m_1)
  P_value_temp =pchisq(u$Deviance[2],(u$Df[2]), lower.tail = FALSE,log.p = F)
  P_value_DA = c(P_value_DA,P_value_temp)
}

P_value_DA = p.adjust(P_value_DA,method = "BH")
names(P_value_DA) = colnames(Table_composition)

#Everything is significant ! 


#III)Looking at the spatial structure of individual cell types


#A) Computing pcf for each cell type

r_vector = seq(0,500,length.out = 100)
N_min_points = 100
List_table_pcf = c()
Clustering = Cell_data$Tier2
for (k in unique(Cell_data$Slice_ID)) {
  print(k)
  Table_pcf=matrix(1,nrow = length(r_vector),ncol = length(unique(Clustering)))
  colnames(Table_pcf) = unique(Clustering)
  for (j in unique(Clustering)) {
    print(j)
    Selected_cells = Cell_data$Slice_ID==k & Clustering==j
    ppp_temp = ppp(Cell_data$x[Selected_cells],
                   Cell_data$y[Selected_cells],
                   window = owin(xrange = range(Cell_data$x[Cell_data$Slice_ID==k]),
                                 yrange = range(Cell_data$y[Cell_data$Slice_ID==k])))
    
    if (sum(Selected_cells)>N_min_points) {
      Kest_temp = Kest(ppp_temp,r = r_vector,rmax = max(r_vector),correction = "isotropic")
      pcf_temp = pcf.fv(Kest_temp,method = "b") #Pcf is computed by taking the derivative of the K function and by setting pcf(0)=0 as two cells cannot be at the same place 
      Table_pcf[,j]=pcf_temp$pcf
    }
  }
  List_table_pcf[[k]] = Table_pcf
}



#B)Fitting the sigmoid model


Table_fits = c()
for (k in 1:length(List_table_pcf)) {
  for (j in 1:ncol(List_table_pcf[[k]])) {
    temp_fit = try(Fit_sigmoid(x = r_vector,y = List_table_pcf[[k]][,j],show_plot = FALSE))
    if (class(temp_fit)=="try-error") {
      temp_fit = rep(NA,5)
    }
    Table_fits = rbind(Table_fits,c(temp_fit,k,j))
  }
}
colnames(Table_fits) = c("tau","p", "C","C_normalised","R2","Sample","Cluster")
Table_fits= as.data.frame(Table_fits)

U = data.frame(Parameter = Table_fits$tau,Cluster =Table_fits$Cluster )
U = U[is.finite(U$Parameter),]
prism.plots(Parameter~Cluster,U,pch=21,bg='red3',cex=1.5)

Table_fits$Sample = unique(Cell_data$Slice_ID)[Table_fits$Sample]
unique(Clustering)[35]


x = strsplit(Table_fits$Sample,split = "_",fixed = TRUE)
x = unlist(lapply(x,FUN = function(x){x[2]}))
Table_fits$Condition = x

#C) Computing mean behavior of each cell type 

Mean_cluster_tau = aggregate(Table_fits$tau,by=list(Table_fits$Cluster),FUN =function(x) {median(x[is.finite(x)])})
Mean_cluster_tau$Clusters = unique(Clustering)[Mean_cluster_tau$Group.1]

Mean_cluster_C = aggregate(Table_fits$C_normalised,by=list(Table_fits$Cluster),FUN =function(x) {median(x[is.finite(x)])})
Mean_cluster_C$Clusters = unique(Clustering)[Mean_cluster_C$Group.1]

Mean_cluster_p = aggregate(Table_fits$p,by=list(Table_fits$Cluster),FUN =function(x) {median(x[is.finite(x)])})
Mean_cluster_p$Clusters = unique(Clustering)[Mean_cluster_p$Group.1]

N_cell_cluster = table(Clustering)
N_cell_cluster = N_cell_cluster[Mean_cluster_p$Clusters]
N_cell_cluster = N_cell_cluster/sum(N_cell_cluster)*100

#D) Plots

par(las=1,bty='l')
plot(Mean_cluster_tau$x,Mean_cluster_p$x,cex=sqrt(N_cell_cluster)*2,
     xlab="Tau parameter (µm)",ylab="p parameter",pch=21,log='',
     bg="red3",cex.lab=1.3)
plot(Mean_cluster_tau$x,Mean_cluster_p$x,cex=sqrt(N_cell_cluster)*2,
     xlab="Tau parameter (µm)",ylab="p parameter",pch=21,log='',
     bg="red3",cex.lab=1.3)
text(Mean_cluster_tau$x,Mean_cluster_p$x,Mean_cluster_p$Clusters)

x = Mean_cluster_tau
x = x[order(x$x,decreasing = TRUE),]
x = x[1:15,]
par(las=1,bty='l',mar=c(4.5,12,2,2))
barplot(x$x,horiz = TRUE,xlim=c(0,max(x$x*1.2)),
        names.arg = x$Clusters,xlab="Tau parameter (µm)",cex.lab=1.3)

x = Mean_cluster_p
x = x[order(x$x,decreasing = TRUE),]
x = x[1:15,]
par(las=1,bty='l',mar=c(4.5,12,2,2))
barplot(x$x,horiz = TRUE,xlim=c(0,max(x$x*1.2)),
        names.arg = x$Clusters,xlab="Tau parameter (µm)",cex.lab=1.3)

Plot_mice_sample_clustering(sample = "062921_D9_m5_1_slice_1",cluster = "B cell 1")
Plot_mice_sample_clustering(sample = "082421_D0_m7_1_slice_2",cluster = "B cell 1")
Plot_mice_sample_clustering(sample = "100221_D9_m5_2_slice_2",cluster = "B cell 1")

Plot_mice_sample_clustering(sample = "062921_D9_m5_1_slice_2",cluster = "Macrophage (Cxcl10+)")
Plot_mice_sample_clustering(sample = "062221_D9_m3_2_slice_3",cluster = "Macrophage (Cxcl10+)")
Plot_mice_sample_clustering(sample = "100221_D9_m5_2_slice_1",cluster = "Macrophage (Cxcl10+)")


Plot_mice_sample_clustering(sample = "082421_D0_m7_1_slice_1",cluster = "Stem cells")
Plot_mice_sample_clustering(sample = "092421_D3_m2_1_slice_1",cluster = "Stem cells")
Plot_mice_sample_clustering(sample = "100221_D9_m3_2_slice_2",cluster = "Stem cells")


Plot_mice_sample_clustering(sample = "062921_D0_m3a_2_slice_1",cluster = "SMC 2")
Plot_mice_sample_clustering(sample = "062921_D9_m5_2_slice_1",cluster = "SMC 2")



Plot_mice_sample_clustering(cluster = "Fibro 1")


#IV) Kinetic analysis


#Quick function to check the stats are good 

Kinetic_cluster_tau = function(Selected_cluster = "SMC 1") {
  Selected_cluster_number = which(unique(Clustering)==Selected_cluster)
  x = Table_fits[Table_fits$Cluster==Selected_cluster_number,]
  x = x[!is.na(x$tau),]
  x$Condition = factor(x$Condition,c("D0","D3","D9","D21"))
  x = x[order(x$Condition),]
  plot(NULL,xlim=c(0.5,4.5),ylim=c(0,max(x$tau)*1.1),xaxs="i",yaxs='i',xaxt='n',
       xlab="Condition",ylab="Tau parameter (µm)",cex.lab=1.3,main=Selected_cluster)
  prism.plots(tau~Condition,x,pch=21,bg=string.to.colors(x$Condition,colors = c("grey","orange2","red2","darkred")),
              cex=1.5,col="black",add=TRUE)
}

Kinetic_cluster_p = function(Selected_cluster = "SMC 1") {
  Selected_cluster_number = which(unique(Clustering)==Selected_cluster)
  x = Table_fits[Table_fits$Cluster==Selected_cluster_number,]
  x = x[!is.na(x$tau),]
  x$Condition = factor(x$Condition,c("D0","D3","D9","D21"))
  x = x[order(x$Condition),]
  plot(NULL,xlim=c(0.5,4.5),ylim=c(min(x$p)*0.9,max(x$p)*1.1),xaxs="i",yaxs='i',xaxt='n',
       xlab="Condition",ylab="p parameter",cex.lab=1.3,main=Selected_cluster)
  prism.plots(p~Condition,x,pch=21,bg=string.to.colors(x$Condition,colors = c("grey","orange2","red2","darkred")),
              cex=1.5,col="black",add=TRUE)
}

Kinetic_cluster_C = function(Selected_cluster = "SMC 1") {
  Selected_cluster_number = which(unique(Clustering)==Selected_cluster)
  x = Table_fits[Table_fits$Cluster==Selected_cluster_number,]
  x = x[!is.na(x$tau),]
  x$Condition = factor(x$Condition,c("D0","D3","D9","D21"))
  x = x[order(x$Condition),]
  plot(NULL,xlim=c(0.5,4.5),ylim=c(0,max(x$C_normalised)*1.1),xaxs="i",yaxs='i',xaxt='n',
       xlab="Condition",ylab="C parameter",cex.lab=1.3,main=Selected_cluster)
  prism.plots(C_normalised~Condition,x,pch=21,bg=string.to.colors(x$Condition,colors = c("grey","orange2","red2","darkred")),
              cex=1.5,col="black",add=TRUE)
}



#A) Kinetic behaviour for tau

List_p_values = c()
for (k in unique(Clustering)) {
  Selected_cluster = k
  Selected_cluster = which(unique(Clustering)==Selected_cluster)
  x = Table_fits[Table_fits$Cluster==Selected_cluster,]
  x = x[!is.na(x$C_normalised),]
  p_temp = NA
  if (length(unique(x$Condition))>1 & nrow(x)>=20) {
    m = anova(aov(x$tau~x$Condition))
    p_temp = m$`Pr(>F)`[1]
  }
  List_p_values = c(List_p_values,p_temp)
}
names(List_p_values) = unique(Clustering)
List_p_values = List_p_values[!is.na(List_p_values)]
List_p_values_tau = p.adjust(List_p_values,method = "BH")
View(data.frame(List_p_values_tau))

#B) Kinetic behaviour for p

List_p_values = c()
for (k in unique(Clustering)) {
  Selected_cluster = k
  Selected_cluster = which(unique(Clustering)==Selected_cluster)
  x = Table_fits[Table_fits$Cluster==Selected_cluster,]
  x = x[!is.na(x$C_normalised),]
  p_temp = NA
  if (length(unique(x$Condition))>1 & nrow(x)>=20) {
    m = anova(aov(x$p~x$Condition))
    p_temp = m$`Pr(>F)`[1]
  }
  List_p_values = c(List_p_values,p_temp)
}
names(List_p_values) = unique(Clustering)
List_p_values = List_p_values[!is.na(List_p_values)]
List_p_values_p = p.adjust(List_p_values,method = "BH")
View(data.frame(List_p_values_p))


#D) Output Plots 

l = unique(c(names(List_p_values_p),names(List_p_values_tau)))
par(las=1,bty='l')
plot(-log10(List_p_values_tau[l]),-log10(List_p_values_p[l]),
     xlab="-Log10(fdr) for tau",ylab="-Log10(fdr) for p",
     xaxs='i',yaxs='i',ylim=c(0,6.5),xlim=c(0,5),cex=2,pch=21,bg='red3')
abline(v=2,lty=2,col="grey")
abline(h=2,lty=2,col="grey")
plot(-log10(List_p_values_tau[l]),-log10(List_p_values_p[l]),
     xlab="-Log10(fdr) for tau",ylab="-Log10(fdr) for p",
     xaxs='i',yaxs='i',ylim=c(0,6.5),xlim=c(0,5),cex=2,pch=21,bg='red3')
abline(v=2,lty=2,col="grey")
abline(h=2,lty=2,col="grey")
text(-log10(List_p_values_tau[l]),-log10(List_p_values_p[l]),labels = l)


par(las=1,bty='l')
Kinetic_cluster_tau("SMC 1")
Kinetic_cluster_p("SMC 1")


par(las=1,bty='l')
Kinetic_cluster_tau("Stem cells")
Kinetic_cluster_p("Stem cells")


par(las=1,bty='l')
Kinetic_cluster_tau("SMC 2")
Kinetic_cluster_p("SMC 2")


par(las=1,bty='l')
Kinetic_cluster_tau("Fibro 1")
Kinetic_cluster_p("Fibro 1")


Plot_mice_sample_clustering(sample = "062921_D0_m3a_1_slice_1",cluster = "SMC 1")
Plot_mice_sample_clustering(sample = "092421_D3_m3_1_slice_3",cluster = "SMC 1")
Plot_mice_sample_clustering(sample = "100221_D9_m5_2_slice_3",cluster = "SMC 1")
Plot_mice_sample_clustering(sample = "062921_D9_m5_1_slice_2",cluster = "Stem cells")
Plot_mice_sample_clustering(sample = "062921_D0_m3a_2_slice_2",cluster = "Stem cells")
Plot_mice_sample_clustering(sample = "092421_D3_m4_1_slice_2",cluster = "Stem cells")
Plot_mice_sample_clustering(sample = "100221_D9_m5_2_slice_2",cluster = "IAF 3")


#V)Looking at the spatial structure of multipke cell types

#A)Reshaping tau, p and C normalised

Table_fit_reshape_tau = c()
for (k in unique(Table_fits$Sample)) {
  x = Table_fits[Table_fits$Sample==k,]
  Table_fit_reshape_tau = rbind(Table_fit_reshape_tau,x$tau)
}

colnames(Table_fit_reshape_tau) = unique(Clustering)
Table_fit_reshape_tau = as.data.frame(Table_fit_reshape_tau)
Table_fit_reshape_tau = Table_fit_reshape_tau[,colSums(!is.na(Table_fit_reshape_tau))>=20]


Table_fit_reshape_p = c()
for (k in unique(Table_fits$Sample)) {
  x = Table_fits[Table_fits$Sample==k,]
  Table_fit_reshape_p = rbind(Table_fit_reshape_p,x$p)
}

colnames(Table_fit_reshape_p) = unique(Clustering)
Table_fit_reshape_p = as.data.frame(Table_fit_reshape_p)
Table_fit_reshape_p = Table_fit_reshape_p[,colSums(!is.na(Table_fit_reshape_p))>=20]


Table_fit_reshape_C = c()
for (k in unique(Table_fits$Sample)) {
  x = Table_fits[Table_fits$Sample==k,]
  Table_fit_reshape_C = rbind(Table_fit_reshape_C,x$C_normalised)
}

colnames(Table_fit_reshape_C) = unique(Clustering)
Table_fit_reshape_C = as.data.frame(Table_fit_reshape_C)
Table_fit_reshape_C = Table_fit_reshape_C[,colSums(!is.na(Table_fit_reshape_C))>=20]


Reshaped_parameters = c()
for (k in unique(Table_fits$Sample)) {
  x = Table_fits[Table_fits$Sample==k,]
  Reshaped_parameters = rbind(Reshaped_parameters,x$Condition)
}
Reshaped_parameters = apply(Reshaped_parameters,MARGIN = 1,FUN = unique)
#B) Correlation analysis 

pheatmap(cor(Table_fit_reshape_tau,use = "pairwise.complete.obs",method = 'pearson'),clustering_method = 'ward')
pheatmap(cor(Table_fit_reshape_p,use = "pairwise.complete.obs",method = 'pearson'),clustering_method = 'ward')
pheatmap(cor(Table_fit_reshape_C,use = "pairwise.complete.obs",method = 'pearson'),clustering_method = 'ward')

par(las=1,bty='l',mar=c(4,4,4,4))
plot(Table_fit_reshape_tau$`SMC 1`,Table_fit_reshape_tau$`SMC 2`,log='')
m = lm(Table_fit_reshape_tau$`SMC 2`~Table_fit_reshape_tau$`SMC 1`)
abline(coef(m))

library(WRS2)
m = pbcor(Table_fit_reshape_tau[,1],Table_fit_reshape_tau[,3])

plot(Table_fit_reshape_tau$`IAE 2`,Table_fit_reshape_tau$`IAE 3`,log='')


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
     cex=as.numeric(N_shared_observation)/15,pch=21,
     bg=string.to.colors(Corrected_p_value>2,colors = c("grey","red3")),
     xlim=c(-1,1),xaxs='i',ylim=c(0,8),yaxs='i',xlab='Robust correlation coefficient',cex.lab=1.3,
     ylab="-Log10(fdr")
abline(h=2,lty=2,col="red3")
plot(as.numeric(Robust_correlation),Corrected_p_value,
     cex=as.numeric(N_shared_observation)/15,pch=21,
     bg=string.to.colors(Corrected_p_value>2,colors = c("grey","red3")),
     xlim=c(-1,1),xaxs='i',ylim=c(0,8),yaxs='i',xlab='Robust correlation coefficient',cex.lab=1.3,
     ylab="-Log10(fdr")
abline(h=2,lty=2,col="red3")

List_interaction[lower.tri(Robust_correlation)]=NA
text(as.numeric(Robust_correlation)[Corrected_p_value>2],Corrected_p_value[Corrected_p_value>2],as.character(List_interaction)[Corrected_p_value>2])
View(data.frame(as.numeric(Robust_correlation),-log10(p.adjust(as.numeric(Robust_correlation_pvalue),method = 'fdr')),as.character(List_interaction)))


Plot_correlation_tau = function(Cell_type_1 = "SMC 1",Cell_type_2 = "SMC 2") {
  plot(Table_fit_reshape_tau[,Cell_type_1],Table_fit_reshape_tau[,Cell_type_2],xaxs='i',yaxs='i',
       xlim=c(0,max(Table_fit_reshape_tau[,Cell_type_1],na.rm = T)*1.1),
       ylim=c(0,max(Table_fit_reshape_tau[,Cell_type_2],na.rm = T)*1.1),
       pch=21,bg=string.to.colors(Reshaped_parameters,colors = c("grey","orange2","red2","darkred")),cex=2,
       xlab= paste("tau parameter in",Cell_type_1,"(µm)"),
       ylab= paste("tau parameter in",Cell_type_2,"(µm)"),cex.lab=1.3)
  m = lm (Table_fit_reshape_tau[,Cell_type_2]~Table_fit_reshape_tau[,Cell_type_1])
  abline(coef(m),lwd=2,lty=2,col="black")
}

par(las=1,bty='l')
Plot_correlation_tau("SMC 1","Lymphatic EC (Ccl21a+)")
Plot_correlation_tau("SMC 1","IAF 2")
Plot_correlation_tau("B cell 1","B cell 2")
Plot_correlation_tau("SMC 1","SMC 2")
Plot_correlation_tau("IAF 2","IAF 3")
Plot_correlation_tau("IAF 2","Neutrophil 1")

#D) 2D plots 

Plot_mice_sample_two_clustering = function(sample="062221_D9_m3_2_slice_2",cluster_1 = "B cell 1",
                                           cluster_2 = "B cell 2",width_FoV=2500) {
  par(las=1,bty='l')
  x_temp = Cell_data[Cell_data$Slice_ID==sample,]
  x_temp$x = x_temp$x-min(x_temp$x)
  x_temp$y = x_temp$y-min(x_temp$y)
  
  plot(x_temp$x,x_temp$y,col="grey",
       cex = 0.2,xaxs='i',yaxs='i',
       xlim = c(0,width_FoV),ylim = c(0,width_FoV),
       xlab="X location (µm)",ylab="Y location (µm)",cex.lab=1.3,pch=16,main = paste(sample,cluster_1,cluster_2))
  points(x_temp$x[x_temp$Tier2==cluster_1],x_temp$y[x_temp$Tier2==cluster_1],pch=16,cex=0.3,col="orange")
  points(x_temp$x[x_temp$Tier2==cluster_2],x_temp$y[x_temp$Tier2==cluster_2],pch=16,cex=0.3,col="darkblue")
}

Plot_mice_sample_two_clustering(sample="062221_D9_m3_2_slice_2",cluster_1 = "B cell 1",cluster_2 = "B cell 2")
Plot_mice_sample_two_clustering(sample="062921_D9_m2a_1_slice_1",cluster_1 = "IAF 2",cluster_2 = "IAF 3")
Plot_mice_sample_two_clustering(sample="062921_D0_m3a_1_slice_1",cluster_1 = "SMC 1",cluster_2 = "SMC 2")
Plot_mice_sample_two_clustering(sample="100221_D9_m5_2_slice_1",cluster_1 = "Neutrophil 1",cluster_2 = "IAF 2")
Plot_mice_sample_two_clustering(sample="062921_D9_m2a_2_slice_2",cluster_1 = "Neutrophil 1",cluster_2 = "IAF 2")
Plot_mice_sample_two_clustering(sample="082421_D21_m1_1_slice_2",cluster_1 = "Neutrophil 1",cluster_2 = "IAF 2")




#V)Checking the overlap with cellular abundance

#A)Correlation between cell abundance and cell parameters

Table_composition = as.data.frame.matrix(Table_composition) 
#Table_composition = Table_composition/rowSums(Table_composition)

Vector_correlation_p = c()
Vector_correlation_tau = c()
Vector_correlation_C = c()

for (k in colnames(Table_fit_reshape_tau)) {
  Vector_correlation_p = c(Vector_correlation_p,cor(Table_composition[,k],Table_fit_reshape_p[,k],use = 'pairwise.complete.obs'))
  Vector_correlation_tau = c(Vector_correlation_tau,cor(Table_composition[,k],Table_fit_reshape_tau[,k],use = 'pairwise.complete.obs'))
  Vector_correlation_C = c(Vector_correlation_C,cor(Table_composition[,k],Table_fit_reshape_C[,k],use = 'pairwise.complete.obs'))
}
plot(Table_composition$`Neutrophil 2`,Table_fit_reshape_tau$`Neutrophil 2`)

par(las=1,bty='l')
boxplot(Vector_correlation_tau,Vector_correlation_p,Vector_correlation_C,ylim=c(-1,1),yaxs='i',
        names=c("Tau","p","C"),ylab="Correlation with cell abundance",cex.lab=1.3,xlab="Parameter")


#B) Comparing the results of differential abundance and differential structure test

x = rownames(Table_composition)
x = strsplit(x,split = "_",fixed = T)
x = unlist(lapply(x,FUN = function(x) {x[2]}))
Time_variable = factor(x)
Total_cell_number = rowSums(Table_composition)

P_value_abundance = c()

for (k in colnames(Table_composition)) {
  x = Table_composition[,k]
  m = glm(x~1+offset(log(Total_cell_number))+Time_variable,family = poisson(link = "log"))
  m =anova(m,test = 'LRT')  
  P_value_abundance = c(P_value_abundance,m$`Pr(>Chi)`[2])
}

names(P_value_abundance) = colnames(Table_composition)
P_value_abundance = p.adjust(P_value_abundance,method = "BH")
View(data.frame(P_value_abundance))

Plot_abundance_condition = function(Selected_cluster="SMC 1") {
  par(las=1,bty='l')
  Condition_temp = factor(Time_variable,c("D0","D3","D9","D21"))
  U = data.frame(Abundance = Table_composition[,Selected_cluster]/Total_cell_number*100,Condition = Condition_temp)
  U = U[order(U$Condition),]
  plot(NULL,xlim=c(0.5,4.5),ylim=c(0,max(U$Abundance)*1.1),xaxs="i",yaxs='i',xaxt='n',
       xlab="Condition",ylab="Cell type proportion (%)",cex.lab=1.3,main=Selected_cluster)
  prism.plots(Abundance~Condition,U,pch=21,bg=string.to.colors(U$Condition,colors = c("grey","orange2","red2","darkred")),
              cex=1.5,col="black",add=TRUE)

}

Plot_abundance_condition("Colonocytes")
Kinetic_cluster_tau("Colonocytes")
Plot_abundance_condition("Pericyte 1")
Kinetic_cluster_tau("Pericyte 1")



Plot_mice_sample_clustering(sample = "062921_D9_m5_2_slice_1","Colonocytes")
Plot_mice_sample_clustering(sample = "082421_D0_m7_1_slice_1","Colonocytes")
