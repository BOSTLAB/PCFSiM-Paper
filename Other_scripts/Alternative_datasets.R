library(spatstat.data)
library(RColorBrewer)
library(fifer)
source("Desktop/Projet_analyse_pcf/List_function_pcf_analysis.R")

#I)Creating a functions to compute pcfs fit sigmoid for ppp object

Get_pcf_sigmoid_fit = function(ppp_object,N_min_points=100,r_vector = seq(0,0.5,length.out = 100),initial_tau=0.1) {
  Table_pcf = c()
  if (is.null(marks(ppp_object))) {
    marks(ppp_object) = rep("No_mark",npoints(ppp_object))
  }
  List_marks = as.character(unique(ppp_object$marks))
  Colnames_vector = c()
  for (k in List_marks) {
    temp_object = ppp_object[ppp_object$marks==k]
    temp_object$marks = NULL
    if (npoints(temp_object)>N_min_points) {
      Kest_temp = Kest(temp_object,r = r_vector,rmax = max(r_vector),correction = "isotropic")
      pcf_temp = pcf.fv(Kest_temp,method = "b") #Pcf is computed by taking the derivative of the K function and by setting pcf(0)=0 as two cells cannot be at the same place 
      Table_pcf = cbind(Table_pcf,pcf_temp$pcf)
      Colnames_vector = c(Colnames_vector,k)
    }
  }
  if (!is.null(ncol(Table_pcf))) {
    colnames(Table_pcf) = Colnames_vector
  }
    
  Table_fits = c()
  for (k in 1:ncol(Table_pcf)) {
    temp_fit = try(Fit_sigmoid(x = r_vector,y = Table_pcf[,k],show_plot = F,initial_tau = initial_tau))
    if (class(temp_fit)=="try-error") {
      temp_fit = rep(NA,5)
    }
    Table_fits = rbind(Table_fits,c(temp_fit))
    
  }
  rownames(Table_fits) = colnames(Table_pcf)
  return(Table_fits)
}

Get_pcf_sigmoid_fit_hyperframe = function(hyperframe_object,N_min_points=100,r_vector = seq(0,0.5,length.out = 100),initial_tau=0.1) {
  N_frames = length(hyperframe_object$pattern)
  List_table_pcf = c()
  for (k in 1:N_frames) {
    List_table_pcf[[k]] = Get_pcf_sigmoid_fit(hyperframe_object$pattern[[k]],N_min_points = N_min_points,r_vector=r_vector,initial_tau = initial_tau)
  }
  return(List_table_pcf)
}

#II) Testing it on different dataset 

#A) Lansing dataset (Lansing Woods Point Pattern)

Lansing_dataset = Get_pcf_sigmoid_fit(lansing)
Lansing_dataset = as.data.frame(Lansing_dataset)
optimal_palette =brewer.pal(length(unique(lansing$marks)), "Spectral")

par(las=1,bty='l',mfrow=c(2,3))
for (k in 1:length(as.character(unique(lansing$marks)))) {
  plot(lansing$x[lansing$marks==as.character(unique(lansing$marks))[k]]*924,
       lansing$y[lansing$marks==as.character(unique(lansing$marks))[k]]*924,
       pch=21,xlim=c(0,924),ylim=c(0,924),xaxs='i',yaxs='i',bg=optimal_palette[k],
       xlab="X location (feet)",ylab = "Y location (feet)",cex.lab=1.3,main=as.character(unique(lansing$marks))[k])
}

par(las=1,bty="l",mar=c(5,8,3,3))
barplot(Lansing_dataset$tau*924,names.arg = rownames(Lansing_dataset),xlab="Tau value (feet)",
        ylab="Tree species",horiz = T,cex.lab=1.5,xlim=c(0,120),col=optimal_palette) 
barplot(Lansing_dataset$p,names.arg = rownames(Lansing_dataset),xlab="p value",
        ylab="Tree species",horiz = T,cex.lab=1.5,xlim=c(0,5),col=optimal_palette) 
barplot(Lansing_dataset$C_normalised,names.arg = rownames(Lansing_dataset),xlab="Cnormalised value",
        ylab="Tree species",horiz = T,cex.lab=1.5,xlim=c(0,0.15),col=optimal_palette) 

#B) Flu dataset ("The influenza virus M2 protein cytoplasmic tail interacts with the M1 protein and influences virus assembly at the site of virus budding")

Flu_analysis = Get_pcf_sigmoid_fit_hyperframe(flu,N_min_points = 30,r_vector = seq(0,1000,length.out = 100),initial_tau = 100)
Merged_flu_analysis = c()
for (k in 1:length(Flu_analysis)) {
  Merged_flu_analysis = rbind(Merged_flu_analysis,cbind(Flu_analysis[[k]],k))
}
colnames(Merged_flu_analysis)[6] = "Sample"
Merged_flu_analysis = as.data.frame(Merged_flu_analysis)
Merged_flu_analysis$Protein = unlist(lapply(strsplit(rownames(Merged_flu_analysis),split = ".",fixed = TRUE),FUN = function(x) {x[1]}))
Merged_flu_analysis$virustype = flu$virustype[Merged_flu_analysis$Sample]
Merged_flu_analysis

U = Merged_flu_analysis[order(Merged_flu_analysis$Protein),]

par(las=1,bty='l')
plot(NULL,xlim=c(0.5,3.5),ylim=c(0,200),xaxt='n',yaxs='i',xlab="Proteins",ylab="Tau parameter value (nm)",cex.lab=1.3)
prism.plots(tau~Protein,U,add=TRUE,pch=21,bg=string.to.colors(U$Protein,brewer.pal(3, "Spectral")),col="black",cex=2)
plotSigBars(tau~Protein,U)

kruskal.test(tau~Protein,U[U$Protein!="M1",])

U_filtered = U[U$Protein=="M2",]
U_filtered = U_filtered[order(U_filtered$virustype),]

par(las=1,bty='l')
plot(NULL,xlim=c(0.5,2.5),ylim=c(1,7),xaxt='n',yaxs='i',xlab="Genotype",ylab="p parameter value",cex.lab=1.3)
prism.plots(p~virustype,U_filtered,add=TRUE,pch=21,bg=string.to.colors(U_filtered$virustype,brewer.pal(3, "Spectral")[1:2]),col="black",cex=2)
plot(NULL,xlim=c(0.5,2.5),ylim=c(0,100),xaxt='n',yaxs='i',xlab="Genotype",ylab="tau parameter value",cex.lab=1.3)
prism.plots(tau~virustype,U_filtered,add=TRUE,pch=21,bg=string.to.colors(U_filtered$virustype,brewer.pal(3, "Spectral")[1:2]),col="black",cex=2)


### Example image of the flu dataset

Plot_flu_M2_sample = function(sample_selected = 1) {
  
}
Example_flu_dataset =flu$pattern[[2]]
par(las=1,bty='l')
plot(Example_flu_dataset$x,Example_flu_dataset$y,xlim=range(Example_flu_dataset$x)*1.05,ylim= range(Example_flu_dataset$y)*1.05,
     xaxs='i',yaxs='i',xlab="X location (nm)",ylab="X location (nm)",cex.lab=1.3,pch=21,
     bg=string.to.colors(Example_flu_dataset$marks,brewer.pal(3, "Spectral")[2:3]),cex=1.3,main="M1 and M2 (WT)")

Example_flu_dataset =flu$pattern[[11]]
plot(Example_flu_dataset$x,Example_flu_dataset$y,xlim=range(Example_flu_dataset$x)*1.05,ylim= range(Example_flu_dataset$y)*1.05,
     xaxs='i',yaxs='i',xlab="X location (nm)",ylab="X location (nm)",cex.lab=1.3,pch=21,
     bg=string.to.colors(Example_flu_dataset$marks,brewer.pal(3, "Spectral")[c(2,1)]),cex=1.3,main="HA and M1 (WT)")
####
Example_flu_dataset =flu$pattern[[8]]
Example_flu_dataset = Example_flu_dataset[Example_flu_dataset$marks=="M2"]
plot(Example_flu_dataset$x,Example_flu_dataset$y,xlim=range(Example_flu_dataset$x)*1.05,ylim= range(Example_flu_dataset$y)*1.05,
     xaxs='i',yaxs='i',xlab="X location (nm)",ylab="X location (nm)",cex.lab=1.2,pch=21,
     bg=string.to.colors(Example_flu_dataset$marks,brewer.pal(3, "Spectral")[3]),cex=1.5,main=" M2 (WT)")


Example_flu_dataset =flu$pattern[[23]]
Example_flu_dataset = Example_flu_dataset[Example_flu_dataset$marks=="M2"]
plot(Example_flu_dataset$x,Example_flu_dataset$y,xlim=range(Example_flu_dataset$x)*1.05,ylim= range(Example_flu_dataset$y)*1.05,
     xaxs='i',yaxs='i',xlab="X location (nm)",ylab="X location (nm)",cex.lab=1.2,pch=21,
     bg=string.to.colors(Example_flu_dataset$marks,brewer.pal(3, "Spectral")[3]),cex=1.5,main=" M2 (mut)")



  
  
