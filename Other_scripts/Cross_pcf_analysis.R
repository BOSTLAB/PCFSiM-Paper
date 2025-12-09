library(pheatmap)
library(fifer)
library(FactoMineR)

#I)Data and function loading

source("Desktop/Projet_analyse_pcf/List_function_pcf_analysis.R")

Cell_data = read.delim("~/Desktop/Projet_analyse_pcf/Mice_colitis_model/cell_properties.csv",sep=",")


r_vector = seq(0,500,length.out = 100)
N_min_points = 100
List_table_pcf = c()
Clustering = Cell_data$Tier2
List_pair_clustering = expand.grid(unique(Clustering),unique(Clustering),stringsAsFactors = F)
List_pair_clustering = matrix(data = paste(List_pair_clustering$Var1,List_pair_clustering$Var2,sep = "_"),ncol = length(unique(Clustering)))
List_pair_clustering = as.character(List_pair_clustering[upper.tri(List_pair_clustering)])

for (k in unique(Cell_data$Slice_ID)) {
  print(k)
  Table_pcf=matrix(1,nrow = length(r_vector),ncol = length(List_pair_clustering))
  colnames(Table_pcf) = List_pair_clustering
  for (j in List_pair_clustering) {
    print(j)
    x = strsplit(j,split = "_")
    Cell_type_1 = unlist(lapply(x,FUN = function(x){x[1]}))
    Cell_type_2 = unlist(lapply(x,FUN = function(x){x[2]}))
    
    Selected_cells = Cell_data$Slice_ID==k & Clustering%in%c(Cell_type_1,Cell_type_2)
    ppp_temp = ppp(Cell_data$x[Selected_cells],
                   Cell_data$y[Selected_cells],
                   window = owin(xrange = range(Cell_data$x[Cell_data$Slice_ID==k]),
                                 yrange = range(Cell_data$y[Cell_data$Slice_ID==k])),marks = factor(Clustering[Selected_cells]))
    
    if (sum(Clustering[Selected_cells]==Cell_type_1)>N_min_points & sum(Clustering[Selected_cells]==Cell_type_2)>N_min_points ) {
      Kest_temp = Kcross(ppp_temp,r = r_vector,rmax = max(r_vector),correction = "isotropic")
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


PCF_computed = c()
for (k in 1:length(List_table_pcf)) {
  for (j in 1:ncol(List_table_pcf[[k]])) {
    u = List_table_pcf[[k]][,j]
    z = identical(u,rep(1,length(u)))
    if (z) {
      PCF_computed = c(PCF_computed,FALSE)
    }
    else {
      PCF_computed = c(PCF_computed,TRUE)
      
    }
  }
}


colnames(Table_fits) = c("tau","p", "C","C_normalised","R2","Sample","Cluster_pairs")
Table_fits= as.data.frame(Table_fits)

#C)Plotting the results

boxplot(Table_fits$R2,ylim=c(0,1),yaxs='i',ylab="R2")

N_sample_fit = sample((1:ncol(Table_pcf))[!is.na(Table_fits$R2)],size = 5,replace = F)

Table_fits_filtered = Table_fits[PCF_computed,]

#24% of the models could note be fitted
#16 % is low quality -> 40% of the pcf could note be fitted correctly 
Failure_fit = is.na(Table_fits_filtered$R2) | Table_fits_filtered$R2<0.6
par(las=1,bty='l')
barplot(table(Failure_fit)/length(Failure_fit)*100,
        names.arg = c("Correct pcf fit","No correct fit"),col=c('grey',"darkred"),ylim=c(0,100),
        ylab="Proportion of models (%)",cex.lab=1.3)

#D)Looking at individual pcfs

N_samples_per_class = 5 #Taking 5 pcfs of interactions that worked and failed
N_sample_no_fit = sample((1:nrow(Table_fits_filtered))[Failure_fit],size = N_samples_per_class,replace = F)
N_sample_fit = sample((1:nrow(Table_fits_filtered))[!Failure_fit],size = N_samples_per_class,replace = F)

U_fit = c()
U_failed_fit = c()

X = Table_fits_filtered[N_sample_fit,]
for (k in 1:N_samples_per_class) {
  u = List_table_pcf[[X[k,"Sample"]]][,X[k,"Cluster_pairs"]]
  U_fit = rbind(U_fit,u)
}

X = Table_fits_filtered[N_sample_no_fit,]
for (k in 1:N_samples_per_class) {
  u = List_table_pcf[[X[k,"Sample"]]][,X[k,"Cluster_pairs"]]
  U_failed_fit = rbind(U_failed_fit,u)
}

par(las=1,bty='l')
plot(NULL,xlim=c(0,500),ylim=c(0,max(U_failed_fit,na.rm = T)*1.05),xaxs='i',yaxs='i',xlab="r (µm)",ylab="Pcf(r)",cex.lab=1.3,main="Non fitted models")
for (k in 1:nrow(U_failed_fit)) {
  points(r_vector,U_failed_fit[k,],type='l',col="darkred")
  
}
abline(h=1,lwd=2,lty=2,col="black")

par(las=1,bty='l')
plot(NULL,xlim=c(0,500),ylim=c(0,max(U_fit,na.rm = T)*1.05),xaxs='i',yaxs='i',xlab="r (µm)",ylab="Pcf(r)",cex.lab=1.3,main="Fitted models")
for (k in 1:nrow(U_fit)) {
  points(r_vector,U_fit[k,],type='l',col="grey")
  
}
abline(h=1,lwd=2,lty=2,col="black")



