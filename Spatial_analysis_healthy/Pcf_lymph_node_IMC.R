#I)Loading function/object
library(RColorBrewer)
#A)List functions

source("Desktop/Projet_analyse_pcf/List_function_pcf_analysis.R")

#B)Sce object
sce= readRDS("Panorama_LN_sce_object.rds") #Can be downloaded from the zenodo associated to the publication Bost et al. 2023 Nature Methods


#II)Computation itself 

#A)Computing pcf
List_pcf = Compute_pcf(sce,r_vector = seq(0,1000,length.out=50),computation_method = "derivative")


#B)Checking the functions on an example
k = 10
x = List_pcf$List_r[[k]]
y = List_pcf$List_pcf[[k]]

m = try(Fit_generalized_gamma(x,y,show_plot = TRUE),silent = TRUE)
m = try(Fit_gamma(x,y,show_plot = TRUE))
Fit_power_law(x,y,show_plot = TRUE)


#C)Checking the functions on an example

Results_generalised_gamma = Fit_parametric_pcf_model(List_pcf,model = "Generalized_gamma")
Results_gamma= Fit_parametric_pcf_model(List_pcf,model = "Gamma")
Results_exponential = Fit_parametric_pcf_model(List_pcf,model = "Exponential")
Results_power_law = Fit_parametric_pcf_model(List_pcf,model = "Power_law")
Results_sigmoid= Fit_parametric_pcf_model(List_pcf,model = "Sigmoid")
Results_beta_prime= Fit_parametric_pcf_model(List_pcf,model = "Beta_prime")

pdf("Desktop/Papier_Bost_2/Plot_Figure_pcf/Pcf_comparison_Lymph_node_IMC.pdf",width = 7,height = 6.5)
par(las=1,bty='l')
boxplot(Results_power_law$R2,Results_exponential$R2,Results_gamma$R2,
        Results_generalised_gamma$R2,Results_sigmoid$R2,Results_beta_prime$R2,ylim=c(0.4,1),
        xaxs='i',yaxs='i',names = c("Power law","Exponential","Gamma","GG","Sigmoid","Beta prime"),
        xlab="Model",ylab="R2",cex.lab=1.3,col=brewer.pal(6, "Spectral"),cex.axis=0.7)
dev.off()


#III) Comparing the Pcf approach with more regular approaches

#A)CE index

CE_index = CE_interaction_tensor(sce)
CE_index_vector = 2^diag((CE_index[,,1])@data)

pdf("Desktop/Papier_Bost_2/Figure_S2//Correlation_CE_IMC_lymph_node.pdf",width = 6.5,height = 6.5,useDingbats = F)
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

dev.off()


#B) Study of the over-dispersion index

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

pdf("Desktop/Papier_Bost_2/Figure_S2//Correlation_NB_fit_IMC_lymph_node.pdf",width = 6.5,height = 6.5,useDingbats = F)
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

dev.off()





