library(irlba)
library(igraph)
library(N2R)
library(randomForest)

List_table = list.files("Desktop/CRC_pcf_analysis/Data/",full.names = TRUE)
Names_samples = list.files("Desktop/CRC_pcf_analysis/Data/",full.names = FALSE)
Names_samples = substr(Names_samples,start = 1,stop = 5)

Merged_IF_data = c()
Merged_metadata = c()

for (k in 1:length(List_table)) {
  x = read.delim(List_table[k],sep = ",")
  Merged_IF_data = rbind(Merged_IF_data,x[,11:36])
  Merged_metadata = rbind(Merged_metadata,cbind(x[,c(37,38,44,45)],rep(Names_samples[k],nrow(x))))
}


hist(log10(Merged_metadata$AREA),100)
Subsampling = sample(1:nrow(Merged_metadata),size = nrow(Merged_metadata)*0.05,replace = F )

Sampled_data = Merged_IF_data[Subsampling,]
Mean_IF = colMeans(Sampled_data)
Var_IF = apply(Sampled_data,MARGIN = 2,FUN = var)

par(las=1,bty="l")
plot(Mean_IF,Var_IF,pch=21,bg="red3",log="xy")
text(Mean_IF,Var_IF,labels = names(Mean_IF))

SVD_IF = irlba(as.matrix(Sampled_data),scale = sqrt(Var_IF),nv = 20)

KNN_matrix = Knn(SVD_IF$u, k = 30, verbose = TRUE, 
                 indexType = "L2", nThreads = 10)
KNN_matrix = t(KNN_matrix) + KNN_matrix
KNN_graph = graph_from_adjacency_matrix(KNN_matrix, weighted = TRUE, 
                                        mode = "undirected")
Clustering = cluster_louvain(KNN_graph)
Clustering = membership(Clustering)
Clustering = as.character(Clustering)

boxplot(Sampled_data$CD20~Clustering,outline=F)
Mean_clustering = aggregate(Sampled_data,FUN = mean,by=list(Clustering))
library(pheatmap)
rownames(Mean_clustering) = Mean_clustering$Group.1
pheatmap(t(Mean_clustering[,-1]),scale="row",clustering_method = "ward")

colnames(Merged_metadata)[5]="Sample"
K = "CRC06"
plot(Merged_metadata$Xt[Merged_metadata$Sample==K],
     Merged_metadata$Yt[Merged_metadata$Sample==K],cex=0.1,
     col=string.to.colors(Final_pred[Merged_metadata$Sample==K]),pch=16)

pred_rf = randomForest::randomForest(as.matrix(Sampled_data),y = factor(Clustering))
Final_pred = predict(pred_rf,as.matrix(Merged_IF_data))
Merged_metadata$Cluster = Final_pred

write.table(Merged_metadata,"Desktop/CRC_pcf_analysis/Annotated_data.txt",sep="\t",quote=F)
write.table(Mean_clustering,"Desktop/CRC_pcf_analysis/Mean_profile_cluster.txt",sep="\t",quote=F)
