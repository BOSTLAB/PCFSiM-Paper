library(TranspaceR)
#-------------------------------------------PATHS AND SAMPLE DEFINITION---------------------------------
path = '~/Spatial_atlas/Datasets/CosMx/Frontal_cortex/'
Output_path = '~/Spatial_atlas/Picture_pipeline/CosMx/Frontal_cortex/'
Method = "CosMx"
Tissue = "Frontal_cortex"
#-------------------------------------------------------------------------------------------------------
# Data loading 
Expression_file = read.csv(paste(path, "CosMx_Frontal_cortex_ExprMat.csv", sep = ""), header = TRUE)
Meta_data = read.csv(paste(path,"CosMx_Frontal_cortex_Metadata.csv",sep =""), header = TRUE) 
colnames(Meta_data)[c(39,40)] <- c("cell_centroid_x", "cell_centroid_y")
# From here, Expression_file should be numeric and count matrix only, cell ID column is deleted
Expression_file = Expression_file[,-c(1,2,3)]

# Quality control
Meta_data = Compute_radius(Meta_data,Method)
QC1_results =  QC_RNA_size_threshold(Expression_file,Meta_data,Method,Tissue,Output_path)
QC2_results =  QC_Gene_threshold(Expression_file,Meta_data,Method,Tissue,Output_path)

# Data filtering with QC results
Data_curated = Curate_data(Expression_file,Meta_data,QC2_resultsmin_lib_size=10,max_lib_size=9500,min_cell_radius=2,max_cell_radius=15)
Expression_file = Data_curated$Expression_file
Meta_data = Data_curated$Meta_data
rm(Data_curated)

# Variable gene Selection
# Based on significative variance 
Variance_computation = Excess_variance_ratio_NB(Expression_file,Output_path,Method,Tissue)
Variance_genes= Variance_computation$Selected_genes
# Based on high zero Proportion 
Zero_score_computation = Excess_zero_score_NB(Expression_file,Output_path,Method,Tissue,P_value_threshold=0.01,Delta_threshold = 0.01)
Zero_genes = Zero_score_computation$Selected_genes
# Based on Spatial variance
# Geary's C
Geary_computation =Geary_C_score(Expression_file,Meta_data,Output_path,Method,Tissue,pvalue=0.01)
Geary_genes = Geary_computation$Selected_genes

# All selected genes from selected score methods
Shared_genes = Select_genes(Selected_objects =list(Variance_genes,Zero_genes,Geary_genes),
                            Selected_names = c('Variance_genes','Zero_genes','Geary_genes'))
# Normalization by cell size and Clustering
Clustering_output = Clusters_maker(Expression_file, Shared_genes, K=50, metric_used="L2", nThreads = 20, resolution = 1)
Clustering = Clustering_output$Clustering  
PCA_data = Clustering_output$PCA_data
Mean_expression = Clustering_output$Mean_expression
Data_correction = Clustering_output$Data_correction
Log2FC_table = Clustering_output$Log2FC_table
Save_heatmap_markers(Expression_file, object = Clustering,Output_path,Method,Tissue)
Save_tissue_visualization(Meta_data,object = Annotation_cells_renamed, Output_path,Method,Tissue,name_object = 'Clustering', scaling_factor=1)
