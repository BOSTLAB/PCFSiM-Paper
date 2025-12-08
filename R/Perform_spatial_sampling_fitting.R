#' Perform spatial subsampling and fit parametric PCF models
#'
#' Perform random spatial sampling on a \code{SingleCellExperiment} object, compute the PCF on each sampled region, 
#' and fit a parametric PCF model (Sigmoid). 
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param N_min_cells Minimum number of cells required to compute the PCF.
#' @param FoV_width Width (in pixels or microns) of the square Field of View used for sampling.
#' @param N_sampling Number of spatial samplings to perform before filtering.
#'
#' @return The function returns the distribution of fitted parameters
#' across all valid sampled regions :
#' \describe{
#'   \item{Tau}{Matrix of \eqn{\tau} values across samplings.}
#'   \item{P}{Matrix of \eqn{p} values across samplings.}
#'   \item{C_normalised}{Matrix of normalised amplitudes.}
#'   \item{R2}{Matrix of model \eqn{R^2} values.}
#' }
#' @importFrom SingleCellExperiment colLabels
#' @importFrom spatstat.explore Kest.fft pcf.fv
#' @export

Perform_spatial_sampling_fitting = function(sce,N_min_cells = 1000,FoV_width=500,N_sampling = 500) {
  List_sampled_cells = c()
  cat("Perfoming sampling....")
  for (k in 1:N_sampling) {
    print(k)
    x = balagan::Random_spatial_sampling(sce,width_FOV = FoV_width,height_FOV = FoV_width,N_samplings = 1,plot_result = FALSE)
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
