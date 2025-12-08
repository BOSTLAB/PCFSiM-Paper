#' Compute the pair correlation function using Fast Fourier Transform 
#' @param sce A \code{SingleCellExperiment} object.
#' @param N_min_points Minimum number of points required to compute the PCF for a ROI × cluster combination.
#' @param r_vector Numeric vector of radii at which to compute the PCF.
#' @param computation_method Either \code{"direct"}  or \code{"derivative"} The direct method calculates the pcf 
#' directly from the cellular point pattern model (pcf.ppp) using a kernel estimator,
#' which is conceptually simpler but generally slower for very dense point patterns. The derivative method uses the Fast Fourier Transform (FFT) 
#' to efficiently calculate Ripley's K function (Kest.fft), and then obtains the pcf by differentiating the K function for optimal performance on large datasets.
#' @param verbose Logical; whether to print progress information.
#' @return A list containing:
#'   \item{List_pcf}{A list of PCF values for each ROI × cluster combination}
#'   \item{List_r}{A list of matching radius values}
#'   \item{Annotation}{A data.frame linking PCF entries to ROI and cluster}
#' @importFrom spatstat.explore pcf.ppp
#' @export
#' 

Compute_pcf = function(sce,N_min_points=50,  r_vector = seq(0,500,length.out = 100),computation_method = "direct",verbose=FALSE) {
 
  List_clusters = unique(colLabels(sce))
  N_ROI = length(unique(sce$ImageNumber))
  List_pcf_functions = list()
  PP_cell_type_annotation = c()
  PP_FoV_annotation = c()
  
  if (!computation_method%in%c("direct","derivative")){
    stop("Please provide a correct pcf computation method: has to be direct of derivative")
  }
  
  i = 1 
  for (j in List_clusters) {
    if (verbose) {
      print(j)
    }

    for (k in unique(sce$ImageNumber)) {
      sce_temp = sce[,sce$ImageNumber==k & colLabels(sce) %in% j]
      #print(k)
      if (ncol(sce_temp)>=N_min_points) {
        
        ppp_temp = ppp(sce_temp$Location_Center_X,sce_temp$Location_Center_Y,window = owin(xrange = range(sce_temp$Location_Center_X),yrange = range(sce_temp$Location_Center_Y)))
        
        if (computation_method=='derivative') {
          Kest_temp = Kest.fft(ppp_temp,r = r_vector,rmax = max(r_vector),sigma = 1)

          if (sum(!is.na(Kest_temp$border))<=4){
            Kest_temp$border[is.na(Kest_temp$border)] = 1
          }
          
          pcf_temp = pcf.fv(Kest_temp,method = "b") #Pcf is computed by taking the derivative of the K function and by setting pcf(0)=0 as two cells cannot be at the same place 
        }
        
        if (computation_method=='direct') {
          pcf_temp = pcf.ppp(ppp_temp,r = r_vector,rmax = max(r_vector),kernel = "rectangular",correction = "translate",divisor = "r")
        }

        List_pcf_functions[[i]] = pcf_temp
        
        PP_cell_type_annotation = c(PP_cell_type_annotation,j)
        PP_FoV_annotation = c(PP_FoV_annotation,k)
        i = i + 1 
        
      }
      
    }
    
  }
  
  if (computation_method=='derivative') {
    List_pcf_values = lapply(List_pcf_functions,FUN = function(x) {x$pcf})
  }
  if (computation_method=='direct') {
    List_pcf_values = lapply(List_pcf_functions,FUN = function(x) {x$trans})
  }
  
  List_r_values = lapply(List_pcf_functions,FUN = function(x) {x$r})
  return(list(List_pcf =List_pcf_values,List_r = List_r_values,Annotation = data.frame(ROI =PP_FoV_annotation,Cluster =PP_cell_type_annotation )))
}

