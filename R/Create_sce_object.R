#' Create a SingleCellExperiment object from a data.frame
#'
#' @param df A data.frame containing at least coordinates and Labels info.
#' @param cell_centroid_x Name of the column for x coordinates (default: "cell_centroid_x").
#' @param cell_centroid_y Name of the column for y coordinates (default: "cell_centroid_y").
#' @param Labels Name of the column for cluster labels (default: "Labels").
#' @return A SingleCellExperiment object with coordinates and labels.
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SingleCellExperiment colLabels
#' @export 

Create_sce_object = function(df, cell_centroid_x=2, cell_centroid_y=3, Labels=4) {
  if (!all(is.numeric(c(cell_centroid_x, cell_centroid_y, Labels)))) {
    stop("Column arguments must be numeric indices (e.g., 2, 3, 4), not column names.")
  }
  if (!all(c(cell_centroid_x, cell_centroid_y, Labels) %in% 1:ncol(df))) {
    stop("One or more specified column indices are out of bounds for the input data.frame.")
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(Raw_intensity = matrix(0, nrow = 2, ncol = nrow(df))),
    metadata = list(dimension = "2D",N_core = 8,Is_nuc_cyt = FALSE))
  
  SingleCellExperiment::colLabels(sce) = as.numeric(df[[Labels]])

  sce$Location_Center_X = df[[cell_centroid_x]]
  sce$Location_Center_Y = df[[cell_centroid_y]]
  sce$ImageNumber = 1
  
  return(sce)
}