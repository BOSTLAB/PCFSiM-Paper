#' Frontal Cortex CosMX Dataset
#'
#' This dataset contains spatial and cellular meta-information derived from a 
#' Nanostring CosMX experiment on a Human Prefrontal Cortex FFPE section.
#'
#' @format A data frame with 188,515 rows (individual cells) and 4 variables:
#' \describe{
#'   \item{Cell_ID}{Unique identifier for each cell (usually a character or factor).}
#'   \item{cell_centroid_x}{The x-coordinate of the cell's centroid (numeric).}
#'   \item{cell_centroid_y}{The y-coordinate of the cell's centroid (numeric).}
#'   \item{Clustering}{A factor or character representing the derived cell type or cluster assignment.}
#' }
#' @source \url{https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-frontal-cortex-ffpe-dataset/}
#' @docType data
#' @keywords datasets
#' @name Frontal_cortex_data 
#' @usage data(Frontal_cortex_data)
#' @examples
#' # Load the dataset
#' data(Frontal_cortex_data)
#'
#' # Display the first few rows
#' head(Frontal_cortex_data)
#'
#' # Check the structure
#' str(Frontal_cortex_data)
NULL