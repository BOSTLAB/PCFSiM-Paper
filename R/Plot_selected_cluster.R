#' Plot cells of a selected cluster
#'
#' @param selected_cluster Cluster name or ID to plot.
#' @param title_show Optional title for the plot.
#' @return A plot of the point pattern for the selected cluster.
#' @importFrom spatstat.geom ppp
#' @importFrom spatstat.geom owin
#' @export

Plot_selected_cluster = function(Meta_data = Meta_data, selected_cluster,title_show=NULL) {
  ppp_temp = ppp(x = Meta_data$cell_centroid_x[Meta_data$Clustering==selected_cluster],
                 y = Meta_data$cell_centroid_y[Meta_data$Clustering==selected_cluster],
                 window = owin(xrange = range(Meta_data$cell_centroid_x),yrange = range(Meta_data$cell_centroid_y)))
  plot(ppp_temp,cex=0.2,main=title_show)
  
}
