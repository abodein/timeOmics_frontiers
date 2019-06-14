#' getSilhouette
#'
#' \code{getSilhouette} is a generic function that compute silhouette coefficient
#' for an object of the type \code{pca}, \code{spca}, \code{pls}, \code{spls},
#' \code{block.pls}, \code{block.spls}.
#'
#' @param X A numeric matrix.
#'
#' @param cluster A data.frame
#' that contains clustering information with molecule and cluster
#'
#' @details
#' This method extract the componant contribution depending on the object,
#' perform the clustering step, and compute the silhouette coefficient.
#'
#' @return
#'
#' @examples
#' data <- get_demo_silhouette()
#' res <- silhouette(X = data$data, cluster = data$cluster)
#'
#' @export
getSilhouette <- function(object){
    UseMethod("getSilhouette")
}


getSilhouette.pca <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")
}
