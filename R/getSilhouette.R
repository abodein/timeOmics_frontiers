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

#' @importFrom magrittr %>%
getSilhouette.pca <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    silhouette.res <- silhouette(X = as.data.frame(object$X), cluster = cluster)
    ComputeSilhouetteAverage(silhouette.res)
}

#' @importFrom magrittr %>%
#' @importFrom dplyr select
getSilhouette.spca <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- object$X %>% as.data.frame() %>%
        dplyr::select(cluster$molecule)
    silhouette.res <- silhouette(X = X, cluster = cluster)
    ComputeSilhouetteAverage(silhouette.res)
}

#' @importFrom magrittr %>%
getSilhouette.mixo_pls <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- object$X %>% as.data.frame()
    Y <- object$Y %>% as.data.frame()
    data <- cbind(X,Y)
    silhouette.res <- silhouette(X = data, cluster = cluster)
    ComputeSilhouetteAverage(silhouette.res)
}

#' @importFrom magrittr %>%
#' @importFrom dplyr select
getSilhouette.mixo_spls <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- object$X %>% as.data.frame()
    Y <- object$Y %>% as.data.frame()
    data <- cbind(X,Y) %>%
        dplyr::select(cluster$molecule)

    silhouette.res <- silhouette(X = data, cluster = cluster)
    ComputeSilhouetteAverage(silhouette.res)
}

#' @importFrom magrittr %>%
getSilhouette.block.pls <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- do.call("cbind", object$X) %>% as.data.frame

    silhouette.res <- silhouette(X = X, cluster = cluster)
    ComputeSilhouetteAverage(silhouette.res)
}


#' @importFrom magrittr %>%
#' @importFrom dplyr select
getSilhouette.block.pls <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- do.call("cbind", object$X) %>% as.data.frame %>%
        dplyr::select(cluster$molecule)

    silhouette.res <- silhouette(X = X, cluster = cluster)
    ComputeSilhouetteAverage(silhouette.res)
}

#' @export
wrapper.getSilhouette.ncomp <- function(X, min = NULL, max = NULL, ...){
    # test for min, max
    if(!is.null(min)){
        stopifnot(is.numeric(min))
        if(!(min > 0) || !(max <= X$ncomp)){
            stop("min should be an integer greater or equal than 1 and lower than ncomp")
        }
    } else {
        min <- 1
    }
    if(!is.null(max)){
        stopifnot(is.numeric(min))
        if(!(max <= X$ncomp) || !(max > 0)){
            stop("max should be an integer greater or equal than 1 and lower than ncomp")
        }
    } else {
        max <- X$ncomp
    }
    # if min and max have been set
    min <- base::min(min, max)
    max <- base::max(min, max)
    UseMethod("wrapper.getSilhouette.ncomp")
}

wrapper.getSilhouette.ncomp.pca <- function(X, min, max, ...){
    # iteratively compute silhouette for min to max; re-run pca for every ncomp
    pca.res <- list()
    for(comp in min:max){
        data <- X$X  # already scaled/centered
        pca.res[[i]] <- mixOmics::pca(X = data, scale = F, ...)
    }
}

wrapper.getSilhouette.ncomp.mixo_pls <- function(X, min, max){

}
