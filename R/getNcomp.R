#' getNcomp
#'
#' @export
getNcomp <- function(object, min.ncomp = NULL, max.ncomp = NULL, ...){
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- object
    stopifnot( is(object, c("pca", "spca", "mixo_pls", "mixo_spls", "block.pls", "block.spls")))

    #-- min.ncomp and max.ncomp
    if(!is.null(min.ncomp)){
        stopifnot(is.numeric(min.ncomp))
        if(!(min.ncomp > 0) || !(max.ncomp <= object$ncomp)){
            stop("min.ncomp should be an integer greater or equal than 1 and lower than ncomp")
        }
    } else {
        min.ncomp <- 1
    }
    if(!is.null(max.ncomp)){
        stopifnot(is.numeric(min.ncomp))
        if(!(max.ncomp <= object$ncomp) || !(max.ncomp > 0)){
            stop("max should be an integer greater or equal than 1 and lower than ncomp")
        }
    } else {
        max.ncomp <- object$ncomp
    }
    # if min.ncomp and max have been set
    min.ncomp <- base::min(min.ncomp, max.ncomp)
    max.ncomp <- base::max(min.ncomp, max.ncomp)

    #-- Iterating ncomp
    mixo.call <- object$call
    silhouette.res <- vector(length = max.ncomp - min.ncomp + 1)

    #-- compute dmatrix using spearman dissimilarity
    dmatrix <- dmatrix.spearman.dissimilarity(object$X)

    #-- iterative
    i <- 1
    for(comp in min.ncomp:max.ncomp){
        mixo.call$ncomp <- comp
        mixo.res <- eval(mixo.call)
        cluster <- getCluster(mixo.res)
        # order of feature is the same as colnames(X)
        stopifnot(cluster$molecule == colnames(object$X))
        silhouette.res[i] <- silhouette(dmatrix, cluster$cluster)$average
        i <- i + 1
    }

    to_return <- list()
    to_return[["dmatrix"]] <- dmatrix
    to_return[["silhouette"]] <- silhouette.res
    to_return[["ncomp"]] <- min.ncomp:max.ncomp
    class(to_return) <- "ncomp.tune.silhouette"
    return(invisible(to_return))
}

#' @export
plot.ncomp.tune.silhouette <- function(object){
    stopifnot(is(object, "ncomp.tune.silhouette"))
    plot(x = object$ncomp, y = object$silhouette, type = "b",
         xlab = "Principal Components", ylab = "Average Silhouette Coefficient")
}
