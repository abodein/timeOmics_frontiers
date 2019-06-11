#' plotLong
#'
#' A more detailed description.
#'
#' @param X A numeric matrix.
#'
#' @param (s)pca, (s)pls, block.(s)pls results
#' that contains clustering information with molecule and cluster
#'
#'
#' @return
#' \describe{
#'   \item{One}{First item}
#'   \item{Two}{Second item}
#' }
#'
#' @examples
#'
#'
#' @export
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import stringr
#' @import ggplot2
#' @importFrom magrittr %>%
plotLong <- function(object, time = NULL, plot = TRUE, center = TRUE, scale = TRUE){
    #allowed_object = c("pca", "spca", "mixo_pls", "mixo_spls", "block.pls", "block.spls")
    #stopifnot(is(object, allowed_object))
    UseMethod("plotLong")
}

plotLong.pca <- function(object, time = NULL, plot = TRUE, center = TRUE, scale = TRUE){
    print("plotLong.pca")

    if(!is.null(time)){
        stopifnot(is.numeric(time)) # works with integer
        stopifnot(length(time) == nrow(object$X))

    } else{ # time IS NULL
        # we rely on rownames and assume it's numeric,  correspond to times
        time <- as.numeric(rownames(object$X))
        # not numeric value can introduce NA
        stopifnot(!is.na(time))
    }
    # unscale and rescale if desired
    data <- rebuild(object) %>% scale(scale, center= TRUE)

    # cluster info
    cluster <- getCluster(object)

    plotLongGGplot(data = data, time = time, cluster = cluster)
}


plotLongGGplot <- function(data, time, cluster, plot = TRUE){
    # graphical call
    stopifnot(is.numeric(time)) # works with integer
    stopifnot(nrow(data) == length(time))

    data.gather <- data %>% as.data.frame() %>%
        mutate(time = time) %>%
        gather(molecule, value, -time) %>%
        left_join(cluster, by = c("molecule"="molecule"))

    gg <- ggplot(data.gather, aes(x = time, y = value, group = molecule)) +
        geom_line(aes(color = cluster)) +
        facet_grid(contribution ~ comp, scales = "free") +


    if(plot){
        print(gg)
    }
    return(invisible(data.gather))
}

