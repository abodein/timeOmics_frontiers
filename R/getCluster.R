#' getCluster
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
#' demo <- suppressMessages(get_demo_cluster())
#' pca.cluster <- getCluster(demo$pca)
#' spca.cluster <- getCluster(demo$spca)
#' pls.cluster <- getCluster(demo$pls)
#' spls.cluster <- getCluster(demo$spls)
#' block.pls.cluster <- getCluster(demo$block.pls)
#' block.spls.cluster <- getCluster(demo$sblock.pls)
#'
#' @export
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import stringr
#' @importFrom magrittr %>%
getCluster <- function(x) UseMethod("getCluster")

get_demo_cluster<- function(){
    X <- matrix(sample(1:1000), nrow = 10)
    rownames(X) <- paste0("l_",1:nrow(X))
    colnames(X) <- paste0("X_",1:ncol(X))

    Y <- matrix(sample(1:100), nrow = 10)
    rownames(Y) <- paste0("l_",1:nrow(Y))
    colnames(Y) <- paste0("Y_",1:ncol(Y))

    Z <- matrix(sample(1:500), nrow = 10)
    rownames(Z) <- paste0("l_",1:nrow(Z))
    colnames(Z) <- paste0("Z_",1:ncol(Z))

    list.res = list()
    list.res$X <- X
    list.res$Y <- Y
    list.res$Z <- Z
    list.res$pca = mixOmics::pca(X = X, ncomp = 5)
    list.res$spca = mixOmics::spca(X = X, ncomp = 5, keepX = c(5, 15, 4,5,6))

    list.res$pls = mixOmics::pls(X = X, Y = Y, ncomp = 5, mode = "canonical")
    list.res$spls = mixOmics::spls(X = X, Y = Y, ncomp = 5, mode = "canonical",
                                   keepX = c(5,6,4,5,6), keepY = c(5,1,4,5,6))

    list.res$block.pls = mixOmics::block.pls(X = list("X" = X, "Y" = Y, "Z" = Z), indY = 1,
                                             ncomp = 5, mode = "canonical")
    list.res$block.spls = mixOmics::block.spls(X = list("X" = X, "Y" = Y, "Z" = Z), indY = 1, ncomp = 3,
                                             mode = "canonical", keepX = list("X" = c(5,6,4), "Y" = c(5,5,5), "Z" = c(4,2,4)))
    return(list.res)
}

#' @import dplyr
#' @import tibble
#' @import stringr
#' @importFrom magrittr %>%
getCluster.pca <- function(X){
    print(class(X))
    # colnames = PC1, PC2...
    loadings.max <- getMaxContrib(X$loadings$X)

    loadings.max %>% rownames_to_column("molecule") %>%
        mutate(cluster = str_remove(comp, "^PC") %>% as.numeric()) %>%
        mutate(cluster = cluster * sign(contrib.max)) %>%
        mutate(block = "X")
}

#' @import dplyr
#' @import tibble
#' @import stringr
#' @importFrom magrittr %>%
getCluster.spca <- function(X){
    print(class(X))
    selected.features.loadings <- X$loadings$X[rowSums(X$loadings$X) != 0,]
    loadings.max <- getMaxContrib(selected.features.loadings)

    loadings.max %>% rownames_to_column("molecule") %>%
        mutate(cluster = str_remove(comp, "^PC") %>% as.numeric()) %>%
        mutate(cluster = cluster * sign(contrib.max)) %>%
        mutate(block = "X")
}

#' @import dplyr
#' @import tibble
#' @import stringr
#' @importFrom magrittr %>%
getCluster.mixo_pls <- function(X){
    print(class(X))
    # block X
    loadings.max.X <- getMaxContrib(X$loadings$X)

    loadings.max.X <- loadings.max.X %>% rownames_to_column("molecule") %>%
        mutate(cluster = str_remove(comp, "^comp ") %>% as.numeric()) %>%
        mutate(cluster = cluster * sign(contrib.max)) %>%
        mutate(block = "X")

    # block Y
    loadings.max.Y <- getMaxContrib(X$loadings$Y)

    loadings.max.Y <- loadings.max.Y %>% rownames_to_column("molecule") %>%
        mutate(cluster = str_remove(comp, "^comp ") %>% as.numeric()) %>%
        mutate(cluster = cluster * sign(contrib.max)) %>%
        mutate(block = "Y")

    rbind(loadings.max.X, loadings.max.Y)
}

#' @import dplyr
#' @import tibble
#' @import stringr
#' @importFrom magrittr %>%
getCluster.mixo_spls <- function(X){
    # note : can not concatenate X and Y
    # because they can have the same features names contrary to block.(s)pls

    print(class(X))
    # block X
    X.selected.features.loadings <- X$loadings$X[rowSums(X$loadings$X) != 0,]
    loadings.max.X <- getMaxContrib(X.selected.features.loadings)

    loadings.max.X <- loadings.max.X %>% rownames_to_column("molecule") %>%
        mutate(cluster = str_remove(comp, "^comp ") %>% as.numeric()) %>%
        mutate(cluster = cluster * sign(contrib.max)) %>%
        mutate(block = "X")

    # block Y
    Y.selected.features.loadings <- X$loadings$Y[rowSums(X$loadings$Y) != 0,]
    loadings.max.Y <- getMaxContrib(Y.selected.features.loadings)

    loadings.max.Y <- loadings.max.Y %>% rownames_to_column("molecule") %>%
        mutate(cluster = str_remove(comp, "^comp ") %>% as.numeric()) %>%
        mutate(cluster = cluster * sign(contrib.max)) %>%
        mutate(block = "Y")

    rbind(loadings.max.X, loadings.max.Y)
}

#' @import purrr
#' @import stringr
#' @import dplyr
#' @import tibble
getCluster.block.pls <- function(X){
    print(class(X))
    # get block info
    block.info <- purrr::imap(X$loadings, function(x,y) rownames(x) %>%
                           as.data.frame %>%
                           set_names("molecule") %>%
                           mutate("block" = y))
    block.info <- do.call("rbind", block.info)

    loadings <- do.call("rbind", X$loadings)
    X.selected.features.loadings <- loadings[rowSums(loadings) != 0,]
    loadings.max <- getMaxContrib(X.selected.features.loadings)

    loadings.max <- loadings.max %>% rownames_to_column("molecule") %>%
        mutate(cluster = str_remove(comp, "^comp ") %>% as.numeric()) %>%
        mutate(cluster = cluster * sign(contrib.max))
    return(loadings.max)
}


#' @import purrr
#' @import stringr
#' @import dplyr
#' @import tibble
getCluster.block.spls <- function(X){

    print(class(X))
    # get block info
    block.info <- purrr::imap(X$loadings, function(x,y) rownames(x) %>%
                                  as.data.frame %>%
                                  set_names("molecule") %>%
                                  mutate("block" = y))
    block.info <- do.call("rbind", block.info)

    loadings <- do.call("rbind", X$loadings)
    loadings.max <- getMaxContrib(loadings)

    loadings.max <- loadings.max %>% rownames_to_column("molecule") %>%
        mutate(cluster = str_remove(comp, "^comp ") %>% as.numeric()) %>%
        mutate(cluster = cluster * sign(contrib.max))
    return(loadings.max)
}


#' @importFrom purrr set_names
#' @importFrom magrittr %>%
getMaxContrib <- function(X){
    # loadings matrix, features in rows, comp in columns
    contrib.max <- apply(X = X, FUN = function(x) { x[which.max( abs(x) )][1]}, MARGIN = 1) %>%
        as.data.frame() %>% purrr::set_names("contrib.max")

    cluster.info <- apply(X = X, FUN = function(x) { colnames(X)[which.max( abs(x) )[1]]}, MARGIN = 1) %>%
        as.data.frame() %>% purrr::set_names("comp")

    stopifnot(rownames(contrib.max) == rownames(cluster.info))
    return(cbind(cluster.info, contrib.max))
}

# absmax <- function(x) { x[which.max( abs(x) )][1]}
# absmax.index <- function(x) { which.max( abs(x) )[1]}





