#' getCluster
#'
#' A more detailed description.
#'
#' @param X A numeric matrix.
#'
#' @param cluster A data.frame
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
#' data <- get_demo_silhouette()
#' res <- silhouette(X = data$data, cluster = data$cluster)
#'
#' @export
# silhoutte general formula
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr %>%
#'
#'
getCluster <- function(x) UseMethod("getCluster")

getCluster.pca <- function(){

    X.block <- do.call("rbind",mixOmics.res$loadings)

    clust.X <- X.block %>% as.data.frame() %>%
        rownames_to_column("molecule") %>%
        gather(comp, value, -molecule) %>%
        mutate(comp = comp %>% str_replace("comp ", "PC")) %>%
        filter_if_sparse(sparse) %>%
        group_by(molecule) %>% mutate(val_abs = abs(value))

    clust.X.abs <- clust.X %>% dplyr::summarise(val_abs = max(val_abs))

    clust.X.2  <- clust.X  %>%
        inner_join(clust.X.abs, by = c("molecule" = "molecule", "val_abs" = "val_abs")) %>%
        filter_cutoff_abs_value(cutoff) %>%
        dplyr::select(-val_abs) %>%
        mutate(comp = comp %>% str_remove("PC") %>% as.numeric) %>%
        mutate(cluster = sign(value)*comp) %>% dplyr::select(c(molecule, cluster))

    return(clust.X.2)
}





getLoadings <- function(X){
    X$loadings %>% reduce("rbind") %>% as.data.frame() %>%
        rownames_to_column("molecule") %>%
        gather(comp, value, -molecule) %>%
        mutate(comp = comp %>% str_replace("comp ", "PC")) %>%
        group_by(molecule) %>% mutate(val_abs = abs(value)) %>%
        mutate(sign = sign(value))
}




loadings.get_cluster <- function(mixOmics.res, sparse = F, cutoff = 0){

    X.block <- do.call("rbind",mixOmics.res$loadings)

    clust.X <- X.block %>% as.data.frame() %>%
        rownames_to_column("molecule") %>%
        gather(comp, value, -molecule) %>%
        mutate(comp = comp %>% str_replace("comp ", "PC")) %>%
        filter_if_sparse(sparse) %>%
        group_by(molecule) %>% mutate(val_abs = abs(value))

    clust.X.abs <- clust.X %>% dplyr::summarise(val_abs = max(val_abs))

    clust.X.2  <- clust.X  %>%
        inner_join(clust.X.abs, by = c("molecule" = "molecule", "val_abs" = "val_abs")) %>%
        filter_cutoff_abs_value(cutoff) %>%
        dplyr::select(-val_abs) %>%
        mutate(comp = comp %>% str_remove("PC") %>% as.numeric) %>%
        mutate(cluster = sign(value)*comp) %>% dplyr::select(c(molecule, cluster))

    return(clust.X.2)
}


filter_if_sparse <- function(df, sparse = F){
    if(sparse){
        df %>% filter(value != 0)
    } else {
        df
    }
}

filter_cutoff_abs_value <- function(df, cutoff = 0){
    if(cutoff){
        df %>% filter(val_abs >= cutoff)
    }else{
        df
    }
}
