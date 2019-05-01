#' silhouette
#'
#' A more detailed description.
#'
#' @param x A number.
#' @param y A number.
#'
#' @return The sum of \code{x} and \code{y}.
#'
#' @import magrittr
#'
#' @export
# silhoutte general formula
silhouette <- function(X, cluster) {
    # validate input
    # consistency between X and cluster

    # cast X in DF or matrix

    # distance
    Silhouette.Obj <- silhouette.distance.spearman(X)

    # add cluster metadata
    Silhouette.Obj <- silhouette.add_cluster(Silhouette.Obj,
                                             cluster_df = cluster)

    Silhouette.Obj <- Slhouette_coef_df(Silhouette.Obj)
}



Valid.Silhouette.Obj <- function(Silhouette.Obj){
  # silhouette.distance.spearman
  # dim, slots, ...

  # Add_Cluster_metadata
  return(TRUE)
}


silhouette.compute <- function(Silhouette.Obj){

  stopifnot(!is.null(Silhouette.Obj$cluster_df))

  sdf <- Silhouette.Obj %>% split( Silhouette.Obj$molecule.x )
  liste_mol <- names(sdf)
  coef.df <- as.data.frame(matrix(data = NA, nrow = length(names(sdf)), ncol = 3 )) %>%
    set_names(c("molecule", "silhouette.coef", "cluster"))

  for( mol in 1:length(names(sdf))){
    mol_name_tmp <- liste_mol[mol]
    sdf_mol_tmp <- sdf[[mol_name_tmp]]
    tmp <- compute_molecule_silhouette(sdf_mol_tmp, mol_name_tmp)

    coef.df[mol,1] <- mol_name_tmp
    coef.df[mol,2] <- tmp
    coef.df[mol,3] <- get_cluster(sdf_mol_tmp, mol_name_tmp)
  }
  return(coef.df)
}

# compute silhoute coef for 1 molecule
compute_molecule_silhouette <- function(sdf, molecule){
    # compute silhouette coef for \molecule
    A = average_inside(df = sdf, molecule = molecule)
    B = min_average_outside(df = sdf, molecule = molecule)
    sc <- silhoutte_coef(A = A, B = B)
    #df <- data.frame(molecule, sc, get_cluster(sdf, molecule)) %>% set_names(c("molecule", "silhouette.coef", "cluster"))
    #return(df)
}

average_inside <- function(df, molecule){
    tmp_inside <- df %>% # filter(molecule.x == molecule) %>% # good df given by split
        filter(cluster.x == cluster.y) %>% # same cluster WARNINGS:: molecule in +1 cluster #WARNINGS 2 : only 1
        group_by(cluster.y) %>%
        summarise("average" =mean(dist)) %>%
        top_n(n=-1, average) %>% dplyr::select(average) %>% pull
    if(is_empty(tmp_inside)){
        return(0)
    } else {
        return(tmp_inside)
    }
}

# silhoutte general formula
silhoutte_coef <- function(A,B){
    # A average dist inside cluster; # B: min average dist outside
    (B - A)/(max(A,B))
}

get_cluster <- function(sdf, molecule){
    tmp <- sdf %>% dplyr::select(cluster.x) %>% unique() %>% pull
    return(tmp[1])
}



# min (average outside cluster)
min_average_outside <- function(df, molecule){
    df %>% # filter(molelcule.x == molecule) %>%
        filter(cluster.x != cluster.y) %>% # outside cluster
        group_by(molecule.x) %>% # 1 molecule
        summarise("average" =mean(dist)) %>%
        top_n(n=-1, average) %>% dplyr::select(average) %>% pull
}

#' @import dplyr
#' @import magrittr
silhouette.add_cluster  <- function(Silhouette.Obj, cluster_df) {
    # cluster_df =  molecule | cluster
    # valiate cluster_df

    # first left_join : add molecule.x cluster
    Silhouette.Obj$distance_df <-
        left_join(Silhouette.Obj$distance_df,
                  cluster_df,
                  by = c("molecule.x" = "molecule")) %>%
        set_names(c("molecule.x", "molecule.y", "dist", "cluster.x"))

    # second left_join : add molecule.y cluster
    Silhouette.Obj$distance_df <-
        left_join(Silhouette.Obj$distance_df,
                  cluster_df,
                  by = c("molecule.y" = "molecule")) %>%
        set_names(c(
            "molecule.x",
            "molecule.y",
            "dist",
            "cluster.x",
            "cluster.y")) %>%
        # discard distance between same molecules.
        filter(molecule.x != molecule.y)

    Silhouette.Obj$cluster_df <- cluster_df

    Valid.Silhouette.Obj(Silhouette.Obj)
    return(Silhouette.Obj)
}



#' Spearman distance from a matrix
#'
#' TODO
#'
#' @param X A matrix
#'
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import magrittr
#'
silhouette.distance.spearman <- function(X) {
    if (!is(X, "Silhouette.Obj")) {
        # time X molecule, dataframe with rownames(time)
        res.cor <- cor(x = X, use = 'pairwise.complete.obs',
                       method = 'spearman') %>%
            as.data.frame() %>%
            rownames_to_column(var = "molecule.x")
        # 3 columns : molecule.x, molecule.y, dist
        res.cor <- gather(res.cor, molecule.y, dist,-molecule.x)
        # cumpute adjusted spearman distance
        # -rho + 1  ## Spearman like -> distance  [0:2]
        res.cor <- mutate(res.cor, dist = -dist + 1)

        silhouette.res <- list()
        silhouette.res$distance_df <- as.data.frame(res.cor)
        silhouette.res$names <- colnames(X)
        silhouette.res$data <- X

        attr(silhouette.res, "class") <- "Silhouette.Obj"
    }

    Valid.Silhouette.Obj(silhouette.res)
    return(silhouette.res)
}
