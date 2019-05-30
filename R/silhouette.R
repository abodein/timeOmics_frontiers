#' silhouette
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
#'   \item{distance_df}{data.frame containing the cumputed distance between features.}
#'   \item{names}{list of chararcter containing the names of the features.}
#'   \item{data}{matrix containing original data.}
#'   \item{cluster_df}{data.frame containing the cluster information by features.}
#'   \item{coef.df}{data.frame containing ths silhouette coefficient by cluster.}
#' }
#'
#' @examples
#' data <- get_demo_silhouette()
#' res <- silhouette(X = data$data, cluster = data$cluster)
#'
#' @export
# silhoutte general formula
silhouette <- function(X, cluster) {
    # validate input
    CheckSilhouetteInput(X, cluster)

    # distance
    # do nothing if X is already a Silhouette.Obj object
    Silhouette.Obj <- silhouette.distance.spearman(X)

    # add cluster metadata
    Silhouette.Obj <- silhouette.add_cluster(Silhouette.Obj,
                                             cluster_df = cluster)

    Silhouette.Obj <- silhouette.compute_silhouette(Silhouette.Obj)
    return(Silhouette.Obj)
}


#' Get data for silhouette demo
#'
#' @return A matrix of expression profile, sample in raws, time in columns.
#'
#' @examples
#' data <- get_demo_silhouette()
#'
#' @export
get_demo_silhouette <- function() {
    readRDS(system.file("extdata/data_silhouette.rds", package="timeOmics",
                        mustWork = TRUE))
}



Valid.Silhouette.Obj <- function(Silhouette.Obj) {
    # silhouette.distance.spearman
    # dim, slots, ...
    stopifnot(is(Silhouette.Obj, "Silhouette.Obj"))
    stopifnot(length(Silhouette.Obj) %in% c(3,4,5))

    stopifnot(ncol(Silhouette.Obj$distance_df) %in% c(3,5))

    #check data
    stopifnot(ncol(Silhouette.Obj$data) == length(Silhouette.Obj$names))
    stopifnot(colnames(Silhouette.Obj$data) == Silhouette.Obj$names)

    # check distance_df
    stopifnot(nrow(Silhouette.Obj$distance_df) ==
                  length(Silhouette.Obj$names)^2 - length(Silhouette.Obj$names))

    if(length(Silhouette.Obj) == 3){
        # first building step
        stopifnot(ncol(Silhouette.Obj$distance_df) == 3)
    } else if(length(Silhouette.Obj) == 4 || length(Silhouette.Obj) == 5){
        # second building step :: add metadata cluster
        stopifnot(ncol(Silhouette.Obj$distance_df) == 5)
        stopifnot(!is.null(Silhouette.Obj$cluster_df))
        stopifnot(ncol(Silhouette.Obj$cluster_df) == 2)
        stopifnot(colnames(Silhouette.Obj$cluster_df) ==
                      c("molecule", "cluster"))
        stopifnot(Silhouette.Obj$cluster_df$molecule %in% Silhouette.Obj$names)
    } else if(length(Silhouette.Obj) == 5){
        stopifnot(!is.null(Silhouette.Obj$coef.df))
        stopifnot(ncol(Silhouette.Obj$coef.df) == 3)
        stopifnot(nrow(Silhouette.Obj$coef.df) == length(Silhouette.Obj$names))
        stopifnot(names(Silhouette.Obj$coef.df) ==
                      c("molecule", "silhouette.coef","cluster"))
    }
}

CheckSilhouetteInput <- function(X, cluster){
    # X
    stopifnot(is(X, "matrix") || is(X, "data.frame"))
    stopifnot(!any(is.na(X)))
    stopifnot(is.numeric(as.matrix(X)))


    # cluster
    stopifnot(is(cluster, "data.frame"))
    # consistency of dim
    stopifnot(ncol(cluster) == 2)
    stopifnot(!any(is.na(cluster)))
    stopifnot(colnames(cluster) == c("molecule", "cluster"))
    stopifnot(nrow(cluster) == nrow(cluster))
    # no duplicated value
    stopifnot(length(unique(cluster$molecule)) == nrow(cluster))
    # consistency of molecule names
    stopifnot(all(cluster$molecule %in% colnames(data$data)))
}



#' Spearman distance from a matrix
#'
#' Compute spearman correlation distance
#'
#' @param X A matrix
#'
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr %>%
#'
silhouette.distance.spearman <- function(X) {
    if (!is(X, "Silhouette.Obj")) {
        # time X molecule, dataframe with rownames(time)
        res.cor <- cor(x = X, use = 'pairwise.complete.obs',
                       method = 'spearman') %>%
            as.data.frame() %>%
            rownames_to_column(var = "molecule.x")
        # 3 columns : molecule.x, molecule.y, dist
        res.cor <- gather(res.cor, molecule.y, dist,-molecule.x) %>%
            # discard distance between same molecules.
            filter(molecule.x != molecule.y)
        # cumpute adjusted spearman distance
        # -rho + 1  ## Spearman like -> distance  [0:2]
        res.cor <- mutate(res.cor, dist = -dist + 1)

        # initialize Silhouette.Obj
        silhouette.res <- list()
        silhouette.res$distance_df <- as.data.frame(res.cor)
        silhouette.res$names <- colnames(X)
        silhouette.res$data <- X

        attr(silhouette.res, "class") <- "Silhouette.Obj"
    }

    Valid.Silhouette.Obj(silhouette.res)
    return(silhouette.res)
}

#' @import dplyr
#' @importFrom magrittr %>%
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
            "cluster.y"))

    Silhouette.Obj$cluster_df <- cluster_df

    Valid.Silhouette.Obj(Silhouette.Obj)
    return(Silhouette.Obj)
}


#' @import purrr
silhouette.compute_silhouette <- function(Silhouette.Obj) {
    # test if Silhouette.Obj recieve cluster data
    stopifnot(!is.null(Silhouette.Obj$cluster_df))

    sdf <- Silhouette.Obj$distance_df %>%
        split(Silhouette.Obj$distance_df$molecule.x)
    liste_mol <- names(sdf)
    coef.df <-
        as.data.frame(matrix(
            data = NA,
            nrow = length(names(sdf)),
            ncol = 3
        )) %>%
        set_names(c("molecule", "silhouette.coef", "cluster"))

    for (mol in 1:length(names(sdf))) {
        mol_name_tmp <- liste_mol[mol]
        sdf_mol_tmp <- sdf[[mol_name_tmp]]
        tmp <- compute_molecule_silhouette(sdf_mol_tmp, mol_name_tmp)

        coef.df[mol, 1] <- mol_name_tmp
        coef.df[mol, 2] <- tmp
        coef.df[mol, 3] <- get_cluster(sdf_mol_tmp, mol_name_tmp)
    }
    Silhouette.Obj$coef.df <- coef.df
    return(Silhouette.Obj)
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

#' @import dplyr
#' @importFrom magrittr %>%
min_average_outside <- function(df, molecule){
    df %>% # filter(molelcule.x == molecule) %>%
        filter(cluster.x != cluster.y) %>% # outside cluster
        group_by(molecule.x) %>% # 1 molecule
        summarise("average" =mean(dist)) %>%
        top_n(n=-1, average) %>% dplyr::select(average) %>% pull
}
