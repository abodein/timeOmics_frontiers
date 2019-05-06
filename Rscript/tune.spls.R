
tune.spls2 <- function(X,Y, ncomp, keepX, keepY) {

  mean_sc <- list()
  grid_min_keepX <- keepX %>% lapply(min) %>% unlist()
  grid_min_keepY <- keepY %>% lapply(min) %>% unlist()

  DF_X <- Spearman_distance(X)
  DF_Y <- Spearman_distance(Y)

  for(comp in 1:ncomp){
    # init
    tmp <- keepY[[comp]] %*% t(keepX[[comp]]) %>% as.data.frame()
    colnames(tmp) <- keepX[[comp]]
    rownames(tmp) <- keepY[[comp]]
    mean_sc[[comp]] <- list( "pos" = tmp,
                             "neg" = tmp,
                             "other" = tmp)
    kkX <- grid_min_keepX
    kkY <- grid_min_keepY


    for(kX in keepX[[comp]]){
      kkX[[comp]] <- kX

      for(kY in keepY[[comp]]){
        kkY[[comp]] <- kY

        spls.res <- mixOmics::spls(X=X, Y=Y, ncomp = ncomp, mode = "canonical", keepX = kkX, keepY = kkY)

        # with X
        tmp_cluster <- spls.res$loadings$X[,comp] %>% sign() %>% as.data.frame() %>% rownames_to_column %>%
          set_names(c("molecule", "cluster"))
        # add cluster metadata and compute silhouette
        tmp_DF_w_clu <- Add_Cluster_metadata(DF_X, tmp_cluster)
        tmp_silhouette <- Slhouette_coef_df(tmp_DF_w_clu)
        # mean by cluster
        tmp_cluster_mean_sc <- tmp_silhouette %>% group_by(cluster) %>% summarise(MSC = mean(silhouette.coef))
        # store information
        tmp_cluster_mean_sc %>% filter(cluster == 1) %>% pull(MSC)

        # pos
        if(!is_empty(tmp_cluster_mean_sc %>% filter(cluster == 1) %>% pull(MSC))){
          mean_sc[[comp]]$pos[ which(rownames(mean_sc[[comp]]$pos)== kY), which(colnames(mean_sc[[comp]]$pos)== kX)] <- tmp_cluster_mean_sc %>% filter(cluster == 1) %>% pull(MSC)
        } else{
          mean_sc[[comp]]$pos[ which(rownames(mean_sc[[comp]]$pos)== kY), which(colnames(mean_sc[[comp]]$pos)== kX)] <- 0
        }
        if(!is_empty(tmp_cluster_mean_sc %>% filter(cluster == -1) %>% pull(MSC))){
          mean_sc[[comp]]$neg[ which(rownames(mean_sc[[comp]]$neg)== kY), which(colnames(mean_sc[[comp]]$neg)== kX)] <- tmp_cluster_mean_sc %>% filter(cluster == -1) %>% pull(MSC)
        } else{
          mean_sc[[comp]]$neg[ which(rownames(mean_sc[[comp]]$neg)== kY), which(colnames(mean_sc[[comp]]$neg)== kX)] <- 0
        }
        if(!is_empty(tmp_cluster_mean_sc %>% filter(cluster == 0) %>% pull(MSC))){
          mean_sc[[comp]]$other[ which(rownames(mean_sc[[comp]]$other)== kY), which(colnames(mean_sc[[comp]]$other)== kX)] <- tmp_cluster_mean_sc %>% filter(cluster == 0) %>% pull(MSC)
        } else{
          mean_sc[[comp]]$other[ which(rownames(mean_sc[[comp]]$other)== kY), which(colnames(mean_sc[[comp]]$other)== kX)] <- 0
        }

        # with Y
        tmp_cluster <- spls.res$loadings$Y[,comp] %>% sign() %>% as.data.frame() %>% rownames_to_column %>%
          set_names(c("molecule", "cluster"))
        # add cluster metadata and compute silhouette
        tmp_DF_w_clu <- Add_Cluster_metadata(DF_Y, tmp_cluster)
        tmp_silhouette <- Slhouette_coef_df(tmp_DF_w_clu)
        # mean by cluster
        tmp_cluster_mean_sc <- tmp_silhouette %>% group_by(cluster) %>% summarise(MSC = mean(silhouette.coef))
        # store information
        if(!is_empty(tmp_cluster_mean_sc %>% filter(cluster == 1) %>% pull(MSC))){
          mean_sc[[comp]]$pos[ which(rownames(mean_sc[[comp]]$pos)== kY), which(colnames(mean_sc[[comp]]$pos)== kX)] <- tmp_cluster_mean_sc %>% filter(cluster == 1) %>% pull(MSC)
        } else{
          mean_sc[[comp]]$pos[ which(rownames(mean_sc[[comp]]$pos)== kY), which(colnames(mean_sc[[comp]]$pos)== kX)] <- 0
        }
        if(!is_empty(tmp_cluster_mean_sc %>% filter(cluster == -1) %>% pull(MSC))){
          mean_sc[[comp]]$neg[ which(rownames(mean_sc[[comp]]$neg)== kY), which(colnames(mean_sc[[comp]]$neg)== kX)] <- tmp_cluster_mean_sc %>% filter(cluster == -1) %>% pull(MSC)
        } else{
          mean_sc[[comp]]$neg[ which(rownames(mean_sc[[comp]]$neg)== kY), which(colnames(mean_sc[[comp]]$neg)== kX)] <- 0
        }
        if(!is_empty(tmp_cluster_mean_sc %>% filter(cluster == 0) %>% pull(MSC))){
          mean_sc[[comp]]$other[ which(rownames(mean_sc[[comp]]$other)== kY), which(colnames(mean_sc[[comp]]$other)== kX)] <- tmp_cluster_mean_sc %>% filter(cluster == 0) %>% pull(MSC)
        } else{
          mean_sc[[comp]]$other[ which(rownames(mean_sc[[comp]]$other)== kY), which(colnames(mean_sc[[comp]]$other)== kX)] <- 0
        }

      } # kY in keepY[[comp]])
    } # end kX in keepX[[comp]]

  } # end comp in 1:ncomp
  return(mean_sc)
} # end tune.spls



deniv_coord_2D_v1 <- function(matrice, coord){
  X_c <- colnames(matrice) %>% as.numeric()
  Y_c <- rownames(matrice) %>% as.numeric()

  X_0 <- which(coord[1] == X_c)
  Y_0 <- which(coord[2] == Y_c)

  matrice <- as.matrix(matrice)
  # coord = vector of length, dim(matrice)
  stopifnot(length(coord) == length(dim(matrice)))

  x.y <- matrice[Y_0, X_0]  # for mor than 2 dim, shoud switch to array
  x_1.y_1 <- matrice[Y_0 - 1, X_0 - 1] # careful of order of row / col

  distance <- sqrt((X_c[X_0] - X_c[X_0-1])^2 + (Y_c[Y_0] - Y_c[Y_0-1])^2)
  deniv <- (x_1.y_1 - x.y) / distance
  return(deniv)
}


.spls.get_inflection_point <- function(tune.spls2.res, comp, contrib){
  C_comp <- tune.spls2.res[[comp]][[contrib]]
  C_comp_tidy <- C_comp %>% rownames_to_column("Y") %>% gather(X, Z, -Y)

  tmp <- vector(mode = "numeric", length = (dim(C_comp)[1]-1) * (dim(C_comp)[2]-1))
  tmp_m <- matrix(nrow = dim(C_comp)[1]-1, ncol = (dim(C_comp)[2]-1))
  k = 0
  for(x in 2:ncol(C_comp)){
    for(y in 2:nrow(C_comp)){
      coord <- c(colnames(C_comp)[x], rownames(C_comp)[y])
      #print(paste(x,y))
      k = k+1
      z <- deniv_coord_2D_v1(C_comp, coord)
      tmp_m[y-1,x-1] <- z
      tmp[k] <- z
    }
  }

  tmp_m %>% as.data.frame() %>%
    set_names(colnames(C_comp)[-length(colnames(C_comp))]) %>% mutate(Y = rownames(C_comp)[-length(rownames(C_comp))])%>%
    gather(X, deniv, -Y) %>% mutate(signif = (abs(deniv) >= 2*sd(deniv) )) %>%
    filter(signif) %>% dplyr::select(X,Y) %>%
    left_join(C_comp_tidy) -> test
  return(test)
}

spls.get_all_inflection_points <- function(tune.spls2.res){
  ncomp = length(tune.spls2.res)
  contrib = c("pos", "neg")
  comp = contribution = X = Y = Z = c()
  for(i in 1:ncomp){
    for(j in contrib){
      res <- .spls.get_inflection_point(tune.spls2.res, i, j)
      comp = c(comp,rep(i, nrow(res)))
      contribution = c(contribution,rep(j, nrow(res)))
      X = c(X, res$X)
      Z = c(Z, res$Z)
      Y = c(Y, res$Y)
    }
  }
  tr <- as.data.frame(list(comp = comp,
                           contrib = contribution,
                           X= as.numeric(as.character(X)),
                           Y=as.numeric(as.character(Y)),
                           Z=as.numeric(as.character(Z))))
  return(tr)
}

spls.get_keepX_keepY<- function(spls.infection){
  wdt <- spls.infection %>% mutate(distance = sqrt(X^2+Y^2))
  wdt_min <- wdt %>% group_by(comp) %>% summarise(min = min(distance))
  wdt %>% left_join(wdt_min, by = c("comp"="comp", "distance"="min"))
  wdt_min %>% left_join(wdt, by = c("comp"="comp", "min"="distance")) %>%
    dplyr::select(comp, X,Y) %>% unique
}


spls.get_cluster <- function(spls.res){

  # X
  pls.clust.X <- spls.res$loadings$X %>% as.data.frame() %>%
    rownames_to_column("molecule") %>%
    gather(comp, value, -molecule) %>%
    filter(value != 0) %>% # sparse
    group_by(molecule) %>% mutate(val_abs = abs(value))

  pls.clust.X.abs <- pls.clust.X %>% summarise(val_abs = max(val_abs))

  pls.clust.X  <- pls.clust.X  %>%
    inner_join(pls.clust.X.abs, by = c("molecule" = "molecule", "val_abs" = "val_abs")) %>%
    dplyr::select(-val_abs) %>%
    mutate(comp = comp %>% str_remove("comp ") %>% as.numeric) %>%
    mutate(cluster = sign(value)*comp) %>% dplyr::select(c(molecule, cluster))%>%
    mutate(block = "X")

  #Y
  pls.clust.Y <- spls.res$loadings$Y %>% as.data.frame() %>%
    rownames_to_column("molecule") %>%
    gather(comp, value, -molecule) %>%
    filter(value != 0) %>% # sparse
    group_by(molecule) %>% mutate(val_abs = abs(value))

  pls.clust.Y.abs <- pls.clust.Y %>% summarise(val_abs = max(val_abs))

  pls.clust.Y  <- pls.clust.Y  %>%
    inner_join(pls.clust.Y.abs, by = c("molecule" = "molecule", "val_abs" = "val_abs")) %>%
    dplyr::select(-val_abs) %>%
    mutate(comp = comp %>% str_remove("comp ") %>% as.numeric) %>%
    mutate(cluster = sign(value)*comp) %>% dplyr::select(c(molecule, cluster))%>%
    mutate(block = "Y")

  # rbind
  pls.clust.final <- rbind(pls.clust.X, pls.clust.Y)

  return(pls.clust.final)
}

spls.plot <- function(spls.Obj, title = "sPLS"){
  cluster.info <- spls.get_cluster(spls.Obj)
  cluster.level <- cluster.info$cluster %>% unique %>% abs %>% sort %>% `*`(c(1,-1))
  X <- spls.Obj$X %>% as.data.frame() %>% rownames_to_column("time") %>% gather(molecule, value, -time)
  Y <- spls.Obj$Y %>% as.data.frame() %>% rownames_to_column("time") %>% gather(molecule, value, -time)

  full <- rbind(X,Y) %>% mutate(time = as.numeric(time)) %>%
    inner_join(cluster.info, by = c("molecule"="molecule")) %>%
    group_by(molecule) %>%
    mutate(cluster = factor(cluster, levels= cluster.level))

  ggplot(data = full, aes(x = time, y = value, group = molecule, col = block)) +
    geom_line() +
    facet_wrap(~ cluster, dir = "v", nrow = 2) +
    ggtitle(title) +
    scale_color_manual(values = color.mixo(1:length(cluster.info$block %>% unique)))
}

