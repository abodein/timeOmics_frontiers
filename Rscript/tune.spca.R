library(tidyverse)
library(mixOmics)
source("~/Documents/timeOmics/silhouette/silhouette.R")
#library(numDeriv)
library(pspline)


tune.spca <- function(X, ncomp, keepX){
  # keepX doit être une liste de longueur ncomp.
  # Should be executed after a first PCA to get ncomp
  mean_sc <- list()
  DF <- Area_distance(X)
  grid_min_keepX <- unlist(lapply(keepX,min))
  choice.keepX <- list()
  # already_selected_molecules <- vector(mode = "character", length = 0)
  
  for(comp in 1:ncomp){
    # keepx grid
    tmp_keepX <- grid_min_keepX
    
    results_mean_silhouette_coef_pos   <- vector(mode = "double", length= length( keepX[[comp]] ))
    results_mean_silhouette_coef_neg   <- vector(mode = "double", length= length( keepX[[comp]] ))
    results_mean_silhouette_coef_other <- vector(mode = "double", length= length( keepX[[comp]] ))
    
    for(kX in keepX[[comp]]){
      # new keepX grd after each loop
      tmp_keepX[comp] <- kX
      # Do sPCA
      spca.res <- mixOmics::spca(X, ncomp = ncomp, keepX = tmp_keepX)
      tmp_cluster <- spca.res$loadings$X[,comp] %>% sign() %>% as.data.frame() %>% rownames_to_column %>%
        set_names(c("molecule","cluster")) # %>% mutate( cluster = ifelse( molecule %in% already_selected_molecules, cluster, 0))
      # already_selected_molecules <- c(already_selected_molecules, tmp_cluster$mulecule[tmp_cluster$cluster != 0])
      # assign cluster to molecule in DF (pairwaise distance between every molecule)
      tmp_DF_w_clu <- Add_Cluster_metadata(DF, tmp_cluster)
      # compute silhouette coef
      tmp_silhouette <- Slhouette_coef_df(tmp_DF_w_clu)
      # mean by cluster
      tmp_cluster_mean_sc <- tmp_silhouette %>% group_by(cluster) %>% summarise(mean(silhouette.coef))
      results_mean_silhouette_coef_pos[ which(kX == keepX[[comp]])[1]] <- if( is_empty(tmp_cluster_mean_sc %>% filter(cluster == 1) %>% pull(`mean(silhouette.coef)`))){
        #tmp_cluster_mean_sc %>% filter(cluster == 0) %>%  pull(`mean(silhouette.coef)`) # other
        0  ## no molecule selected => mean sc in this cluster = 0 àor -1, why not ? otherwise boost global mean
      } else {
        tmp_cluster_mean_sc %>% filter(cluster == 1) %>%  pull(`mean(silhouette.coef)`)
      }
      results_mean_silhouette_coef_neg[ which(kX == keepX[[comp]])[1]] <- if(is_empty(tmp_cluster_mean_sc %>% filter(cluster == -1) %>% pull(`mean(silhouette.coef)`))){
        #tmp_cluster_mean_sc %>% filter(cluster == 0) %>% pull(`mean(silhouette.coef)`) # other
        0
      } else {
        tmp_cluster_mean_sc %>% filter(cluster == -1) %>%  pull(`mean(silhouette.coef)`)
      }
      results_mean_silhouette_coef_other[ which(kX == keepX[[comp]])[1]] <- tmp_cluster_mean_sc %>% filter(cluster == 0) %>% pull(`mean(silhouette.coef)`)
    }
    
    # attribute results
    mean_sc[[comp]] <- list("pos" = results_mean_silhouette_coef_pos, 
                            "neg" = results_mean_silhouette_coef_neg,
                            "other" = results_mean_silhouette_coef_other)
    choice.keepX[[comp]] <- cbind( keepX[[comp]], rowMeans(cbind(results_mean_silhouette_coef_pos,results_mean_silhouette_coef_neg))) %>% 
      as.tibble %>% set_names("X", "value") %>% top_n(n = 1, wt = value) %>% pull(X) %>% .[1]
  }
  
  results <- list("choice.keepX" = choice.keepX, "mean_sc" = mean_sc, "keepX" = keepX, ncomp = ncomp,
                  "data" = X)
  return(results)
  #return( list())
} 
spca.tuning.plot <- function(tune.spca.Obj, comp = 1){
  keepX_c <- tune.spca.Obj$keepX[[comp]]
  mean_sc_c <- tune.spca.Obj$mean_sc[[comp]] %>% 
    as.tibble %>% mutate(keepX = keepX_c) %>%
    #mutate(mena = cbind(pos,neg,other) %>% rowMeans ) %>% 
    #mutate(mean_pos_neg = cbind(pos,neg) %>% rowMeans ) %>% 
    gather(cluster, value, -keepX)
  ggplot(mean_sc_c, aes(keepX, value, group = cluster, color = cluster)) + geom_line() + theme_minimal()
}

spca.tuning.plot.allcomp <- function(tune.spca.Obj){
  tib <- tibble()
  for(comp in 1:tune.spca.Obj[["ncomp"]]){
    keepX_c <- tune.spca.Obj$keepX[[comp]]
    mean_sc_c <- tune.spca.Obj$mean_sc[[comp]] %>% 
      as.tibble %>% mutate(keepX = keepX_c) %>%
      #mutate(mena = cbind(pos,neg,other) %>% rowMeans ) %>% 
      #mutate(mean_pos_neg = cbind(pos,neg) %>% rowMeans ) %>% 
      gather(cluster, value, -keepX) %>% mutate(comp = comp) 
    #print(ggplot(mean_sc_c, aes(keepX, value, group = cluster, color = cluster)) + geom_line() + theme_minimal())
    tib <- tib %>% rbind(mean_sc_c)
  }
  ggplot(tib, aes(keepX, value, group = cluster, color = cluster)) + geom_line() + facet_wrap(~comp) + theme_bw() + 
    labs(title = "sPCA - keepX optim. using mean silhouette by sign contrib. cluster", ylab = "Mean Silhouette Coef.")
  #return(tib)
}

.tune.spca.optimal.choice.keepX <- function(tune.spca.Obj){
  ncomp <- tune.spca.Obj$ncom
  grid_min_keepX <- list()
  for(comp in 1:ncomp){
    grid_min_keepX[[comp]] <- min(tune.spca.Obj$keepX[[comp]])
  }
  selected_molecules <- list()
  unlist(grid_min_keepX) -> grid_min_keepX
  for(comp in 1:tune.spca.Obj$ncomp){
    keepX_abs <- tune.spca.Obj$keepX[[comp]]
    selected_molecules[[comp]] <- list()
    for(contrib in c("pos","neg")){
      mean_sc <- tune.spca.Obj$mean_sc[[comp]][[contrib]]
      # first derivative  ~ quite useless
      der_1 <- predict(sm.spline(keepX_abs, mean_sc), keepX_abs, 1)  # too smooth
      # second derivative
      der_2 <- predict(sm.spline(keepX_abs, mean_sc), keepX_abs, 2)
      mab <- cbind(keepX_abs, mean_sc, der_1, der_2, abs(der_1), abs(der_2))
      # der_2_max <- mab[which.max(mab[,6]),1]
      # der_1_max <- mab[which.max(mab[,5]),1]
      nb_molecule_group <- detec_peaks(mab[,4], keepX_abs)[1]  #first significative drop
      tmp_keepX <- grid_min_keepX %>% replace(comp, nb_molecule_group)
      spca.res <- mixOmics::spca(tune.spca.Obj$data, ncomp = ncomp, keepX = tmp_keepX)
      tmp_cluster <- spca.res$loadings$X[,comp] %>% sign() %>% as.data.frame() %>% rownames_to_column %>%
        set_names(c("molecule","cluster"))
      if(contrib == "pos"){
        selected_molecules[[comp]][[contrib]] <- tmp_cluster %>% 
          filter(cluster == 1) %>% 
          pull(molecule) 
      } else {  # contrib = neg 
        selected_molecules[[comp]][[contrib]] <- tmp_cluster %>% 
          filter(cluster == -1) %>% 
          pull(molecule) 
      }
    }
  }
  return(selected_molecules)
}

.plot.tune.spca.optimal.choice.keepX<- function(selected_molecules, tune.spca.Obj){
  # to modify S4 class object plot(...)
  data = tune.spca.Obj$data
  # transforme selected moleduls list to cluster (molecule | cluster)
  ta <- imap_dfr(selected_molecules,~ tibble(comp = rep(.y, length(.x$pos)+length(.x$neg)),
                                             molecules = c(.x$pos, .x$neg),
                                             contrib = c(rep("pos", length(.x$pos)), rep("neg", length(.x$neg))))) %>%
    mutate(cluster = paste0(comp,"_",contrib)) %>% dplyr::select(molecules, cluster) %>% mutate(cluster = as.factor(cluster))
  
  # order cluster
  lev <- vector(mode = "character", length = nlevels(ta$cluster))
  for(i in 1:(length(lev)/2)){
    # lev[2*i -1] <- paste0(i,"_pos")  
    # lev[2*i] <- paste0(i,"_neg")
    lev[i + (length(lev)/2)] <- paste0(i,"_neg")
    lev[i] <- paste0(i,"_pos")
  }
  
  tune.spca.Obj$data %>% as.data.frame() %>% rownames_to_column("time") %>% gather(molecules, value, -time) %>%
    filter(molecules %in% ta$molecules) %>% left_join(ta) %>% mutate(time = as.numeric(time)) %>%
    mutate(cluster = factor(cluster,levels = lev)) %>%
    ggplot(aes(x = time, y = value, fill = molecules, col = cluster)) + geom_line() + facet_wrap(~cluster) + theme_bw()
  
}

detec_peaks <- function(X, x.axis){
  # der2
  X.sd <- sd(abs(X))
  X.diff <- diff(X)
  X.1 <- X.diff[1: (length(X.diff)-1)] %>% sign()
  X.2 <- X.diff[2:length(X.diff)] %>% sign()
  peaks <- X[2:length(X.diff)][X.1 != X.2]
  peaks <- peaks[abs(peaks) >= X.sd]  # check if peaks are greater than standard dev.
  x.axis[which(X %in% peaks)]  # what about no peaks is detected ?
}

detec_drop <- function(X, x.axis){  # negative significant slopes
  # compute slope
  deniv <- function(xb, xa, yb, ya){return((yb-ya)/(xb-xa))}
  D <- vector(mode = "numeric", length =length(X)-1)
  for(i in 1 : length(D)){
    D[i] <- (X[i+1]-X[i])/(x.axis[i+1]+x.axis[i])
  }
  c1 <- abs(D) > sd(D[D < 0])  # On regarde parmis les pentes négatives le SD. et on regarde si on a des choses au dessus.
  c2 <- D < 0 # only drop
  return(x.axis[c1 & c2][1]) # first value
}


pca.get_cluster <- function(pca.Obj){
  # should verfied if pca.Obj == pca.Obj
    pca.clust <- pca.Obj$loadings$X %>% as.data.frame() %>% 
      rownames_to_column("molecule") %>% 
      gather(comp, value, -molecule) %>%
      group_by(molecule) %>% mutate(val_abs = abs(value))
    
    pca.clust.abs <- pca.clust %>% summarise(val_abs = max(val_abs))
    pca.clust <- pca.clust %>% 
      inner_join(pca.clust.abs, by = c("molecule" = "molecule", "val_abs" = "val_abs")) %>%
      dplyr::select(-val_abs) %>% 
      mutate(comp = comp %>% str_remove("PC") %>% as.numeric) %>%
      mutate(cluster = sign(value)*comp) %>% dplyr::select(c(molecule, cluster)) 
  return(pca.clust)
}

pca.plot <- function(pca.Obj, title = "PCA"){
  cluster.info <- pca.get_cluster(pca.Obj)
  #cluster.level <- cluster.info$cluster %>% unique %>% abs %>% sort %>% `*`(c(1,-1))
  cluster.level <- .get.cluster.info(cluster.info$cluster)
  X <- pca.Obj$X
  # norm_profile <- function(profile){
  #   profile*length(profile)/(sum(abs(profile), na.rm = T))
  # }
  # X <- lapply(as.data.frame(X), norm_profile) %>% as.data.frame()
  # X <- scale(X, center = T, scale = T)
  
  data <- X %>% as.data.frame() %>% rownames_to_column("time") %>%
    gather(molecule, value, -time) %>%
    left_join(cluster.info, by = c("molecule"="molecule")) %>%
    group_by(molecule) %>%
    mutate(time = as.numeric(time)) %>%
    mutate(cluster = factor(cluster, levels= cluster.level))
  
  ggplot(data = data, aes(x = time, y = value, group = molecule, col = cluster)) + 
    geom_line() + 
    facet_wrap(~ cluster, dir = "v", nrow = 2) +
    ggtitle(title) + 
    scale_color_manual(values = color.mixo(1:length(cluster.level)))
}

norm_profile <- function(profile){
  profile*length(profile)/(sum(abs(profile), na.rm = T))
}

spca.plot <- function(pca.Obj, title = "sPCA"){
  cluster.info <- pca.get_cluster(pca.Obj) %>% filter(cluster != 0)
  # cluster.level <- cluster.info$cluster %>% unique %>% abs %>% sort %>% `*`(c(1,-1))
  cluster.level <- .get.cluster.info(cluster.info$cluster)
  
  X <- pca.Obj$X
  # norm_profile <- function(profile){
  #   profile*length(profile)/(sum(abs(profile), na.rm = T))
  # }
  # X <- lapply(as.data.frame(X), norm_profile) %>% as.data.frame()
  # X <- scale(X, center = T, scale = T)
  
  data <- X %>% as.data.frame() %>% rownames_to_column("time") %>%
    gather(molecule, value, -time) %>%
    left_join(cluster.info, by = c("molecule"="molecule")) %>%
    filter(!is.na(cluster)) %>%  ## filter cluster != 0, NA introduced
    group_by(molecule) %>%
    mutate(time = as.numeric(time)) %>% 
    mutate(cluster = factor(cluster, levels= cluster.level))
  
  ggplot(data = data, aes(x = time, y = value, group = molecule, col = cluster)) + 
    geom_line() + 
    facet_wrap(~ cluster, dir = "v", nrow = 2) +
    ggtitle(title) + 
    scale_color_manual(values = color.mixo(1:length(cluster.level)))
}

.get.cluster.info <- function(cluster) {
  c1 <- cluster %>% abs %>% unique
  cc <- c(c1,c1) %>% sort %>% `*`(c(1,-1))
  cc[cc %in% unique(cluster)]
}


## sparse v2
tune.spca.choice.keepX <- function(tune.spca.Obj, draw = TRUE) {
  ncomp = tune.spca.Obj$ncomp
  choice.keepX = vector(mode = "numeric", length = ncomp)
  #choice.keepX <- tune.spca.Obj$keepX %>% lapply(max) %>%  unlist  # init to max and reduce with beter value
  choice.keepX.all <- list()
  for(k in 1:ncomp){
    V <- tune.spca.Obj$mean_sc[[k]]
    choice.keepX.all[[k]] <- list()
    for(contrib in names(V)){
      choice.keepX.all[[k]][[contrib]] <- detec_drop(V[[contrib]], tune.spca.Obj$keepX[[k]])
    }
    choice.keepX[k] = min(choice.keepX.all[[k]][["pos"]], choice.keepX.all[[k]][["neg"]])
  }
  
  if(draw){
    nrow = 0
    for(k in 1:ncomp){
      nrow = nrow + length(tune.spca.Obj$keepX[[k]])*3
    }
    plot_df <- matrix(ncol = 4, nrow = nrow) %>% as.data.frame() %>% set_names(c("comp", "kX", "contrib", "MSC"))
    r <- 1
    for(k in 1:ncomp){
      for(kX in 1:length(tune.spca.Obj$keepX[[k]])) {
        for( contrib in c("pos", "neg", "other")){
          plot_df[r,] <- c(k, tune.spca.Obj$keepX[[k]][kX], contrib, tune.spca.Obj$mean_sc[[k]][[contrib]][kX])
          r <- r+1
        }
      }
    }
    plot_df <- plot_df %>% 
      mutate(kX = as.numeric(kX)) %>%
      mutate(MSC = as.numeric(MSC)) %>%
      mutate(contrib = factor(contrib, levels = c("pos", "neg", "other")))%>%
      mutate(comp = ifelse(contrib == "other", paste0(comp, "_other"), comp) %>% factor(levels = paste0(rep(1:ncomp, each=2), c("", "_other")) ))  
      # filter(contrib != "other") 
    # plot_df <- split(plot_df, paste(plot_df$comp, plot_df$contrib)) %>% map_df(~mutate(.x, MSC = norm_profile(MSC)))
    gg <- ggplot(data = plot_df, aes(x=kX, y=MSC, group = contrib, col = contrib)) + geom_line() + 
      scale_color_manual(values = color.mixo(1:3)) +
      facet_wrap(.~comp, scales = "free", dir = "v") + 
      ggtitle("Tuning sPCA")
    print(gg)
    
    # ggplot(data = plot_df, aes(x=kX, y=MSC, group = contrib, col = contrib)) + geom_line() + 
    #   scale_color_manual(values = color.mixo(1:3)) +
    #   geom_vline(data = cbind(kX = rep(choice.keepX, each = 2), 
    #                           comp = paste0(rep(1:ncomp, each=2), c("", "_other")),
    #                           contrib = rep("other", ncomp*2),
    #                           MSC = rep(0, ncomp*2)) %>% as.data.frame(),
    #              aes(xintercept = kX, col = contrib)) +
    #   facet_wrap(.~comp, scales = "free", dir = "v")
    #     
  }
  return(choice.keepX)
}

wrapper.silhouette.pca <- function(X, ...){
  X <- as.data.frame(X)
  X.pca <- pca(X = X, ...)
  X.pca.cluster <- pca.get_cluster(X.pca)
  
  X.DF <- Spearman_distance(X)
  X.DF_clu <- Add_Cluster_metadata(X.DF,  X.pca.cluster)
  X.SC <- Slhouette_coef_df(X.DF_clu)
  return(mean(X.SC$silhouette.coef))
}

wrapper.silhouette.spca <- function(X, keepX, plot.t = FALSE, ...){
  X <- as.data.frame(X)
  X.spca <- spca(X = X, keepX = keepX, ...)
  X.spca.cluster <- pca.get_cluster(X.spca) %>% filter(cluster != 0)
  X.filter <- X %>% dplyr::select(X.spca.cluster$molecule)
  
  X.DF <- Spearman_distance(X.filter)
  X.DF_clu <- Add_Cluster_metadata(X.DF,  X.spca.cluster)
  X.SC <- Slhouette_coef_df(X.DF_clu)
  if(plot.t){
    print(plot_silhouette_order_color(X.SC))
  }
  return(mean(X.SC$silhouette.coef))
}

wrapper.silhouette.spca.paper <- function(X, keepX, plot.t = FALSE, ...){
  X <- as.data.frame(X)
  X.spca <- spca(X = X, keepX = keepX, ...)
  X.spca.cluster <- pca.get_cluster(X.spca) %>% filter(cluster != 0)
  X.filter <- X %>% dplyr::select(X.spca.cluster$molecule)
  
  X.DF <- Spearman_distance(X.filter)
  X.DF_clu <- Add_Cluster_metadata(X.DF,  X.spca.cluster)
  X.SC <- Slhouette_coef_df(X.DF_clu)
  if(plot.t){
    print(plot_fig.paper(X.SC))
  }
  return(mean(X.SC$silhouette.coef))
}
