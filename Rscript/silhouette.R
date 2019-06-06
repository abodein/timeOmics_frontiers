# silhouette.R
library(tidyverse)
library(geiger)

Spearman_distance <- function(df){
  # time X molecule, dataframe with rownames(time)
  df %>% cor(use='pairwise.complete.obs', method='spearman') %>%
    as.data.frame %>%
    rownames_to_column("molecule.x") %>%
    gather(molecule.y, dist, -molecule.x) %>%
    mutate(dist = -dist+1) # -rho + 1  ## Spearman like -> distance  [0:2]
}

Spearman_distance_KA <- function(df){
  # time X molecule, dataframe with rownames(time)
  df %>% cor(use='pairwise.complete.obs', method='spearman') %>%
    as.data.frame %>%
    rownames_to_column("molecule.x") %>%
    gather(molecule.y, dist, -molecule.x) %>%
    mutate(dist = (dist-1)*(dist-1)) # -rho + 1  ## Spearman like -> distance  [0:2]
}


Area_distance <- function(df){
  lapply(df %>% as.data.frame, norm_profile) %>% as.tibble() %>% as.list() %>%
    pairwise.aera.under.curve(df %>% rownames() %>% as.double()) %>% as.data.frame() %>%
    rownames_to_column("molecule.x") %>%
    gather(molecule.y, dist, -molecule.x)
}

Add_Cluster_metadata <- function(Spearman_dist_df, cluster_df ){
  # cluster_df =  molecule | cluster
  Spearman_dist_df %>%
  left_join(cluster_df, by = c("molecule.x" = "molecule")) %>% set_names(c("molecule.x", "molecule.y", "dist", "cluster.x")) %>%
  left_join(cluster_df, by = c("molecule.y" = "molecule")) %>% set_names(c("molecule.x", "molecule.y", "dist", "cluster.x", "cluster.y")) %>%
  filter(molecule.x != molecule.y) # discard distance between same molecules.
}

# min (average inside cluster)
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

# min (average outside cluster)
min_average_outside <- function(df, molecule){
  df %>% # filter(molelcule.x == molecule) %>%
    filter(cluster.x != cluster.y) %>% # outside cluster
    group_by(molecule.x) %>% # 1 molecule
    summarise("average" =mean(dist)) %>%
    top_n(n=-1, average) %>% dplyr::select(average) %>% pull
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

# compute silhoute coef for 1 molecule
compute_molecule_silhouette <- function(sdf, molecule){
  # compute silhouette coef for \molecule
  A = average_inside(df = sdf, molecule = molecule)
  B = min_average_outside(df = sdf, molecule = molecule)
  sc <- silhoutte_coef(A = A, B = B)
  #df <- data.frame(molecule, sc, get_cluster(sdf, molecule)) %>% set_names(c("molecule", "silhouette.coef", "cluster"))
  #return(df)
}

Slhouette_coef_df <- function( spearman_df_w_cluster){
  sdf <- spearman_df_w_cluster %>% split( spearman_df_w_cluster$molecule.x )
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


.plot_silhouette <- function(coef.df){
  cluster_level <- coef.df %>% group_by(cluster) %>% 
    summarise(max = max(silhouette.coef)) %>% 
    arrange(desc(max)) %>% 
    dplyr::select(cluster) %>% pull %>% 
    factor(levels = .) %>% levels()
  
  mean_silhouette = mean( coef.df$silhouette.coef)
  
  C <-  coef.df  %>% 
    mutate(cluster = factor(cluster, levels = cluster_level)) %>% 
    arrange(desc(cluster),silhouette.coef) %>% 
    mutate(molecule = factor(molecule, levels = molecule) ) 
  
  ggplot(C, aes(molecule, silhouette.coef)) + 
    geom_bar(aes(fill = cluster), position = "dodge", stat = "identity") + 
    coord_flip() + geom_hline(aes(yintercept = mean_silhouette)) + 
    labs(title = paste0("Silhouette Graph, mean = ",round(mean_silhouette, digits = 2))) +
    theme_minimal() + 
    theme(axis.text.y=element_blank()) + 
    ylim(-1, 1)
}

.plot_silhouette_mean_median <- function(coef.df){
  cluster_level <- coef.df %>% group_by(cluster) %>% 
    summarise(max = max(silhouette.coef)) %>% 
    arrange(desc(max)) %>% 
    dplyr::select(cluster) %>% pull %>% 
=======
plot_silhouette <- function(coef.df){
#  cluster_level <- coef.df %>% group_by(cluster) %>%
#    summarise(max = max(silhouette.coef)) %>%
#    arrange(desc(max)) %>%
#    dplyr::select(cluster) %>% pull %>%
#    factor(levels = .) %>% levels()

  cluster_level <- coef.df$cluster %>% unique %>% as.numeric() %>% as.data.frame() %>% set_names("cluster") %>%
      mutate(abs = abs(cluster)) %>% arrange(abs, desc(cluster)) %>% pull(cluster) %>%
      as.character()

  mean_silhouette = mean( coef.df$silhouette.coef)

  C <-  coef.df  %>%
    mutate(cluster = factor(cluster, levels = cluster_level)) %>%
    arrange(desc(cluster),silhouette.coef) %>%
    mutate(molecule = factor(molecule, levels = molecule) )

  ggplot(C, aes(molecule, silhouette.coef)) +
    geom_bar(aes(fill = cluster), position = "dodge", stat = "identity") +
    coord_flip() + geom_hline(aes(yintercept = mean_silhouette)) +
    labs(title = paste0("Silhoutte Graph, mean = ",round(mean_silhouette, digits = 2))) +
    theme_minimal() +
    theme(axis.text.y=element_blank()) +
    ylim(-1, 1) +
    scale_fill_manual(values = color.mixo(1:length(cluster_level)))
}

plot_silhouette_mean_median <- function(coef.df){
  cluster_level <- coef.df %>% group_by(cluster) %>%
    summarise(max = max(silhouette.coef)) %>%
    arrange(desc(max)) %>%
    dplyr::select(cluster) %>% pull %>%
    factor(levels = .) %>% levels()

  mean_silhouette = mean( coef.df$silhouette.coef)
  median_silhouette = median( coef.df$silhouette.coef)

  C <-  coef.df  %>%
    mutate(cluster = factor(cluster, levels = cluster_level)) %>%
    arrange(desc(cluster),silhouette.coef) %>%
    mutate(molecule = factor(molecule, levels = molecule) )

  ggplot(C, aes(molecule, silhouette.coef)) +
    geom_bar(aes(fill = cluster), position = "dodge", stat = "identity") +
    coord_flip() + geom_hline(aes(yintercept = mean_silhouette)) +
    coord_flip() + geom_hline(aes(yintercept = median_silhouette), color="blue") +
    labs(title = paste0(c("Silhoutte Graph, mean = ",mean_silhouette))) +
    theme_minimal() +
    theme(axis.text.y=element_blank())
}

plot_silhouette_order_color <- function(coef.df){
  cluster_level <- coef.df %>% group_by(cluster) %>%
    summarise(max = max(silhouette.coef)) %>%
    arrange(desc(max)) %>%
    dplyr::select(cluster) %>% pull %>%
    factor(levels = .) %>% levels()

  mean_silhouette = mean( coef.df$silhouette.coef)

  C <-  coef.df  %>%
    mutate(cluster = factor(cluster, levels = cluster_level)) %>%
    arrange(desc(cluster),silhouette.coef) %>%
    mutate(molecule = factor(molecule, levels = molecule) ) %>% mutate(color = as.character(cluster))
  
  ggplot(C, aes(molecule, silhouette.coef)) + 
    geom_bar(aes(fill = color), position = "dodge", stat = "identity") + 
    coord_flip() + geom_hline(aes(yintercept = mean_silhouette)) + 
    labs(title = paste0("Silhouette Graph, mean = ",round(mean_silhouette, digits = 2))) +
    theme_minimal() + 
    theme(axis.text.y=element_blank()) + 
    scale_fill_manual(values = color.mixo(1:4)) +
    ylim(-0.5,1) +
    guides(fill=guide_legend(title="Cluster"))
}

norm_profile <- function(profile){
  profile*length(profile)/sum(profile, na.rm = T)
}

area.between.curves <- function(P1,P2,time){
  geiger:::.area.between.curves(time, P1, P2, xrange = c(min(time),max(time))) %>% abs()
}

pairwise.aera.under.curve <- function(X, Y)
  # X = list
  # Y = vector of x time
{
  robj <- matrix(ncol=length(X), nrow = length(X)) %>% as.data.frame %>% set_names(names(X)) %>%
    mutate( row = names(X)) %>% column_to_rownames("row")
  for(x in names(X)){
    px <- X[[x]]
    for(y in names(X)){
      py <- X[[y]]
      robj[x,y] <- area.between.curves(px, py, time = Y)
    }
  }
  robj
}

# plot(x=c(1,2,3), y=X$V1, ylim=c(0,10), type = "l")
# lines(x=c(1,2,3), y=X$V2, type = "l")
# lines(x=c(1,2,3), y=X$V3, type = "l")

## WORKFLOW
# DF <- Spearman_distance(df)
# DF2 <- Add_Cluster_metadata(DF, cluster_info)
# DF3 <- Slhouette_coef_df(DF2)
# plot_silhouette(DF3)

plot_curves <- function(raw_df, cluster){
  color_palette <- c("red", "black", "blue")
  raw_df %>% as.tibble %>%
    mutate(time = rownames(mat_mRNA) %>% as.numeric()) %>%
    gather(molecule, value, -time) %>%
    left_join(cluster) %>%   # by molecule
    ggplot(aes(x=time, y=value, group=molecule, color = as.factor(cluster))) +
    geom_line() + facet_wrap(~as.factor(cluster)) + theme(legend.position="none")
}

sPCA_feature_contrib_cluster <- function(spca.res, comp, to_rm=c()){
  c <- mixOmics::selectVar(spca.res, comp = comp)$name  # selected features
  cluster_label <- spca.res$loadings$X[,comp] %>% sign() %>% as.data.frame() %>% rownames_to_column %>%
    set_names(c("molecule","cluster")) %>% mutate( cluster = ifelse( molecule %in% c, cluster, 0))
  if(is_empty(to_rm)){
    return(cluster_label)
  } else {
    return(cluster_label %>% mutate(cluster = ifelse(molecule %in% to_rm, 0, cluster)))
  }
}

plot_curves <- function(raw_df, cluster){
  color_palette <- c("red", "black", "blue")
  raw_df %>% as.tibble %>%
    mutate(time = rownames(raw_df) %>% as.numeric()) %>%
    gather(molecule, value, -time) %>%
    left_join(cluster) %>%   # by molecule
    ggplot(aes(x=time, y=value, group=molecule, color = as.factor(cluster))) +
    geom_line() + facet_wrap(~as.factor(cluster)) + theme(legend.position="none")
}

#######
# Optim sPCA feature selection
#######



### figure paper
plot_fig.paper <- function(coef.df, title = ""){
  cluster_level <- coef.df %>% group_by(cluster) %>%
    summarise(max = max(silhouette.coef)) %>%
    arrange(desc(max)) %>%
    dplyr::select(cluster) %>% pull %>%
    factor(levels = .) %>% levels()

  mean_silhouette = mean( coef.df$silhouette.coef)

  C <-  coef.df  %>%
    mutate(cluster = factor(cluster, levels = cluster_level)) %>%
    arrange(desc(cluster),silhouette.coef) %>%
    mutate(molecule = factor(molecule, levels = molecule) ) %>% mutate(color = as.character(cluster))

  ggplot(C, aes(molecule, silhouette.coef)) +
    geom_bar(aes(fill = color), position = "dodge", stat = "identity") +
    coord_flip() + geom_hline(aes(yintercept = mean_silhouette)) +
    labs(title = paste0(title, "Silhouette Graph, mean = ",round(mean_silhouette, digits = 2))) +
    theme_minimal() +
    theme(axis.text.y=element_blank()) +
    scale_fill_manual(values = color.mixo(1:4)) +
    ylim(-0.5,1)
  #    ylim(min(coef.df$silhouette.coef), 1) +
  #    guides(fill=guide_legend(title="cluster"))
}

plot_fig.paper2 <- function(coef.df, title = ""){
  cluster_level <- coef.df %>% group_by(cluster) %>%
    summarise(max = max(silhouette.coef)) %>%
    arrange(desc(max)) %>%
    dplyr::select(cluster) %>% pull %>%
    factor(levels = .) %>% levels()

  mean_silhouette = mean( coef.df$silhouette.coef)

  C <-  coef.df  %>%
    mutate(cluster = factor(cluster, levels = cluster_level)) %>%
    arrange(desc(cluster),silhouette.coef) %>%
    mutate(molecule = factor(molecule, levels = molecule) ) %>% mutate(color = as.character(cluster))

  ggplot(C, aes(molecule, silhouette.coef)) +
    geom_bar(aes(fill = color), position = "dodge", stat = "identity") +
    coord_flip() + geom_hline(aes(yintercept = mean_silhouette)) +
    labs(title = paste0(title, "Silhouette Graph, mean = ",round(mean_silhouette, digits = 2))) +
    theme_minimal() +
    theme(axis.text.y=element_blank()) +
    scale_fill_manual(values = color.mixo(1:4)) +

    ylim(min(coef.df$silhouette.coef), 1) +
    guides(fill=guide_legend(title="cluster"))
}
