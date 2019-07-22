
loadings.get_cluster <- function(mixOmics.res, sparse = F, cutoff = 0){
  # if sparse == TRUE remove molecule with cluster == 0.
  # X

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

  X.block <- do.call("rbind",mixOmics.res$loadings)

  block.info <- imap(mixOmics.res$loadings,
                       function(x,y) rownames(x) %>%
                           as.data.frame %>%
                           set_names("molecule") %>%
                           mutate("block" = y))
  block.info <- do.call("rbind",block.info)


  clust.X <- X.block %>% as.data.frame() %>%
    rownames_to_column("molecule") %>%
    gather(comp, value, -molecule) %>%
    mutate(comp = comp %>% str_replace("comp ", "PC")) %>%
    mutate(comp = comp %>% str_replace("comp", "PC")) %>%
    filter_if_sparse(sparse) %>%
    group_by(molecule) %>% mutate(val_abs = abs(value))

  clust.X.abs <- clust.X %>% dplyr::summarise(val_abs = max(val_abs))

  clust.X.2  <- clust.X  %>%
    inner_join(clust.X.abs, by = c("molecule" = "molecule", "val_abs" = "val_abs")) %>%
    filter_cutoff_abs_value(cutoff) %>%
    dplyr::select(-val_abs) %>%
    mutate(comp = comp %>% str_remove("PC") %>% as.numeric) %>%
    mutate(cluster = sign(value)*comp) %>% dplyr::select(c(molecule, cluster))

  clust.X.2 <- suppressWarnings(clust.X.2 %>% left_join(block.info, by = c("molecule"= "molecule")))
  clust.X.2<-  clust.X.2 %>% mutate(comp = abs(cluster)) %>%
      mutate(contrib = sign(cluster))

  return(clust.X.2)
}
