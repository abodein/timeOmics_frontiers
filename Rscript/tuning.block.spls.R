
tune.block.spls <- function(data, list_keepX, ncomp, design){
  
  DF <- lapply(data, Spearman_distance)  # DF$omic_1; DF$omic2, ...
  # All combination of keepX to compute
  tmp_list_keepX <- expand.grid(list_keepX, stringsAsFactors = F,  # should not has string anyway
                                KEEP.OUT.ATTRS = F)  # and remove (useless) attributes
  SC <- matrix(NA, ncol = (length(data) + 2 + 3), nrow = (nrow(tmp_list_keepX)*length(data)*ncomp)) %>% 
    as.data.frame() %>% set_names(c(names(data), "comp", "omic", "c_pos", "c_neg", "c_other"))
  compteur <- 0
  
  for(comp in 1:ncomp){
    # print(comp)
    for( ligne in 1:nrow(tmp_list_keepX)){
      # 1) identified keepX to compute
      kX <- tmp_list_keepX[ligne,] %>% as.list()
      # print(kX)
      # 2) run block.spls w selected keepX
      block.spls.res <- mixOmics::block.spls(X = data, indY = 1,  # To be verified
                                             keepX = kX, ncomp = ncomp,
                                             design=design, mode = "canonical")
      # 3) by omic clusters : extract contrib.   (pos, neg, others)
      for( omic in names(data)){
        tmp_cluster <- block.spls.res[["loadings"]][[omic]][,comp] %>% sign() %>% as.data.frame() %>% rownames_to_column %>%
          set_names(c("molecule", "cluster"))
        # 4) compute silhouette  (ADD metadata to SIlhouette preformated object + Mean SC)
        tmp_DF_w_clu <- Add_Cluster_metadata(DF[[omic]], tmp_cluster)
        tmp_silhouette <- Slhouette_coef_df(tmp_DF_w_clu)
        tmp_cluster_mean_sc <- tmp_silhouette %>% group_by(cluster) %>% summarise(MSC = mean(silhouette.coef))
        # store information
        contrib.pos <- tmp_cluster_mean_sc %>% filter(cluster == 1) %>% pull(MSC) %>% ifelse(test = is_empty(.), yes = 0, no = .)
        contrib.neg <- tmp_cluster_mean_sc %>% filter(cluster == -1) %>% pull(MSC) %>% ifelse(test = is_empty(.), yes = 0, no = .)
        contrib.other <- tmp_cluster_mean_sc %>% filter(cluster == 0) %>% pull(MSC) %>% ifelse(test = is_empty(.), yes = 0, no = .)
        # 5) Store results
        SC[compteur <- compteur+1,] <- c(kX %>% as.data.frame() %>% dplyr::select(names(data)) %>% unlist() %>% as.vector(), 
                                         comp, omic, contrib.pos, contrib.neg, contrib.other)
      }
    }
  }
  return(SC)
}


all_neighbours <- function(all_position){
  
  neihbourhood <- function(all_position, set_p){
    # from a set, find in data, all the set[][n+1] combinaision
    index <- imap(set_p, ~ which(all_position[[.y]] == .x))
    # find [index + 1] // all combinaison
    index_2 <- list()
    value <- list()
    for( x in names(set_p)){
      value[[x]] = list()
      i <- (index[x] %>% unlist() %>% as.vector()) 
      idx_0 <- all_position[[x]][i]
      idx_1 <- all_position[[x]][i +1]
      #value[[x]] <- ifelse(is.na(idx), all_position[[x]][(index[x] %>% unlist() %>% as.vector())], idx)
      value[[x]][[x]] <- ifelse(is.na(idx_1), idx_0, idx_1)
      for(y in names(set_p)[-match(x,names(set_p))]){
        value[[x]][[y]] <- all_position[[y]][(index[y] %>% unlist() %>% as.vector())]
      }
    }
    #index_2 <- imap(value, ~ which(all_position[[.y]] == .x))
    return(value)
  }
  
  all_set_p <- expand.grid(all_position, stringsAsFactors = F,  # should not has string anyway
                           KEEP.OUT.ATTRS = F)
  
  AAA <- matrix(NA, ncol = 2*ncol(all_set_p), nrow= nrow(all_set_p)*ncol(all_set_p)) %>%
    as.data.frame() %>% 
    purrr::set_names(names(all_set_p), paste0(names(all_set_p), ".1"))
  w <- 1
  for(i in 1:nrow(all_set_p)){
    #AAA[[x]] <- neihbourhood(all_position, all_set_p[x,])
    a <- neihbourhood(all_position, all_set_p[i,])
    for(z in names(a)){
      for(j in names(a)){
        AAA[w, match(j, names(AAA))] <- all_set_p[i, match(j, names(all_position))]
        AAA[w, match(paste0(j, ".1"), names(AAA))] <- a[[z]][[j]]
      }
      w <- w + 1
    }
  }
  return(AAA)
}

remove_duplicate_neighbours <- function(AAA, first_name){
  X <- apply( AAA[ , first_name ] , 1 , paste , collapse = "-" )
  Y <- apply( AAA[ , paste0(first_name,'.1') ] , 1 , paste , collapse = "-" )
  return(X != Y)
}


get_direction <- function(AAA, first_name){
  X <- AAA[ , first_name ]
  Y <- AAA[ , paste0(first_name,'.1') ]
  return(apply( X!=Y, 1, function(x) names(data)[x]) %>% as.character())
}

get_origin <- function(x) lapply(x, min)

get_distance_w_origin <- function(DF, first_name, origin){
  u <- lapply(DF[, first_name], function(x) as.character(x) %>% as.numeric) 
  u <- purrr::map(seq_along(u[[1]]), function(x) map_dbl(u, `[`,x))
  r <- lapply(u, norme,y = origin) %>% unlist()
  return(r)
}


norme <- function(x,y){
  # x and y are numeric coordinate vectors and have both the same size
  stopifnot(length(x) == length(y))
  stopifnot(is.numeric(x) && is.numeric(y))
  stopifnot(all(names(x) %in% names(y)))
  y <- y[names(x)]  # same things
  return(sqrt(sum((x - y)^2)))
}

deniv <- function(A.coord, B.coord, A.value, B.value){
  # A.coord and B.coord are numeric coordinate vectors and have both the same size
  # A.coord, B.coord, A.value, B.value are numeric
  stopifnot(length(A.coord) == length(B.coord))
  stopifnot(is.numeric(A.coord) && is.numeric(B.coord))
  stopifnot(is.numeric(A.value) && is.numeric(B.value))
  return( (B.value - A.value)/norme(A.coord, B.coord) )
}

wrapper.deniv <- function(df.A.coord, df.B.coord, vc.A.value, vc.B.value){
  # TODO check dim and lenght
  # col order
  den <- numeric(length = length(vc.A.value))
  for(i in 1:nrow(df.A.coord)){  # == nrow df.B.coord, == length valueS
    A.coord <- df.A.coord[i,] %>% unlist()
    B.coord <- df.B.coord[i,] %>% unlist() %>% .[paste0(names(A.coord), ".1")] %>% `names<-`(names(A.coord))
    A.value <- vc.A.value[i]
    B.value <- vc.B.value[i]
    den[i] <- deniv(A.coord, B.coord, A.value, B.value)
  }
  return(den)
}
#undebug(wrapper.deniv)

detect_drop_out <- function(DF, comp, block, contrib = c("pos", "neg", "other"), first_name){
  contrib = match.arg(contrib)
  a <- DF %>% ungroup %>% 
    .[.[["comp"]] == comp,] %>% 
    .[.[["omic"]] == block,] %>%
    .[.[[paste0("signif_", contrib)]],] %>% 
    .[.[["distance_from_origin"]] == min(.[["distance_from_origin"]]),] %>%
    .[1, first_name] %>% unlist
  return(a)
}

compute_slopes <- function(tuning.obj, data){
  
  all_position <- lapply(tuning.obj[names(data)], unique) %>% lapply(as.numeric) %>% lapply(sort) # keepX
  
  AAA <- all_neighbours(all_position = all_position) %>% lapply(as.character) %>% 
    as.data.frame() #%>% unique %>% .[ remove_duplicate_neighbours(., names(data)),]
  
  full_tuning <- AAA %>% full_join(tuning.obj) %>%  # to get every 'first' points contribution
    rename_at(vars(c_pos, c_neg, c_other), ~paste0(c("c_pos", "c_neg", "c_other"), ".bis")) %>%  # rename 'first' point contribution
    rename_at(vars(names(data)), ~ paste0(names(data), ".bis")) %>%  # rename 'first' point
    rename_at(vars(paste0(names(data),".1")), ~names(data)) %>%  # rename neigbourhood in 'first' point
    left_join(tuning.obj) %>%  # get neighboorhoud contrib.
    ### WARNINGS
    rename_at(vars(names(data)), ~ paste0(names(data), ".1")) %>%  # rename neigbourhood
    rename_at(vars(paste0(names(data), ".bis")), ~ names(data)) %>%   # rename 'first' point
    rename_at(vars(c_pos, c_neg, c_other), ~paste0(c("c_pos", "c_neg", "c_other"), ".1")) %>%  # rename neigbourhoud contrib
    rename_at(vars(paste0(c("c_pos", "c_neg", "c_other"), ".bis")), ~c("c_pos", "c_neg", "c_other")) %>%  # rename 'first' point contrib
    mutate_at(vars(-omic),as.numeric) %>%
    mutate(slope_pos = wrapper.deniv(df.A.coord = dplyr::select(., names(data)),
                                     df.B.coord = dplyr::select(., paste0(names(data), ".1")),
                                     vc.A.value = c_pos,
                                     vc.B.value = c_pos.1)) %>%
    mutate(slope_neg = wrapper.deniv(df.A.coord = dplyr::select(., names(data)),
                                     df.B.coord = dplyr::select(., paste0(names(data), ".1")),
                                     vc.A.value = c_neg,
                                     vc.B.value = c_neg.1)) %>%
    mutate(slope_other = wrapper.deniv(df.A.coord = dplyr::select(., names(data)),
                                       df.B.coord = dplyr::select(., paste0(names(data), ".1")),
                                       vc.A.value = c_other,
                                       vc.B.value = c_other.1)) %>%
    na.omit() # %>%# remove duplicates
  #unique %>%  .[ remove_duplicate_neighbours(., names(data)),]
  return(full_tuning)
}

detect_signif_slopes <- function(slopes.obj, data, all_position){
  full_tuning <- slopes.obj
  
  all_position <- lapply(full_tuning[names(data)], unique) %>% lapply(as.numeric) %>% lapply(sort)
  
  SD <- full_tuning %>% mutate(direction = get_direction(., names(data))) %>%
    group_by(comp, omic, direction) %>% 
    summarise(sd_pos = sd(slope_pos), sd_neg = sd(slope_neg), sd_other = sd(slope_other)) %>%
    right_join(full_tuning %>% mutate(direction = get_direction(., names(data)))) %>% 
    # mutate(signif_pos = slope_pos > sd_pos) %>%
    # mutate(signf_neg = slope_neg > sd_neg) %>%
    # mutate(signf_other = slope_other > sd_other) %>% 
    ungroup %>%
    mutate(distance_from_origin = get_distance_w_origin(., names(data), origin= get_origin(all_position)%>% unlist)) %>%
    group_by(comp, omic, direction)
  return(SD)
}

detect_drop_out <- function(DF, comp, block, contrib = c("pos", "neg", "other"), first_name, treshold = 2){
  
  # define short function
  # sd2 <- function(x) 2*sd(x)
  # sd3 <- function(x) 3*sd(x)
  
  contrib = match.arg(contrib)
  # treshold = match.arg(treshold)
  # sdi = match.fun(treshold)
  
  a <- DF %>% ungroup %>% 
    .[.[["comp"]] == comp,] %>% 
    .[.[["omic"]] == block,] # %>%  ## don't change anything for few omic, first approche more sensible
    # mutate(stddev = sd(slope_pos))
  
  b <- a[a[[paste0("slope_", contrib)]] *-1 >= a[[paste0("sd_", contrib)]] * treshold,] %>%  # signif delta  - deniv(B,A) >= 2sd
    .[.[["distance_from_origin"]] == min(.[["distance_from_origin"]]),] %>%
    .[1, c("comp", "omic",first_name, "c_pos", "slope_pos", "c_neg", "slope_neg")] #%>% unlist
  
  # b <- a[a[[paste0("slope_", contrib)]] *-1 >= a[["stddev"]] * treshold,] %>%  # signif delta  - deniv(B,A) >= 2sd
  #   .[.[["distance_from_origin"]] == min(.[["distance_from_origin"]]),] %>%
  #   .[1, first_name] %>% unlist
  
  return(b)
}

wrapper.detect_drop_out <- function(DF, first_name, threshold = 2){
  drop = detect_drop_out(DF, DF[["comp"]][1], DF[["omic"]][1], contrib = c("pos", "neg", "other"), first_name, threshold) %>%
    t %>% as.data.frame() %>% rownames_to_column() %>% set_names(c('rowname',"base"))

  for( comp in DF[["comp"]] %>% unique %>% sort){
    for( block in  DF[["omic"]] %>% unique %>% sort){
      for(contrib in c("pos", "neg")){
        a <- detect_drop_out(DF, comp, block, contrib = contrib, first_name, threshold) %>%
          t %>% as.data.frame() %>% rownames_to_column() %>% set_names(c("rowname", paste(c(comp,block,contrib), collapse = "-")))
        drop <- left_join(drop, a, by=c("rowname" = "rowname"))
      }
    }
  }
  drop <- drop %>% dplyr::select(-base) %>% column_to_rownames("rowname") %>% t %>% as.data.frame() %>%
    rownames_to_column("contrib") %>% mutate(contrib = contrib %>% str_replace(".*-",""))
  return(drop)
}




