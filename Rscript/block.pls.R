spls.plot <- function(spls.Obj, title = "sPLS"){
    cluster.info <- loadings.get_cluster(spls.Obj)
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

plot.block.spls <- function(block.spls.Obj, title = "block.sPLS"){
    stopifnot(is(block.spls.Obj, c("block.spls", "block.pls")))
    stopifnot(is(title, "character"))
    
    # get cluster
    cluster.info <- loadings.get_cluster(block.spls.Obj) %>%
        mutate(comp = paste0("Comp ", comp))
    cluster.level <- cluster.info %>% pull(cluster) %>% abs %>% unique %>% 
        sort %>% rep(each=2) %>% `*`(c(1,-1)) 
    
    
    # raw data in block.spls.Obj$X
    X <- lapply(block.spls.Obj$X, function(x) as.data.frame(t(x))) %>% 
        purrr::reduce(rbind) %>% t
    # need to cast in dataframe because rbind check colnames consistensy 
    # and cbind does not
    
    # gather, add cluster info
    X.gather <- X %>% as.data.frame %>% rownames_to_column("time") %>%
        mutate(time = as.numeric(time)) %>%
        gather(molecule, value, -time) %>%
        left_join(cluster.info, by = c("molecule" = "molecule")) %>%
        mutate(cluster = factor(cluster, levels = cluster.level))
    
    ggplot(X.gather, aes(time, value, col = block, group = molecule)) + geom_line() +
        facet_grid(contrib~comp) +
        scale_color_manual(values = color.mixo(1:length(X.gather$block %>% unique))) + 
        ggtitle(title)
}

    
