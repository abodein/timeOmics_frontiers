# graphlan bioreactor

library(tidyverse)
library(mixOmics)

# bioreactor <- readxl::read_excel("~/Documents/timeOmics/Examples/propr_biowaste (1).xls") %>%
#     dplyr::select(Taxonomy, `Cluster ID`) %>% rename(Cluster= `Cluster ID`) %>%
#     na.omit()

#save(bioreactor, file = "../Data/bior_spls.RData")
load("../Data/bior_spls.RData")

# taxo tree
bioreactor <- bioreactor %>% mutate(Taxonomy = Taxonomy %>% str_replace_all(";", "\\.")) %>%
    mutate(Taxonomy = Taxonomy %>% str_replace("sp\\.", "sp")) #%>%
    #mutate(Taxonomy = Taxonomy %>% str_remove("\\.unknown.*")) %>%
    #unique
bioreactor %>% pull(Taxonomy) %>% as.data.frame %>% write_csv('~/ls3x/graphlan/tree.txt', quote = F, col_names = F)

bioreactor %>% dplyr::select(Taxonomy) %>%
    mutate(option = "clade_marker_size") %>%
    mutate(value = "100") %>%
    write_tsv("~/ls3x/graphlan/annot_size_leaf.txt",  col_names = F)

bioreactor %>% dplyr::select(Taxonomy, Cluster) %>%
    mutate(Cluster = factor(Cluster) %>% color.mixo()) %>%
    mutate(option = "annotation_background_color") %>%
    dplyr::select(Taxonomy, option, Cluster) %>%
    write_tsv("~/ls3x/graphlan/annot_color_cluster.txt",  col_names = F)

# leaf color
bioreactor %>% mutate(Cluster = factor(Cluster) %>% color.mixo()) %>%
    mutate(option = "clade_marker_color") %>%
    dplyr::select(Taxonomy, option, Cluster) %>%
    write_tsv('~/ls3x/graphlan/annot.txt', col_names = F)

# annotation color
annot <- bioreactor %>%
    mutate(L1 = Taxonomy %>% str_split("\\.") %>% map_chr(~.x[1])) %>%
    mutate(L2 = Taxonomy %>% str_split("\\.") %>% map_chr(~.x[2])) %>%
    mutate(L1_L2 = paste(L1,L2, sep = "."))

annot %>%    dplyr::select(L1_L2) %>% mutate(option = "annotation_background_color") %>%
    mutate(color = "#DBDBDB") %>%
    write_tsv("~/ls3x/graphlan/annot_background_1.txt",  col_names = F)

# annotation
annot %>% group_by(L1_L2) %>% filter(n() == 1)
annot %>% mutate(option = "annotation") %>%
    dplyr::select(L1_L2, option, L2) %>%
    write_tsv("~/ls3x/graphlan/annot_L2.txt",  col_names = F)

