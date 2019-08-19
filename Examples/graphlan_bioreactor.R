# graphlan bioreactor
library(tidyverse)
library(mixOmics)

bioreactor <- readxl::read_excel("~/Documents/timeOmics/Examples/propr_biowaste (1).xls") %>%
    dplyr::select(Taxonomy, `Cluster ID`) %>% rename(Cluster= `Cluster ID`) %>%
    na.omit()

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


### c-section
csec.spca <- read_csv("~/Documents/timeomics_analysis/CS_microbiome/Milk/Final/C-scection_spca.csv") %>%
    dplyr::select(molecule, cluster)
csec.pca <- read_csv("~/Documents/timeomics_analysis/CS_microbiome/Milk/Final/c-section_pca_all.csv")
csec.pca <- csec.pca %>% mutate(taxo = taxo %>% str_remove_all("D_.__") %>% str_replace("sp\\.", "sp") %>% str_replace_all(";","\\."))  %>%
    arrange(taxo) %>%
    mutate(cluster = factor(cluster, levels = c(1,-1,2,-2))) %>%
    mutate(selected = molecule %in% csec.spca$molecule) %>% na.omit() %>%
    mutate(taxo = str_replace(taxo, "Clostridium sensu stricto .", "Clostridium") %>%
               str_replace("Clostridiaceae 1", "Clostridiaceae") %>%
               str_replace(" Clostridium", "Clostridium")) %>%
    filter(molecule != "F_2072")

csec.pca <- csec.pca %>% mutate(taxoo = str_count(taxo,"\\.")) %>%
    mutate(taxo = paste0(taxo,str_dup(".unknown", 6-taxoo)))

csec.pca <- csec.pca %>%
    mutate(L1 = str_split(taxo, "\\.") %>% map_chr(~.x[1])) %>%
    mutate(L2 = str_split(taxo, "\\.") %>% map_chr(~.x[2])) %>%
    mutate(L3 = str_split(taxo, "\\.") %>% map_chr(~.x[3])) %>%
    mutate(L4 = str_split(taxo, "\\.") %>% map_chr(~.x[4])) %>%
    mutate(L5 = str_split(taxo, "\\.") %>% map_chr(~.x[5])) %>%
    mutate(L6 = str_split(taxo, "\\.") %>% map_chr(~.x[6]))

csec.pca %>% write_tsv("~/ls3x/graphlan/csection/csection_pca.tsv", col_names = F)
#csec.pca %>% dplyr::select(taxo) %>% write_tsv("~/ls3x/graphlan/csection/tree.txt", col_names = F)


# annotation L2
csec.pca %>% mutate(L1_L2 = paste(L1,L2, sep = ".")) %>% dplyr::select(L1_L2) %>% unique %>%
    filter(!str_detect(L1_L2, "unknown")) %>%
    mutate(option = "annotation") %>%
    mutate(annot = str_remove(L1_L2, ".*\\.")) %>%
    write_tsv("~/ls3x/graphlan/csection/annot_L2.tsv", col_names = F)

# annotation bacground_color L2
csec.pca %>% mutate(L1_L2 = paste(L1,L2, sep = ".")) %>% dplyr::select(L1_L2) %>% unique %>%
    mutate(option = "annotation_background_color") %>%
    mutate(color = "#E8E8E8") %>%
    write_tsv("~/ls3x/graphlan/csection/annot_bacground_color_L2.tsv", col_names = F)

# annotation bacground_color L3
csec.pca %>% mutate(L1_L3 = paste(L1,L2,L3, sep = ".")) %>% dplyr::select(L1_L3) %>% unique %>%
    filter(!str_detect(L1_L3, "unknown")) %>%
    mutate(option = "annotation_background_color") %>%
    mutate(color = "#DADADA") %>%
    write_tsv("~/ls3x/graphlan/csection/annot_bacground_color_L3.tsv", col_names = F)

# annotation bacground_color L4
csec.pca %>% mutate(L1_L4 = paste(L1,L2,L3,L4, sep = ".")) %>% dplyr::select(L1_L4) %>% unique %>%
    filter(!str_detect(L1_L4, "unknown")) %>%
    mutate(option = "annotation_background_color") %>%
    mutate(color = "#CCCCCC") %>%
    write_tsv("~/ls3x/graphlan/csection/annot_bacground_color_L4.tsv", col_names = F)

# annotation bacground_color L5
csec.pca %>% mutate(L1_L5 = paste(L1,L2,L3,L4,L5, sep = ".")) %>% dplyr::select(L1_L5) %>% unique %>%
    filter(!str_detect(L1_L5, "unknown")) %>%
    mutate(option = "annotation_background_color") %>%
    mutate(color = "#BEBEBE") %>%
    write_tsv("~/ls3x/graphlan/csection/annot_bacground_color_L5.tsv", col_names = F)

# annotation bacground_color L6
csec.pca %>% mutate(L1_L6 = paste(L1,L2,L3,L4,L5,L6, sep = ".")) %>% dplyr::select(L1_L6) %>% unique %>%
    filter(!str_detect(L1_L6, "unknown")) %>%
    mutate(option = "annotation_background_color") %>%
    mutate(color = "#B0B0B0") %>%
    write_tsv("~/ls3x/graphlan/csection/annot_bacground_color_L6.tsv", col_names = F)

# selected star
csec.pca %>% filter(selected) %>% dplyr::select(taxo) %>%
    mutate(option = "clade_marker_shape") %>%
    mutate(color = "*") %>%
    write_tsv("~/ls3x/graphlan/csection/annot_clade_marker_shape.tsv", col_names = F)

# selected size
csec.pca %>% filter(selected) %>% dplyr::select(taxo) %>%
    mutate(option = "clade_marker_size") %>%
    mutate(color = "400") %>%
    write_tsv("~/ls3x/graphlan/csection/annot_clade_marker_size_1.tsv", col_names = F)

# color leaf
csec.pca %>%  mutate(color = color.mixo(cluster)) %>% mutate(option = "annotation_background_color") %>%
    dplyr::select(taxo, option, color) %>%
    write_tsv("~/ls3x/graphlan/csection/annot_cluster_color.tsv", col_names = F)

csec.pca %>%  mutate(color = color.mixo(cluster)) %>% mutate(option = "clade_marker_color") %>%
    dplyr::select(taxo, option, color) %>%
    write_tsv("~/ls3x/graphlan/csection/annot_cluster_color_leaf.tsv", col_names = F)
