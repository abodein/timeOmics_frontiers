library(tidyverse)
#setwd("/home/antoine/Documents/timeOmics/CS_microbiome/Milk/Final/")

phylochip <- read_tsv("../Data/data_phylochip.txt")
other_metadata <- read_csv("../Data/metadata.csv")

design <- colnames(phylochip) %>% as.data.frame() %>% set_names("sample") %>% filter( sample %>% str_detect("-")) %>%
  mutate(sample = sample %>% str_replace("D-", "Day ")) %>%
  mutate(sample = sample %>% str_replace("Vag", "Vag Day 0")) %>%
  mutate(sample = sample %>% str_replace("Milk", "Milk Day 0")) %>%
  mutate(BABY = sample %>% str_remove(" - .*") %>% factor(levels = unique(.))) %>%
  mutate(time = sample %>% str_remove(".* - ")) %>%
  mutate(time_1 = time %>% str_extract("[0123456789\\.]+") %>% as.numeric()) %>%
  mutate(time_2 = time %>% str_extract("Day|Week|Month|Year")) %>%
  mutate(time_3 = case_when(time_2 == "Day" ~ "1",
                            time_2 == "Week" ~ "7",
                            time_2 == "Month" ~ "30",
                            time_2 == "Year" ~ "365") %>% as.numeric()) %>% 
  mutate(TIME = time_1 * time_3) %>%
  # resolve data type
  mutate(TYPE = case_when(str_detect(sample, "Vag") ~ "Vag_swap",
                          str_detect(sample, "Dad") ~ "Dad_stool",
                          str_detect(sample, "Milk") ~ "Milk",
                          str_detect(sample, "Mom") ~ "Mom_stool",
                          TRUE ~ "Baby_stool"
  ))

# duplicated TIME = 5, antibiotics
# baby = 4,6,8,9,9

de_short <- design %>% filter(TYPE == "Baby_stool") %>% dplyr::select(BABY, TIME) 
#design <- design %>% filter(TYPE == "Baby_stool") %>% .[!(duplicated(de_short)),] 


design <- design %>% filter(TYPE =="Baby_stool") %>%
  mutate(BABY = as.character(BABY) %>% as.numeric()) %>%
  filter(BABY %in% c(2,3,5,7,9,11,13,14, 1,8,12)) %>%
  filter(TIME <= 100) %>%
  dplyr::select(sample, BABY, TIME) %>% 
  left_join(other_metadata) %>%
  mutate(BABY = as.character(BABY))

OTUref <- phylochip[c(1,2)] %>% na.omit %>% mutate(Feature = paste0("F_",Feature))

OTU <- phylochip[-1,-c(2,3)] %>% gather(sample, value, -Feature) %>%  
  mutate(Feature = paste0("F_",Feature)) %>%
  filter(sample %in% design$sample) %>% right_join(design) %>% 
  mutate(s.short = paste0(BABY, "_", TIME)) %>%
  dplyr::select(Feature, s.short, value) %>% 
  spread(Feature, value) %>% as.data.frame() %>% column_to_rownames("s.short")

save(OTU, OTUref, design, file = "./milk_data.RData")


