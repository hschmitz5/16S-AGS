rm(list = ls())
source("./code/R/01_load_ps.R")

sam_name_full <- c("40A", "40B", "40C", "20A", "20B", "20C", "14A", "14B", "14C", 
                   "10A", "10B", "10C", "7A", "7B", "7C", "5A", "5B", "5C")

# load phyloseq object for all sample sizes
ps_full <- readRDS("./data/ps_genus_full.rds") 


rel_wide <- get_rel(ps_full) %>%
  dplyr::select(Genus, Sample, Abundance) %>%  
  pivot_wider(
    names_from = Sample,
    values_from = Abundance
  ) %>%
  dplyr::select(Genus, all_of(sam_name_full)) %>%
  filter(str_detect(Genus, "Unk")) %>%
  column_to_rownames(var = "Genus")

test <- tibble(
  Sample = names(colSums(rel_wide)),
  sum_na = colSums(rel_wide)
)

