rm(list = ls())
library(writexl)
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")

fname_out <- "./data/ASV_names.xlsx"

# Define taxa in which at least one sample has abundance > rel_ab_cutoff
high_ab_taxa <- get_rel_ASV(ps) %>%
  filter(Abundance > rel_ab_cutoff) %>%
  distinct(OTU) %>%
  pull(OTU)

# Define DA Taxa
ancom_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel = FALSE)
DA_taxa <- ancom_taxa$high_ab 

taxonomy <- get_taxonomy(ps) %>%
  filter(OTU %in% high_ab_taxa) %>%
  mutate(
    DA = ifelse(OTU %in% DA_taxa, "T", "F")
  ) %>%
  arrange(desc(DA)) %>%
  select(-OTU, -Kingdom, -Class, -Phylum) %>%
  select(rev(names(.))) %>%
  rename(ASV_name = Species_updated)

write_xlsx(taxonomy, path = fname_out)