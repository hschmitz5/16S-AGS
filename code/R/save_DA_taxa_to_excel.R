rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")
library(writexl)

fname_excel <- "./data/DA_taxa.xlsx"

# Choose DA Taxa
ancom_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel = FALSE)
DA_taxa <- ancom_taxa$high_ab 

DA_genera <- taxonomy %>%
  filter(Species_updated %in% DA_taxa) %>%
  mutate(
    Species_short = sub("_(?=[^_]*$).*", "", Species_updated, perl = TRUE)
  )

write_xlsx(taxa_names, path = fname_excel)