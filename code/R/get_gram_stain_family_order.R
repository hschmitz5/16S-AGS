rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")
library(BacDive)
library(purrr)

fname_out_partial <- "./data/bd_gram_stain_partial.rds"
fname_out_final   <- "./data/bd_gram_stain.rds"

# Choose DA Taxa
ancom_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel = FALSE)
DA_taxa <- ancom_taxa$high_ab 

DA_order <- taxonomy %>%
  filter(Species_updated %in% DA_taxa) %>%
  mutate(
    g_suffix = if_else(
      str_detect(Species_updated, "_o-[0-9]+$"),
      str_replace(Species_updated, "_o-[0-9]+$", ""),
      NA_character_
    )
  ) %>%
  pull(g_suffix) %>%
  na.omit() %>%
  unique()

DA_family <- taxonomy %>%
  filter(Species_updated %in% DA_taxa) %>%
  mutate(
    g_suffix = if_else(
      str_detect(Species_updated, "_f-[0-9]+$"),
      str_replace(Species_updated, "_f-[0-9]+$", ""),
      NA_character_
    )
  ) %>%
  pull(g_suffix) %>%
  na.omit() %>%
  unique()

#### Family

n <- 2
fname <- paste0("./data/Bacdive/", DA_family[n], "_f.csv")

df <- read.csv(
  fname,
  skip = 2,           # skip the first 2 lines of notes
  stringsAsFactors = FALSE
)

# remove anything not matching family
df_filt <- df %>%
  filter(df$family %in% DA_family[n])

ids_int <- as.integer(df_filt$ID)

bd <- open_bacdive("hannah.schmitz@northwestern.edu", "QP4,@,_nunfQEY6")

bd_gram <- map_df(ids_int, ~ {
  fetch(object = bd, ids = .x, fields = "Gram_stain")
})

recs <- fetch(object = bd, ids = df_filt$ID)
