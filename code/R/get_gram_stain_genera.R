# This code only works for genus or species

rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")
library(BacDive)
library(purrr)

fname_out_partial <- "./data/bd_gram_stain_partial.rds"
fname_out_final   <- "./data/bd_genera.rds"

# Choose DA Taxa
ancom_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel = FALSE)
DA_taxa <- ancom_taxa$high_ab 

DA_genera <- taxonomy %>%
  filter(Species_updated %in% DA_taxa) %>%
  mutate(
    g_suffix = if_else(
      str_detect(Species_updated, "_g-[0-9]+$"),
      str_replace(Species_updated, "_g-[0-9]+$", ""),
      NA_character_
    )
  ) %>%
  pull(g_suffix) %>%
  na.omit() %>%
  unique()

#####

bd <- open_bacdive("hannah.schmitz@northwestern.edu", "QP4,@,_nunfQEY6")

safe_retrieve <- function(g, max_tries = 3, wait = 2) {
  for (i in seq_len(max_tries)) {
    res <- tryCatch(
      retrieve(object = bd, query = g, search = "taxon"),
      error = function(e) NULL
    )
    
    if (!is.null(res)) return(res)
    
    message("  ⏳ Retry ", i, " for ", g)
    Sys.sleep(wait)
  }
  NULL
}

bd_gram_stain <- tibble()

for (g in DA_genera) {
  
  message("Processing: ", g)
  
  records <- safe_retrieve(g)
  
  if (is.null(records)) {
    message("  ❌ Giving up on ", g)
    next
  }
  
  gram_df <- map_dfr(records, ~ {
    tibble(
      ID = as.character(pluck(.x, "General", "BacDive-ID", .default = NA)),
      gram_stain = pluck(
        .x,
        "Morphology", "cell morphology", "gram stain",
        .default = NA_character_
      ),
      query = g
    )
  })
  
  bd_gram_stain <- bind_rows(bd_gram_stain, gram_df)
  Sys.sleep(1)
  
  saveRDS(bd_gram_stain, fname_out_partial)
}

saveRDS(bd_gram_stain, fname_out_final)

gram_counts <- bd_gram_stain %>%
  filter(!is.na(gram_stain)) %>%
  group_by(query, gram_stain) %>%
  summarise(n = n(), .groups = "drop")
