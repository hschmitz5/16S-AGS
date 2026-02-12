# This code only works for genus or species

rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")
library(BacDive)
library(purrr)

fname_out   <- "./data/bacdive/bd_genera.rds"
fname_out2  <- "./data/bacdive/bd_genera_summary.rds"

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
    
    cm <- pluck(.x, "Morphology", "cell morphology", .default = NULL)
    
    # If NULL → no morphology
    if (is.null(cm)) {
      gram_values <- character(0)
      
    } else {
      
      # If it is a single morphology entry (not list of entries)
      if (!is.list(cm[[1]])) {
        cm <- list(cm)   # wrap into list
      }
      
      gram_values <- purrr::map_chr(
        cm,
        ~ pluck(.x, "gram stain", .default = NA_character_)
      ) %>%
        na.omit()
    }
    
    tibble(
      ID = as.character(pluck(.x, "General", "BacDive-ID", .default = NA)),
      gram_stain = list(gram_values),   # <- key change
      query = g
    )
  })
  
  bd_gram_stain <- bind_rows(bd_gram_stain, gram_df)
  Sys.sleep(1)
}

# Majority calculations
bd_gram_stain <- bd_gram_stain %>%
  rowwise() %>%
  mutate(
    majority_gram = if (length(gram_stain) == 0) NA_character_ else names(which.max(table(gram_stain))),
    majority_percent = if (length(gram_stain) == 0) NA_real_ else max(table(gram_stain)) / length(gram_stain) * 100
  ) %>%
  ungroup()

saveRDS(bd_gram_stain, fname_out)

genera_summary <- bd_gram_stain %>%
  # Keep only IDs with a majority_gram (exclude IDs with no data)
  filter(!is.na(majority_gram)) %>%
  
  # Group by genus
  group_by(query) %>%
  summarise(
    counts = list(table(majority_gram)),  # count how many IDs have each majority
    n_IDs = n(),                           # number of IDs contributing
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    majority_gram = names(which.max(counts)),
    majority_percent = as.numeric(max(counts) / n_IDs * 100)
  ) %>%
  select(
    genus = query,
    majority_gram,
    majority_percent
  ) %>%
  ungroup()

saveRDS(genera_summary, fname_out2)
