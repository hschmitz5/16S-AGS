# csv files acquired from advanced search on Bac Dive
rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")
library(BacDive)
library(purrr)
library(tibble)

parent_fname <- "./data/Bacdive/"
fname_out <- "./data/bd_family.rds"

# length of chunk
L <- 100 

# Choose DA Taxa
ancom_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel = FALSE)
DA_taxa <- ancom_taxa$high_ab 

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

#### Functions

load_bd_csv <- function(taxon) {
  fname <- paste0(parent_fname, taxon, "_f.csv")
  
  if (!file.exists(fname)) {
    message("File not found for family: ", taxon, " — skipping")
    return(NULL)  
  }
  
  message("Processing family: ", taxon)
  
  # read in data
  df_raw <- read.csv(
    fname, 
    skip = 2, # skip header lines
    stringsAsFactors = FALSE
  )
  
  # filter only the current family
  df <- df_raw %>%
    filter(family %in% taxon) %>%
    arrange(ID)
  
  removed <- nrow(df_raw) - nrow(df)
  cat("Rows removed:", removed, "\n")
  
  return(df)
}

# define safe fetch function
safe_fetch <- function(id, max_tries = 5, wait = 2) {
  for (i in seq_len(max_tries)) {
    res <- tryCatch(
      fetch(object = bd, ids = id),
      error = function(e) {
        message("    Failed for ID ", id, ": ", e$message)
        return(NULL)
      }
    )
    
    if (!is.null(res)) return(res)
    
    message("  ⏳ Retry ", i, " for ", id)
    Sys.sleep(wait)
  }
  NULL
}

bd <- open_bacdive("hannah.schmitz@northwestern.edu", "QP4,@,_nunfQEY6")

gs_list <- list()  # list to hold results per DA_family

for (taxon in DA_family) {
  # load Bac Dive csv file
  df <- load_bd_csv(taxon)
  if (is.null(df)) next  # skip to next family
  
  # df$ID will be broken into chunks
  total_strains <- length(df$ID)
  num_chunks <- ceiling(total_strains/L)
  
  # Initialize rows list to collect results
  rows <- vector(mode = "list", length = total_strains)
  row_idx <- 1
  
  for (i in seq_len(num_chunks)) {
    # define indices for strains in chunk
    start_idx <- 1 + L*(i-1)
    end_idx <- min(L*i, total_strains)
    chunk_ids <- as.character(df$ID[start_idx:end_idx])
    
    message("Processing chunk: ", i)
    batch <- safe_fetch(chunk_ids)
    if (is.null(batch)) next
    
    for (strain_id in chunk_ids) {
      morph <- batch$results[[strain_id]]$Morphology
      
      # CASE 1: no morphology → NA row
      if (is.null(morph) || length(morph) == 0) {
        rows[[row_idx]] <- list(
          id = strain_id,
          ref = NA,
          stain = NA
        )
        row_idx <- row_idx + 1
        next
      }
      
      # CASE 2: one or more references → collect all refs/stains into lists
      cell_morphs <- morph$`cell morphology`
      # if only one reference, wrap it in a list
      if (!is.list(cell_morphs[[1]])) {
        cell_morphs <- list(cell_morphs)
      }
      
      all_refs <- lapply(cell_morphs, function(cm) if ("@ref" %in% names(cm)) cm$`@ref` else NA)
      all_stains <- lapply(cell_morphs, function(cm) if ("gram stain" %in% names(cm)) cm$`gram stain` else NA)
      
      rows[[row_idx]] <- list(
        id = strain_id,
        ref = all_refs,    # list of refs
        stain = all_stains # list of stains
      )
      row_idx <- row_idx + 1
    }
  }
  
  # convert list to tibble
  gs <- tibble::tibble(
    id = sapply(rows, `[[`, "id"),
    ref = lapply(rows, `[[`, "ref"),
    stain = lapply(rows, `[[`, "stain")
  )
  
  # Store results per family
  gs_list[[taxon]] <- gs
  message("Finished family: ", taxon, " — ", nrow(gs), " rows collected")
}

saveRDS(gs_list, fname_out)

# Summarize stain data
summarize_strain <- function(stain_list) {
  
  # If NA or empty
  if (is.null(stain_list) || all(is.na(stain_list))) {
    return(list(majority_stain = NA, percentage = NA))
  }
  
  # Flatten in case of nested lists
  stains <- unlist(stain_list)
  stains <- stains[!is.na(stains)]
  
  if (length(stains) == 0) {
    return(list(majority_stain = NA, percentage = NA))
  }
  
  tab <- table(stains)
  majority <- names(tab)[which.max(tab)]
  percent <- max(tab) / sum(tab)
  
  return(list(
    majority_stain = majority,
    percentage = percent
  ))
}

gs_summary_list <- list()

for (taxon in names(gs_list)) {
  
  gs <- gs_list[[taxon]]
  
  majority_vec <- character(nrow(gs))
  percent_vec  <- numeric(nrow(gs))
  
  for (i in seq_len(nrow(gs))) {
    
    result <- summarize_strain(gs$stain[[i]])
    
    majority_vec[i] <- result$majority_stain
    percent_vec[i]  <- result$percentage
  }
  
  gs_summary <- data.frame(
    id = gs$id,
    majority_stain = majority_vec,
    percentage = percent_vec,
    stringsAsFactors = FALSE
  )
  
  gs_summary_list[[taxon]] <- gs_summary
}



family_summary <- data.frame(
  family = character(),
  majority_stain = character(),
  percentage = numeric(),
  stringsAsFactors = FALSE
)

for (taxon in names(gs_summary_list)) {
  
  df <- gs_summary_list[[taxon]]
  
  # Remove NA strains
  stains <- df %>%
    pull(majority_stain) %>%
    na.omit()
  
  if (length(stains) == 0) {
    family_summary <- rbind(
      family_summary,
      data.frame(
        family = taxon,
        majority_stain = NA,
        percentage = NA,
        stringsAsFactors = FALSE
      )
    )
    next
  }
  
  tab <- table(stains)
  majority <- names(tab)[which.max(tab)]
  percent <- max(tab) / sum(tab)
  
  family_summary <- rbind(
    family_summary,
    data.frame(
      family = taxon,
      majority_stain = majority,
      percentage = percent,
      stringsAsFactors = FALSE
    )
  )
}

