# csv files acquired from advanced search on Bac Dive
rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")
library(BacDive)
library(purrr)

fname_out <- "./data/bd_family.rds"

# Choose DA Taxa
ancom_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel = FALSE)
DA_taxa <- ancom_taxa$high_ab 

# DA_order <- taxonomy %>%
#   filter(Species_updated %in% DA_taxa) %>%
#   mutate(
#     g_suffix = if_else(
#       str_detect(Species_updated, "_o-[0-9]+$"),
#       str_replace(Species_updated, "_o-[0-9]+$", ""),
#       NA_character_
#     )
#   ) %>%
#   pull(g_suffix) %>%
#   na.omit() %>%
#   unique()

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

bd <- open_bacdive("hannah.schmitz@northwestern.edu", "QP4,@,_nunfQEY6")

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

gs_list <- list()  # list to hold results per DA_family

for (n in seq_along(DA_family)) {
  fname <- paste0("./data/Bacdive/", DA_family[n], "_f.csv")
  
  # check if the file exists
  if (!file.exists(fname)) {
    message("File not found for family: ", DA_family[n], " — skipping")
    next  # skip to the next family
  } else {
    message("Processing family: ", DA_family[n])
  }
  
  # read in data
  df_raw <- read.csv(
    fname, 
    skip = 2, # skip header lines
    stringsAsFactors = FALSE
  )
  
  # filter only the current family
  df <- df_raw %>%
    filter(family %in% DA_family[n]) %>%
    arrange(ID)
  
  removed <- nrow(df_raw) - nrow(df)
  cat("Rows removed:", removed, "\n")
  
  
  # initialize data frame for this family
  gs <- data.frame(
    id    = character(),
    ref   = character(),
    stain = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(df))) {
    strain_id <- as.character(df$ID[i])
    message("Processing: ", strain_id)
    strain <- safe_fetch(strain_id)
    
    if (is.null(strain)) next
    
    morph <- strain$results[[strain_id]]$Morphology
    
    # CASE 1: no morphology → NA row
    if (is.null(morph) || length(morph) == 0) {
      gs[nrow(gs) + 1, ] <- list(
        id    = strain_id,
        ref   = NA,
        stain = NA
      )
      next
    }
    
    # CASE 2: one or more references → one row per reference
    for (j in seq_along(morph$`cell morphology`)) {
      cm <- morph$`cell morphology`[[j]]
      
      gs[nrow(gs) + 1, ] <- list(
        id    = strain_id,
        ref   = if ("@ref" %in% names(cm)) cm$`@ref` else NA,
        stain = if ("gram stain" %in% names(cm)) cm$`gram stain` else NA
      )
    }
  }
  
  # store this family’s result in the list
  gs_list[[DA_family[n]]] <- gs
}




####

gs_list <- list()  # list to hold results per DA_family

for (n in seq_along(DA_family)) {
  fname <- paste0("./data/Bacdive/", DA_family[n], "_f.csv")
  
  # check if the file exists
  if (!file.exists(fname)) {
    message("File not found for family: ", DA_family[n], " — skipping")
    next  # skip to the next family
  } else {
    message("Processing family: ", DA_family[n])
  }
  
  # read in data
  df_raw <- read.csv(
    fname, 
    skip = 2, # skip header lines
    stringsAsFactors = FALSE
  )
  
  # filter only the current family
  df <- df_raw %>%
    filter(family %in% DA_family[n]) %>%
    arrange(ID)
  
  removed <- nrow(df_raw) - nrow(df)
  cat("Rows removed:", removed, "\n")
  

  # initialize data frame for this family
  gs <- data.frame(
    id    = character(),
    ref   = character(),
    stain = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(df))) {
    strain_id <- as.character(df$ID[i])
    strain <- safe_fetch(strain_id, bd)
    Sys.sleep(0.5)
    
    morph <- strain$results[[strain_id]]$Morphology
    
    # CASE 1: no morphology → NA row
    if (is.null(morph) || length(morph) == 0) {
      gs[nrow(gs) + 1, ] <- list(
        id    = strain_id,
        ref   = NA,
        stain = NA
      )
      next
    }
    
    # CASE 2: one or more references → one row per reference
    for (j in seq_along(morph$`cell morphology`)) {
      cm <- morph$`cell morphology`[[j]]
      
      gs[nrow(gs) + 1, ] <- list(
        id    = strain_id,
        ref   = if ("@ref" %in% names(cm)) cm$`@ref` else NA,
        stain = if ("gram stain" %in% names(cm)) cm$`gram stain` else NA
      )
    }
  }
  
  # store this family’s result in the list
  gs_list[[DA_family[n]]] <- gs
}


fname_out_final <- "./data/bd_family.rds"
saveRDS(gs_list, fname_out_final)

# initialize a summary table
gram_summary <- data.frame(
  stain = c("Gram-positive", "Gram-negative"),
  stringsAsFactors = FALSE
)

# loop over families in gs_list
for (fam in names(gs_list)) {
  
  df_family <- gs_list[[fam]]
  
  # count positives and negatives (ignore NA)
  counts <- table(df_family$stain, useNA = "no")
  
  # add counts to summary table
  gram_summary[[fam]] <- c(
    counts["Gram-positive"] %||% 0,  # %||% is a convenient way to replace NULL with 0
    counts["Gram-negative"] %||% 0
  )
}

