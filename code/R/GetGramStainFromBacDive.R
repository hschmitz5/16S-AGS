###############################################################################
# Lung Microbiome Genera Oxygen Tolerance Analysis via BacDive
#
# Description:
#   This script processes lung microbiome data from a phyloseq object. Genera
#   are collapsed to the genus level, contaminants are removed, and BacDive is
#   queried for oxygen tolerance information. The final output is a summary
#   table of oxygen tolerance profiles for lung-associated genera.
###############################################################################

# ============================
# 0. Load required packages
# ============================
rm(list = ls())
suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(purrr)
  library(BacDive)
})

# ============================
# 1. Collapse to Genus level and clean phyloseq object
# ============================
# Input: phyloseq object `asv_table_cp_20240625$Amplicon`
# Goal: remove contaminants, collapse to genus, and filter unwanted taxa

ps_fname    <- "./data/ps_ASV_subset.rds"

# absolute counts
ps <- readRDS(ps_fname)

df_genus <- ps

# Extract abundance matrix
otu <- as.data.frame(otu_table(df_genus))
if (taxa_are_rows(df_genus)) otu <- t(otu)
otu <- as.data.frame(otu)

# define taxonomy
tax <- as.data.frame(tax_table(df_genus))

# Find OTUs with non-NA species
keep_otus <- !is.na(tax$Species) & !grepl("^midas", tax$Species, ignore.case = TRUE)
# keep_otus <- !is.na(tax$Species)

# Delete OTU for NA values
otu <- otu[, keep_otus, drop = FALSE]
tax <- tax[keep_otus, , drop = FALSE]

# Add genus names
colnames(otu) <- tax$Species # or Genus

all_species <- gsub("_", " ", unique(tax$Species))

# ============================
# 2. Metadata alignment
# ============================
# Align sample metadata with abundance table
source("./code/R/02_process_ps.R")
meta <- get_metadata(ps)
otu$Sample <- rownames(otu)
#otu_meta <- left_join(otu, meta, by = "Sample") $dislikes duplicate species

# ============================
# 3. Build genus list
# ============================
# Summarize genera by group (segment) and compute abundance metrics

# groups <- unique(otu_meta$segment)
# 
# summarize_group <- function(df, group_name) {
#   df_group <- df %>% filter(segment == group_name)
#   df_counts <- df_group %>%
#     select(-SampleID, -segment) %>%
#     mutate(across(everything(), as.numeric)) %>%
#     as.matrix()
#   total_counts <- colSums(df_counts, na.rm = TRUE)
#   avg_counts   <- colMeans(df_counts, na.rm = TRUE)
#   rel_abund <- sweep(df_counts, 1, rowSums(df_counts), "/")
#   avg_rel_abund <- colMeans(rel_abund, na.rm = TRUE)
#   tibble(
#     Genus = colnames(df_counts),
#     RawTotalCounts = as.numeric(total_counts),
#     AvgCountsPerSample = as.numeric(avg_counts),
#     AvgRelativeAbundance = as.numeric(avg_rel_abund)
#   )
# }
# 
# genus_lists <- lapply(groups, function(g) {
#   summarize_group(otu_meta, g) %>%
#     arrange(desc(AvgRelativeAbundance)) %>%
#     pull(Genus)
# })
# names(genus_lists) <- groups
# 
# # Union of all genera across lung groups
# lung_genera <- unique(unlist(genus_lists))
# 
# message("=== Lung genera list constructed ===")
# message("Total lung genera: ", length(lung_genera))

# ============================
# 4. BacDive oxygen tolerance pipeline
# ============================
# Query BacDive for oxygen tolerance information for each genus

bd <- open_bacdive("hannah.schmitz@northwestern.edu", "QP4,@,_nunfQEY6")
#bd <- open_bacdive("your.email@domain.com", "your_password")

extract_gram_stain <- function(out) {
  mor <- out[["Morphology"]][["cell morphology"]][["gram stain"]]
  print(mor)
  gs  <- if (!is.null(mor)) mor[["Gram stain"]] else NULL
  #phys <- out[["Physiology and metabolism"]]
  #oxy  <- if (!is.null(phys)) phys[["oxygen tolerance"]] else NULL
  if (is.null(gs)) return(NA_character_)
  if (is.list(gs) && all(sapply(gs, is.list))) {
    vals <- vapply(
      gs,
      function(x) if (!is.null(x[["Gram stain"]])) x[["Gram stain"]] else NA_character_,
      character(1)
    )
  } else {
    vals <- gs[["Gram stain"]]
  }
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(NA_character_)
  paste(unique(vals), collapse = "; ")
}

master_summary <- tibble()

#for (g in lung_genera) {
#for (g in all_species) {
g <- all_species[1]
  message("=== Processing genus: ", g, " ===")
  
  recs <- try(retrieve(object = bd, query = g, search = "taxon"), silent = TRUE)
  if (inherits(recs, "try-error") || length(recs) == 0) {
    message("No records found for ", g, " — skipping")
    next
  }
  
  strain_tbl <- purrr::map_dfr(recs, function(out) {
    nm <- out[["Name and taxonomic classification"]]
    tibble(
      Genus      = nm$genus,
      Species    = nm$species,
      BacDiveID  = out[["General"]][["BacDive-ID"]],
      Gram_Stain = extract_gram_stain(out)  # !!!
    )
  })
  
  # summary_tbl <- strain_tbl %>%
  #   mutate(Oxygen = ifelse(is.na(Oxygen), "Unknown", Oxygen)) %>%
  #   separate_rows(Oxygen, sep = ";\\s*") %>%
  #   count(Genus, Oxygen, name = "StrainCount") %>%
  #   pivot_wider(names_from = Oxygen, values_from = StrainCount, values_fill = 0) %>%
  #   mutate(TotalStrains = rowSums(across(-Genus)))
  # 
  # master_summary <- bind_rows(master_summary, summary_tbl)
#}

# ============================
# 5. Final result
# ============================
# Save lung-only summary
saveRDS(master_summary, "master_summary_lung.rds")

message("=== Lung-only pipeline complete ===")
print(master_summary)

