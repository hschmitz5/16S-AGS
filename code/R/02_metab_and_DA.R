library(readxl)

# ANCOM-BC2 differential abundance data
ancom_fname <- "./data/ancombc2_ASV.rds"
# Metabolism input file
metab_fname <- "./data/metabolism_midas.xlsx"

get_metabolism <- function(df, metab_fname) {
  m <- read_excel(metab_fname, sheet = "input") # tibble
  
  metab <- df %>%
    dplyr::select(Genus, OTU) %>%
    distinct() %>%
    left_join(., m, by = "Genus") %>%
    dplyr::select(-Genus) %>%
    column_to_rownames("OTU")
}

get_ancom_taxa <- function(fname_in, ps, p_threshold) {
  # load differential abundance data
  output <- readRDS(fname_in)
  
  taxonomy <- get_taxonomy(ps)
  
  # Subset based on significance and passing the sensitivity analysis
  all_sig_taxa <- output$res %>%
    rename(OTU = taxon) %>%
    left_join(taxonomy, join_by(OTU)) %>%
    filter(!is.na(Phylum)) %>%
    # Combines diff_size* and passed_ss* together
    pivot_longer(
      cols = matches("q_size\\.name|passed_ss_size\\.name"),
      names_to = c(".value","size"),
      names_pattern = "(q|passed_ss)_size\\.name(.*)"
    ) %>%
    filter(q < p_threshold & passed_ss == TRUE) %>%   
    distinct(OTU) %>%
    pull(OTU)
}

get_rel_agglom <- function(ps, ancom_fname, rel_ab_cutoff, p_threshold) { 
  # define taxa in which at least one sample has abundance > rel_ab_cutoff
  high_ab_taxa <- get_rel_ASV(ps) %>%
    filter(Abundance > rel_ab_cutoff) %>%
    distinct(OTU) %>%
    pull(OTU)
  
  # define differentially abundant (DA) taxa
  DA_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold)
  
  # define relative abundance
  rel_wide <- get_rel_ASV(ps) %>%
    filter(OTU %in% high_ab_taxa) %>%
    dplyr::select(Genus, OTU, Sample, Abundance) %>%  
    pivot_wider(
      names_from = Sample,
      values_from = Abundance
    ) %>%
    dplyr::select(Genus, OTU, all_of(sam_name)) %>%
    mutate(
      DA = ifelse(OTU %in% DA_taxa, "T", "F")
    )
  
  # agglomerate significant and insignificant taxa separately
  data_mat <- rel_wide %>%
    # define name_prefix - just the core OTU name without "-#" at the end
    mutate(
      name_prefix = if_else(
        grepl("[0-9]$", OTU),       # last character is a number
        sub("-[^-]+$", "", OTU),    # remove last dash and everything after
        OTU                         # last character is a letter → keep full OTU
      )
    ) %>%
    # sum relative abundance of groups with same name_prefix and DA
    group_by(DA, Genus, name_prefix) %>%
    summarise(
      n = n(), # number of rows in group
      OTU_orig = first(OTU),  # only used when n = 1
      across(all_of(sam_name), \(x) sum(x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    # create new name
    mutate(
      OTU = if_else(
        n == 1,
        OTU_orig,
        paste0(name_prefix, "_", DA)
      )
    ) %>%
    select(Genus, OTU, DA, all_of(sam_name))
}

# Only the difference between bias corrected data is meaningful.
# It is not a substitute for relative abundance.
get_bc_abund <- function(fname) {
  
  output <- readRDS(fname)
  
  metadata <- get_metadata(ps)
  
  bc_long <- output$bias_correct_log_table %>%
    rownames_to_column("OTU") %>%
    filter(!is.na(OTU)) %>%
    pivot_longer(-OTU, names_to = "Sample", values_to = "bc_abund") %>%
    left_join(metadata, by = "Sample") %>%
    dplyr::select(Sample, size.mm, size.name, OTU, bc_abund) 
  
  return(bc_long)
}