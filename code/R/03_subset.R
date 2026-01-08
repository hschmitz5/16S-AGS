get_ancom_taxa <- function(fname_in, ps, p_threshold, rel_ab_cutoff, write2excel, fname_out) {
  
  output <- readRDS(fname_in)
  
  taxonomy <- get_taxonomy(ps)
  
  # Subset based on significance
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
  
  # Rename ASVs based on taxonomy
  sig_taxa_tbl <- taxonomy %>%
    filter(OTU %in% all_sig_taxa) %>%
    mutate(
      Species_tmp = case_when(
        !is.na(Species) & !startsWith(Species, "midas") ~ Species,
        !is.na(Genus)   & !startsWith(Genus, "midas")   ~ paste0(Genus, "_sp_g"),
        !is.na(Family)  & !startsWith(Family, "midas")  ~ paste0(Family, "_sp_f"),
        !is.na(Order)   & !startsWith(Order, "midas")   ~ paste0(Order, "_sp_o"),
        !is.na(Class)   & !startsWith(Class, "midas")   ~ paste0(Class, "_sp_c"),
        !is.na(Phylum)  & !startsWith(Phylum, "midas")  ~ paste0(Phylum, "_sp_p"),
        .default = Species 
      )
    ) %>%
    group_by(Species_tmp) %>%
    mutate(
      Species_updated = if (n() == 1) {
        Species_tmp
      } else {
        paste0(Species_tmp, "-", row_number())  
      }
    ) %>%
    ungroup() %>%
    select(-Species_tmp)
  
  # ----- Subset abundance ------
  # At least one sample has Abundance > rel_ab_cutoff
  high_ab_taxa <- get_rel_ASV(ps) %>%
    filter(OTU %in% all_sig_taxa) %>%
    left_join(sig_taxa_tbl, join_by(OTU)) %>%
    filter(Abundance > rel_ab_cutoff) %>%
    distinct(Species_updated) %>%
    pull(Species_updated)
  
  low_ab_taxa <- setdiff(sig_taxa_tbl$Species_updated, high_ab_taxa)
  
    
  # --- Write Data to Excel
  if (isTRUE(write2excel)) {
    # Make a single data frame with two columns
    taxa_df <- data.frame(
      high_abundance = c(sort(high_ab_taxa), rep(NA, length(low_ab_taxa) - length(high_ab_taxa))),
      low_abundance  = sort(low_ab_taxa)
    )
    # Write to Excel
    write_xlsx(taxa_df, path = fname_out)
  } 
  return(list(
    name_tbl = sig_taxa_tbl,
    high_ab = high_ab_taxa,
    low_ab  = low_ab_taxa
  ))
}

get_aldex_taxa <- function(fname_in, ps, p_threshold, effect_threshold, rel_ab_cutoff, write2excel, fname_out) {
  
  output <- readRDS(fname_in) 
  
  taxonomy <- get_taxonomy(ps)
  
  all_sig_taxa <- output %>%
    mutate(
      sig_taxa = map(res, ~ .x %>%
                       rownames_to_column("OTU") %>%
                       filter(wi.eBH < p_threshold & abs(effect) > effect_threshold) %>%
                       pull(OTU)
      )
    ) %>%
    unnest(sig_taxa) %>%
    distinct(OTU = sig_taxa) %>%
    left_join(taxonomy, join_by(OTU)) %>%
    filter(!is.na(Phylum)) %>%
    pull(OTU)
  
  # Rename ASVs based on taxonomy
  sig_taxa_tbl <- taxonomy %>%
    filter(OTU %in% all_sig_taxa) %>%
    mutate(
      Species_tmp = case_when(
        !is.na(Species) & !startsWith(Species, "midas") ~ Species,
        !is.na(Genus)   & !startsWith(Genus, "midas")   ~ paste0(Genus, "_sp_g"),
        !is.na(Family)  & !startsWith(Family, "midas")  ~ paste0(Family, "_sp_f"),
        !is.na(Order)   & !startsWith(Order, "midas")   ~ paste0(Order, "_sp_o"),
        !is.na(Class)   & !startsWith(Class, "midas")   ~ paste0(Class, "_sp_c"),
        !is.na(Phylum)  & !startsWith(Phylum, "midas")  ~ paste0(Phylum, "_sp_p"),
        .default = Species 
      )
    ) %>%
    group_by(Species_tmp) %>%
    mutate(
      Species_updated = if (n() == 1) {
        Species_tmp
      } else {
        paste0(Species_tmp, "-", row_number())  
      }
    ) %>%
    ungroup() %>%
    select(-Species_tmp)
  
  # ----- Subset abundance ------
  # At least one sample has Abundance > rel_ab_cutoff
  high_ab_taxa <- get_rel_ASV(ps) %>%
    filter(OTU %in% all_sig_taxa) %>%
    left_join(sig_taxa_tbl, join_by(OTU)) %>%
    filter(Abundance > rel_ab_cutoff) %>%
    distinct(Species_updated) %>%
    pull(Species_updated)
  
  low_ab_taxa <- setdiff(sig_taxa_tbl$Species_updated, high_ab_taxa)
  
  
  # --- Write Data to Excel
  if (isTRUE(write2excel)) {
    # Make a single data frame with two columns
    taxa_df <- data.frame(
      high_abundance = c(sort(high_ab_taxa), rep(NA, length(low_ab_taxa) - length(high_ab_taxa))),
      low_abundance  = sort(low_ab_taxa)
    )
    # Write to Excel
    write_xlsx(taxa_df, path = fname_out)
  }
  return(list(
    name_tbl = sig_taxa_tbl,
    high_ab = high_ab_taxa,
    low_ab  = low_ab_taxa
  ))
}

get_common_taxa <- function(ancom_fname, aldex_fname, p_threshold, effect_threshold, high_ab_OTUs) {
  # names of significant taxa
  all_taxa_ancom <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel, fname_out) 
  high_ancom <- all_taxa_ancom$high_ab
  
  all_taxa_aldex <- get_aldex_taxa(aldex_fname, ps, p_threshold, effect_threshold, rel_ab_cutoff, write2excel, fname_out)
  high_aldex <- all_taxa_aldex$high_ab
  
  common_taxa <- intersect(high_ancom, high_aldex) 
  
  return(common_taxa)
}