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
  
  # ----- Subset abundance ------
  # At least one sample has Abundance > rel_ab_cutoff
  high_ab_taxa <- get_rel_ASV(ps) %>%
    filter(OTU %in% all_sig_taxa) %>%
    filter(Abundance > rel_ab_cutoff) %>%
    distinct(OTU) %>%
    pull(OTU)
  
  low_ab_taxa <- setdiff(all_sig_taxa, high_ab_taxa)
  
    
  # --- Write Data to Excel
  if (isTRUE(write2excel)) {
    if (is.null(fname_out)) {
      # Generate a default filename if not provided
      fname_out <- "./data/ANCOM.xlsx"
      message("No fname_out provided. Writing Excel file to default: ", fname_out)
    }
    
    # Make a single data frame with two columns
    taxa_df <- data.frame(
      high_abundance = c(sort(high_ab_taxa), rep(NA, length(low_ab_taxa) - length(high_ab_taxa))),
      low_abundance  = sort(low_ab_taxa)
    )
    # Write to Excel
    write_xlsx(taxa_df, path = fname_out)
  } 
  
  return(list(
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
  
  # ----- Subset abundance ------
  # At least one sample has Abundance > rel_ab_cutoff
  high_ab_taxa <- get_rel_ASV(ps) %>%
    filter(OTU %in% all_sig_taxa) %>%
    filter(Abundance > rel_ab_cutoff) %>%
    distinct(OTU) %>%
    pull(OTU)
  
  low_ab_taxa <- setdiff(all_sig_taxa, high_ab_taxa)
  
  # --- Write Data to Excel
  if (isTRUE(write2excel)) {
    if (is.null(fname_out)) {
      # Generate a default filename if not provided
      fname_out <- "./data/ALDEx.xlsx"
      message("No fname_out provided. Writing Excel file to default: ", fname_out)
    }
    
    # Make a single data frame with two columns
    taxa_df <- data.frame(
      high_abundance = c(sort(high_ab_taxa), rep(NA, length(low_ab_taxa) - length(high_ab_taxa))),
      low_abundance  = sort(low_ab_taxa)
    )
    # Write to Excel
    write_xlsx(taxa_df, path = fname_out)
  }
  
  return(list(
    high_ab = high_ab_taxa,
    low_ab  = low_ab_taxa
  ))
}

get_common_taxa <- function(ancom_fname, aldex_fname, p_threshold, effect_threshold) {
  # names of significant taxa
  all_taxa_ancom <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel = NULL, fname_out) 
  high_ancom <- all_taxa_ancom$high_ab
  
  all_taxa_aldex <- get_aldex_taxa(aldex_fname, ps, p_threshold, effect_threshold, rel_ab_cutoff, write2excel = NULL, fname_out)
  high_aldex <- all_taxa_aldex$high_ab
  
  common_taxa <- intersect(high_ancom, high_aldex) 
  
  return(common_taxa)
}