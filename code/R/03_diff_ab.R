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