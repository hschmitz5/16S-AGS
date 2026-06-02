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