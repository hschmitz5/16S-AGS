get_avg_genus <- function(ps){
  # Average per Genus across all samples
  genus_avg <- get_rel_genus(ps) %>%
    filter(!is.na(Genus)) %>%
    group_by(Genus) %>%
    summarise(
      mean_ab = mean(Abundance),
      sd_ab = sd(Abundance),
      .groups = "drop") %>%
    filter(mean_ab > rel_ab_cutoff) %>%
    arrange(desc(mean_ab)) 
}

get_OTUs_above_cutoff <- function(ps){
  # Many genera have multiple ASVs in which one ASV is almost zero
  # Few have ASVs with a non-negligible abundance
  OTU_subset <- get_rel_ASV(ps) %>%
    filter(Abundance > rel_ab_cutoff) %>%
    # rename OTUs
    mutate(
      Species_tmp = case_when(
        !is.na(Species) & !startsWith(Species, "midas") ~ Species,
        !is.na(Genus)   & !startsWith(Genus, "midas")   ~ paste0(Genus, "_sp_g"),
        !is.na(Family)  & !startsWith(Family, "midas")  ~ paste0(Family, "_sp_f"),
        !is.na(Order)   & !startsWith(Order, "midas")   ~ paste0(Order, "_sp_o"),
        !is.na(Class)   & !startsWith(Class, "midas")   ~ paste0(Class, "_sp_c")
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
    select(-Species_tmp) %>% # delete column
    arrange(desc(Abundance)) 
}


# Define genera above rel_ab_cutoff
genus_avg <- get_avg_genus(ps)
high_ab_genera <- genus_avg$Genus
rm(genus_avg)

# Define OTUs above rel_ab_cutoff
OTU_subset <- get_OTUs_above_cutoff(ps)
high_ab_OTUs <- OTU_subset$Species_updated 
rm(OTU_subset)