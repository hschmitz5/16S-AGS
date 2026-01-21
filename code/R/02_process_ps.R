get_metadata <- function(ps){
  metadata <- as.data.frame(as.matrix(ps@sam_data)) %>%
  rownames_to_column("Sample") %>%
  mutate(size.mm = factor(size.mm, levels = size$ranges),
         size.name = factor(size.name, levels = size$name))
}

get_taxonomy <- function(ps){
  taxonomy <- as.data.frame(as.matrix(ps@tax_table)) %>%
  rownames_to_column("OTU") 
}

# get_rel_genus <- function(ps) {
#   metadata <- get_metadata(ps)
#   taxonomy <- get_taxonomy(ps)
#   
#   # define relative abundance
#   ps_rel_ASV <- phyloseq::transform_sample_counts(ps, function(x) x*100/sum(x))
#   
#   table_rel <- as.data.frame(as.matrix(ps_rel_ASV@otu_table)) %>%
#     rownames_to_column(var = "OTU")
#   
#   table_rel_long <- table_rel %>%
#     pivot_longer(cols = !OTU, names_to = "Sample", values_to = "tmp_abund") %>% 
#     left_join(taxonomy, join_by(OTU)) %>%
#     left_join(metadata, join_by(Sample)) %>%
#     group_by(Sample, size.name, size.mm, Genus) %>%
#     # Adds multiple ASVs in the same Genus
#     summarise(Abundance = sum(tmp_abund), .groups = "drop")
#   
#   return(table_rel_long)
# }

get_rel_ASV <- function(ps) {
  # define relative abundance
  ps_rel_ASV <- phyloseq::transform_sample_counts(ps, function(x) x*100/sum(x))
  
  phyloseq::psmelt(ps_rel_ASV)
}