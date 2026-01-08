get_bc_abund <- function(fname) {
  
  output <- readRDS(fname)
  
  metadata <- get_metadata(ps)
  
  bc_long <- output$bias_correct_log_table %>%
    rownames_to_column("Genus") %>%
    filter(!is.na(Genus)) %>%
    pivot_longer(-Genus, names_to = "Sample", values_to = "bc_abund") %>%
    left_join(metadata, by = "Sample") %>%
    dplyr::select(Sample, size.mm, size.name, Genus, bc_abund) 
  
  return(bc_long)
}
