get_bc_abund <- function(fname) {
  
  output <- readRDS(fname)
  
  metadata <- get_metadata(ps, sam_name)
  
  bc_long <- output$bias_correct_log_table %>%
    rownames_to_column("OTU") %>%
    filter(!is.na(OTU)) %>%
    pivot_longer(-OTU, names_to = "Sample", values_to = "bc_abund") %>%
    left_join(metadata, by = "Sample") %>%
    dplyr::select(Sample, size.mm, size.name, OTU, bc_abund) 
  
  return(bc_long)
}
