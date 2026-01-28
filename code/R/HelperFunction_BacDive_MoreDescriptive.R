# ============================================================
# BacDive Oxygen Tolerance Extraction & Genus-Level Summary
#
# This script:
#   1. Queries BacDive for each genus in `lung_genera`
#   2. Extracts oxygen tolerance information for each strain
#   3. Counts oxygen tolerance categories per genus
#   4. Produces a genus-level summary table
# ============================================================


# ----------------------------
# Connect to BacDive
# ----------------------------
# Opens a session to the BacDive database using user credentials.
# The returned object (`bd`) is required for all subsequent queries.
bd <- open_bacdive("HannahSchmitz2026@u.northwestern.edu", "UgQ2ac7fULX4SpC/")
# bd <- open_bacdive("your.email@domain.com", "your_password")

# ----------------------------
# Helper function: extract oxygen tolerance from a single BacDive record
# ----------------------------
# This function:
#   - Safely navigates the nested list
#   - Handles both simple and multi-entry oxygen annotations
#   - Returns a single semicolon-delimited character string
#   - Returns NA if no oxygen information is available
extract_oxygen <- function(out) {
  
  # Pull the "Physiology and metabolism" section if it exists
  phys <- out[["Physiology and metabolism"]]
  
  # Extract oxygen tolerance from physiology section if present
  oxy  <- if (!is.null(phys)) phys[["oxygen tolerance"]] else NULL
  
  # If oxygen tolerance is completely missing, return NA
  if (is.null(oxy)) return(NA_character_)
  
  # Case 1: oxygen tolerance is stored as a list of sub-lists
  # (multiple structured entries)
  if (is.list(oxy) && all(sapply(oxy, is.list))) {
    
    # Extract oxygen tolerance from each sub-entry
    vals <- vapply(
      oxy,
      function(x) {
        if (!is.null(x[["oxygen tolerance"]])) {
          x[["oxygen tolerance"]]
        } else {
          NA_character_
        }
      },
      character(1)
    )
    
  } else {
    
    # Case 2: oxygen tolerance is stored in a simpler structure
    vals <- oxy[["oxygen tolerance"]]
  }
  
  # Remove missing values
  vals <- vals[!is.na(vals)]
  
  # If nothing remains after filtering, return NA
  if (length(vals) == 0) return(NA_character_)
  
  # Collapse unique oxygen tolerance values into a single string
  # (e.g., "aerobic; facultative anaerobic")
  paste(unique(vals), collapse = "; ")
}


# ----------------------------
# Initialize master summary table
# ----------------------------
# This tibble will store genus-level oxygen tolerance summaries
# and will be appended to during the loop below.
master_summary <- tibble()


# ----------------------------
# Loop over all genera of interest
# ----------------------------
for (g in lung_genera) {
  
  # Print progress to console
  message("=== Processing genus: ", g, " ===")
  
  # Query BacDive for all records matching this genus
  # Wrapped in try() to prevent the script from failing
  # if BacDive returns an error for a given genus
  recs <- try(
    retrieve(object = bd, query = g, search = "taxon"),
    silent = TRUE
  )
  
  # Skip genus if query fails or returns no records
  if (inherits(recs, "try-error") || length(recs) == 0) {
    message("No records found for ", g, " — skipping")
    next
  }
  
  # ----------------------------
  # Build strain-level table
  # ----------------------------
  # Each BacDive record corresponds to a strain/species entry.
  # We extract taxonomy, BacDive ID, and oxygen tolerance.
  strain_tbl <- purrr::map_dfr(recs, function(out) {
    
    # Pull taxonomic classification
    nm <- out[["Name and taxonomic classification"]]
    
    # Return one row per BacDive record
    tibble(
      Genus     = nm$genus,
      Species   = nm$species,
      BacDiveID = out[["General"]][["BacDive-ID"]],
      Oxygen    = extract_oxygen(out)
    )
  })
  
  
  # ----------------------------
  # Summarize oxygen tolerance at the genus level
  # ----------------------------
  summary_tbl <- strain_tbl %>%
    
    # Replace missing oxygen tolerance values with explicit label
    mutate(Oxygen = ifelse(is.na(Oxygen), "Unknown", Oxygen)) %>%
    
    # Split multiple oxygen categories into separate rows
    # (e.g., "aerobic; microaerophilic" → two rows)
    separate_rows(Oxygen, sep = ";\\s*") %>%
    
    # Count number of strains per genus and oxygen category
    count(Genus, Oxygen, name = "StrainCount") %>%
    
    # Convert from long to wide format:
    # one column per oxygen tolerance category
    pivot_wider(
      names_from  = Oxygen,
      values_from = StrainCount,
      values_fill = 0
    ) %>%
    
    # Calculate total counts across all oxygen categories
    # NOTE: Strains with multiple oxygen annotations contribute
    # to multiple categories.
    mutate(
      TotalStrains = rowSums(across(-Genus))
    )
  
  
  # ----------------------------
  # Append to master summary
  # ----------------------------
  master_summary <- bind_rows(master_summary, summary_tbl)
}


# At the end of the script:
# - `master_summary` contains genus-level oxygen tolerance counts
# - Each row corresponds to a genus
# - Columns represent oxygen tolerance categories + TotalStrains
