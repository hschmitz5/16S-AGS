rm(list = ls())

library(phyloseq)
library(ANCOMBC)

fname_out <- "./data/ancombc2_ASV.rds"

ps_ASV <- readRDS("./data/ps_ASV_subset.rds") 

# Remove taxa not seen more than 3 times (reads) in at least 20% of the samples.
# This protects against an OTU with small mean & trivially large C.V.
ps_filt = filter_taxa(ps_ASV, function(x) sum(x > 3) >= (0.2*length(x)), TRUE)

contrast_mats = list(
  # monotonically increasing
  matrix(c(1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1),
         nrow = 4, byrow = TRUE),
  # monotonically decreasing
  matrix(c(-1, 0, 0, 0, 1, -1, 0, 0, 0, 1, -1, 0, 0, 0, 1, -1),
         nrow = 4, byrow = TRUE)
)

set.seed(123)
output <- ancombc2(
  data = ps_filt,
  fix_formula = "size.name",    
  group = "size.name",
  struc_zero = TRUE,
  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE #, #TRUE,  
  # trend_control = list(contrast = contrast_mats,
  #                      node = list(4, 4),
  #                      solver = "ECOS",
  #                      B = 100)
)

saveRDS(output, file = fname_out)
