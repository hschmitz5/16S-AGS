rm(list = ls())

library(phyloseq)
library(ANCOMBC)

fname_out <- "./data/ancombc2_genus.rds"

ps <- readRDS("./data/ps_genus_subset.rds") 

set.seed(123)
output <- ancombc2(
  data = ps,
  fix_formula = "size.name",    
  group = "size.name",
  struc_zero = TRUE,
  global = TRUE, dunnet = TRUE, 
  pairwise = FALSE, trend = FALSE   
)

saveRDS(output, file = fname_out)

## use if trend = TRUE
# 
# contrast_mats = list(
#   # monotonically increasing
#   matrix(c(1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1),
#          nrow = 4, byrow = TRUE),
#   # monotonically decreasing
#   matrix(c(-1, 0, 0, 0, 1, -1, 0, 0, 0, 1, -1, 0, 0, 0, 1, -1),
#          nrow = 4, byrow = TRUE)
# )
# 
# trend_control = list(contrast = contrast_mats,
#                      node = list(4, 4),
#                      solver = "ECOS",
#                      B = 100)