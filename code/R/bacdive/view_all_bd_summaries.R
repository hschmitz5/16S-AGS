parent_fname <- "./data/bacdive/"

fname_g <- "bd_genera_summary.rds"
fname_f <- "bd_family_summary.rds"
fname_o <- "bd_order_summary.rds"

g <- readRDS(paste0(parent_fname,fname_g))
f <- readRDS(paste0(parent_fname,fname_f))
o <- readRDS(paste0(parent_fname,fname_o))
