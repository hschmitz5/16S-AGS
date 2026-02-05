#rm(list = ls())
library(phyloseq)

# Relative Abundance Cutoff (%) used to subset high abundance taxa
rel_ab_cutoff <- 0.5
# p-value used for filtering taxa (alpha)
p_threshold   <- 0.05
effect_threshold <- 1 # used for ALDEx DA

# Color Palettes (MetBrewer)
size_pal <- "Java"
taxa_pal <- "Hiroshige" #"Hokusai2"

# File Names
ps_fname    <- "./data/ps_ASV_subset.rds"
ancom_fname <- "./data/ancombc2_ASV.rds"
aldex_fname <- "./data/aldex_ASV.rds"
mech_fname  <- "./data/EPS_moduli.xlsx"

# absolute counts
ps <- readRDS(ps_fname)

# Change OTU to Species_updated
taxonomy <- as.data.frame(as.matrix(ps@tax_table))
rownames(ps@otu_table) <- taxonomy$Species_updated
rownames(ps@tax_table) <- taxonomy$Species_updated
ps@phy_tree$tip.label <- taxonomy$Species_updated

# define dimensions of sample grouping
n_replicates  <- 3
n_sizes <- length(levels(ps@sam_data$size.mm))

# For ordering
sam_name <- c("20A", "20B", "20C", "14A", "14B", "14C", "10A", "10B", "10C",
              "7A", "7B", "7C", "5A", "5B", "5C")
# For correlation
size_midpoint <- c(1.125, 1.7, 2.4, 3.4, 4.5)

# define sample names
size <- data.frame(
  ranges = levels(ps@sam_data$size.mm),
  name = levels(ps@sam_data$size.name)
)

# Other Data
eps <- read_excel(mech_fname, range = cell_cols("A:E"))
eps$size.name = factor(eps$size.name, levels = size$name)
mu <- read_excel(mech_fname, range = cell_cols("G:H"))
mu$size.name = factor(mu$size.name, levels = size$name)

modulus <- read_excel(mech_fname, range = cell_cols("K:N")) |>
  as.data.frame()

rownames(modulus) <- mu$size.name
