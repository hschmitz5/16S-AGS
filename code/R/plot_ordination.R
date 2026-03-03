rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")

fname_ord <- "./figures/ordination-PCoA.png"
fname_ord2 <- "./figures/ordination-PCoA-mu.png"

# load phyloseq object for all sample sizes
ps_full <- readRDS("./data/ps_ASV_full.rds")

# ps_full: all sample groups
ps.ord <- ordinate(ps_full, "PCoA", "wunifrac")

# colors
cols <- c("gray", met.brewer(size_pal, n_sizes))

# symbol
# 16 = filled circle, 17 = triangle, 15 = square, 18 = diamond, etc.
shapes <- c(16, 17, 15, 18, 3, 7)

p <- plot_ordination(ps_full, ps.ord, type="samples", color="size.name", shape = "size.name") +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  labs(color = "Size", shape = "Size") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 12)
  ) 

ordination_plot <- p

ggsave(fname_ord, plot = ordination_plot, width = 5, height = 3, dpi = 600)