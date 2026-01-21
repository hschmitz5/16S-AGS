rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")

fname_ord <- "./figures/ordination-PCoA.png"
fname_ord2 <- "./figures/ordination-PCoA-mu.png"

# load phyloseq object for all sample sizes
ps_full <- readRDS("./data/ps_ASV_full.rds")

# ps_full: all sample groups
ps.ord <- ordinate(ps_full, "PCoA", "wunifrac")

cols <- c("gray", met.brewer(size_pal, n_sizes))

p <- plot_ordination(ps_full, ps.ord, type="samples", color="size.name") +
  scale_color_manual(values = cols) +
  labs(title="PCoA (wunifrac)", color = "Size") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 12)
  ) 

ordination_plot <- p

ggsave(fname_ord, plot = ordination_plot, width = 5, height = 3, dpi = 300)

# Add mu to ps
ps@sam_data$mu <- mu$mu[ match(ps@sam_data$size.name, mu$size.name) ]

ps2.ord <- ordinate(ps, "PCoA", "wunifrac")

p2 <- plot_ordination(ps, ps2.ord, type="samples", color="mu") +
  labs(title="PCoA (wunifrac)", color = "Mu") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 12)
  ) 

ordination_plot_mu <- p2

ggsave(fname_ord2, plot = ordination_plot_mu, width = 5, height = 3, dpi = 300)

# library(vegan)
# ord <- metaMDS(ps, trace = FALSE)
