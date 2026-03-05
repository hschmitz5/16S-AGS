rm(list = ls())
library(ggtext)
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")

fname_p4 <- "./figures/moduli_corr_taxa.png"

bs <- 12 # base size
ar <- 0.7 # aspect ratio

corr_taxa <- c("Ardenticatenales_o-13", "Ca_Competibacter_g-75", "Lysobacter_g-3", "Terrimonas_g-8")

# 4 taxa it corr_taxa
n_display <- 4

shapes <- c(16, 17, 15, 18)

# Get relative abundance corresponding to corr_taxa
ASV_size <- get_rel_ASV(ps) %>%
  filter(OTU %in% corr_taxa) %>%
  mutate(OTU = factor(OTU, levels = corr_taxa)) %>%
  group_by(size.mm, size.name, OTU) %>%
  summarise(
    mean_ab = mean(Abundance),
    std_ab = sd(Abundance),
    .groups = "drop"
  ) 

taxa_order <- ASV_size %>%
  filter(size.name == "S") %>%
  arrange(desc(mean_ab)) %>%
  pull(OTU)

ASV_size <- ASV_size %>%
  mutate(OTU = factor(OTU, levels = taxa_order))
  
# italicize some labels
otu_levels <- levels(taxa_order)

italic_rows <- grepl("_(g|s)(?:_|$|-)", otu_levels)

otu_labels <- ifelse(
  italic_rows,
  paste0("<i>", otu_levels, "</i>"),
  otu_levels
)

#### Size vs Modulus

p1 <- ggplot(modulus, aes(x = size.name, y = G1.avg)) +
  geom_point(size = 2) + 
  geom_errorbar(
    aes(ymin = G1.avg - G1.sd, ymax = G1.avg + G1.sd),
    width = 0.1 
  ) +
  labs(
    x = "Size",
    y = "Storage Modulus [Pa]"
  ) +
  theme_minimal(base_size = bs) +
  theme(aspect.ratio = ar)
            
#### Size vs Taxa that correlate to modulus by size
       
p2 <- ggplot(ASV_size, aes(x = size.name, y = mean_ab, colour = OTU, shape = OTU)) +
  geom_point(position = position_dodge(width = 0.2), size = 2) + 
  geom_errorbar(
    aes(ymin = mean_ab - std_ab, ymax = mean_ab + std_ab),
    width = 0.2,
    position = position_dodge(width = 0.2)
  ) +
  scale_colour_manual(
    values = met.brewer(taxa_pal, n_display),
    labels = otu_labels,
    name = "ASV"
  ) +
  scale_shape_manual(
    values = shapes,
    labels = otu_labels,
    name = "ASV"
  ) +
  labs(
    x = "Size",
    y = "Relative Abundance [%]"
  ) +
  theme_minimal(base_size = bs) +
  theme(
    legend.position = "right",
    legend.justification = "left",
    aspect.ratio = ar,
    legend.text = element_markdown()
  ) 

# PN/PS
source("./code/R/plot_pnps.R")

p3 <- ggplot(data = eps_long, aes(x = size.name, y = pnps.avg, color = extract.type)) +
  geom_point(position = position_dodge(width = 0.2), size = 2) +
  geom_errorbar(
    aes(ymin = pnps.avg - pnps.sd, ymax = pnps.avg + pnps.sd),
    width = 0.2,
    position = position_dodge(width = 0.2)
  ) +
  scale_color_manual(
    values = met.brewer(taxa_pal, 2),
    labels = c(
      LB = "LB (p.adj = 0.517)",
      TB = "TB (p.adj = 0.389)"
    )
  ) +
  labs(
    x = "Size",
    y = "Mean PN/PS",
    color = "Extract Type"
  ) +
  theme_minimal(base_size = bs) +
  theme(
    legend.position = "right",
    legend.justification = "left",
    aspect.ratio = ar
  ) 

# rel scatter
source("./code/R/plot_rel_scatter.R")

p4 <- ggplot(ASV_size, aes(x = size.name, y = mean_ab, colour = OTU, shape = OTU)) +
  geom_point(position = position_dodge(width = 0.2), size = 2) + 
  geom_errorbar(
    aes(ymin = mean_ab - std_ab, ymax = mean_ab + std_ab),
    width = 0.2,
    position = position_dodge(width = 0.2)
  ) +
  scale_colour_manual(
    values = met.brewer(taxa_pal, n_display),
    labels = otu_labels,
    name = "ASV"
  ) +
  scale_shape_manual(
    values = shapes,
    labels = otu_labels,
    name = "ASV"
  ) +
  labs(
    x = "Size",
    y = "Relative Abundance [%]"
  ) +
  theme_minimal(base_size = bs) +
  theme(
    legend.position = "right",
    legend.justification = "left",
    aspect.ratio = ar,
    legend.text = element_markdown()
  ) 

combined_plot <- p1 / p2 / p3 / p4

ggsave(fname_p4, plot = combined_plot, width = 6, height = 10, dpi = 600)

