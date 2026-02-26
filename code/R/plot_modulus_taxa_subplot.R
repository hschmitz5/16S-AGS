rm(list = ls())
library(ggtext)
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")

fname <- "./figures/moduli_corr_taxa.png"

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

italic_rows <- 
  !grepl("_(f|o|c|p)(?:_|$|-)", otu_levels) 

otu_labels <- ifelse(
  italic_rows,
  paste0("<i>", otu_levels, "</i>"),
  otu_levels
)

#### Size vs Modulus

modulus$f_label <- factor(modulus$f, 
                          levels = 0.1, 
                          labels = "0.1 rad/s")

p1 <- ggplot(modulus, aes(x = size.name, y = G1.avg, colour = f_label)) +
  geom_point(size = 2) + 
  geom_errorbar(
    aes(ymin = G1.avg - G1.sd, ymax = G1.avg + G1.sd),
    width = 0.1 
  ) +
  scale_colour_manual(
    values = met.brewer(taxa_pal, n_display),
    name = "Frequency"
  ) +
  labs(
    x = "Size",
    y = "Storage Modulus [Pa]"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "right",
    legend.justification = "left"
  )
            
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
  theme_minimal(base_size = 10) +
  theme(
    legend.text = element_markdown()
  )

combined_plot <- p1 / p2

ggsave(fname, plot = combined_plot, width = 6, height = 6, dpi = 300)
