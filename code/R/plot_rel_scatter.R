rm(list = ls())
library(ggtext)
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")

fname <- "./figures/abund_size_scatter.png"

shapes <- c(1, 17, 0, 16, 2, 15, 1) 

# names of significant taxa
# write2excel = FALSE
# sig_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel, fname_excel)
# plot_taxa <- sig_taxa$high_ab[1:n_display]
# filter(OTU %in% plot_taxa) 

plot_genus <- "Ca_Competibacter"

# 1) Taxa passing cutoff in any group
plot_taxa <- get_rel_ASV(ps) %>%
  filter(startsWith(OTU, plot_genus),
         Abundance > rel_ab_cutoff) %>%
  distinct(OTU) %>%
  pull(OTU)

# 2) Summarise once
ASV_size <- get_rel_ASV(ps) %>%
  filter(OTU %in% plot_taxa) %>%
  group_by(size.mm, size.name, OTU) %>%
  summarise(
    mean_ab = mean(Abundance),
    std_ab  = sd(Abundance),
    .groups = "drop"
  )

# 3) Derive ordering from that summary
taxa_order <- ASV_size %>%
  filter(size.name == "S") %>%
  arrange(desc(mean_ab)) %>%
  pull(OTU)

# 4) Apply factor levels
ASV_size <- ASV_size %>%
  mutate(OTU = factor(OTU, levels = taxa_order))

n_display <- length(plot_taxa)

# italicize some labels
otu_levels <- levels(ASV_size$OTU)

italic_rows <- 
  !grepl("_(f|o|c|p)(?:_|$|-)", otu_levels) 

otu_labels <- ifelse(
  italic_rows,
  paste0("<i>", otu_levels, "</i>"),
  otu_levels
)

names(otu_labels) <- otu_levels

p <- ggplot(ASV_size, aes(x = size.name, y = mean_ab, colour = OTU, shape = OTU)) +
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

# Save plot
ggsave(fname, plot = p, width = 6, height = 4, dpi = 300)
