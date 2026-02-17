library(ggtext)
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")

# display top n most abundant genera
n_display <- 5

shapes <- c(16, 17, 15, 18, 3)

fname <- "./figures/abund_size_scatter.png"

# names of significant taxa
write2excel = FALSE
sig_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel, fname_excel)

ASV_size <- get_rel_ASV(ps) %>%
  dplyr::select(Sample, size.mm, size.name, OTU, Abundance) %>%
  group_by(Sample, size.mm, size.name, OTU) %>%
  summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  ) %>%
  filter(OTU %in% sig_taxa$high_ab[1:n_display]) %>%
  mutate(OTU = factor(OTU, levels = sig_taxa$high_ab[1:n_display])) 

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

p <- ggplot(ASV_size, aes(x = Sample, y = Abundance, colour = OTU, shape = OTU)) +
  # group samples by size
  facet_grid(. ~ size.name, scales = "free_x", space = "free_x", switch = "x") +
  geom_point() + 
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
    axis.text.x = element_blank(),    # hide sample labels
    axis.ticks.x = element_blank(),
    strip.placement = "outside",      # place strips below the panel
    strip.text.x = element_text(size = 10, margin = margin(t = 5)),
    legend.text = element_markdown()  # italicize labels
  ) 

# Save plot
ggsave(fname, plot = p, width = 6, height = 4, dpi = 300)
