source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")

# display top n most abundant genera
n_display <- 10

fname <- "./figures/abund_size_stacked.png"

# names of significant taxa
write2excel = FALSE
sig_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel, fname_excel)

ASV_size <- get_rel_ASV(ps) %>%
  dplyr::select(Sample, size.mm, size.name, OTU, Abundance) %>%
  # Rename genera outside top n_display as "Other"
  mutate(
    OTU = ifelse(is.na(OTU) | !(OTU %in% sig_taxa$high_ab[1:n_display]), "Other", OTU)
  ) %>%
  group_by(Sample, size.mm, size.name, OTU) %>%
  summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  ) %>%
  mutate(OTU = factor(OTU, levels = c(sig_taxa$high_ab[1:n_display], "Other")))
  

p <- ggplot(ASV_size, aes(x = Sample, y = Abundance, fill = fct_rev(OTU))) +
  # group samples by size
  facet_grid(. ~ size.name, scales = "free_x", space = "free_x", switch = "x") +  
  geom_col(width = 0.95) + # space between same size
  scale_fill_manual(
    values = c("gray", met.brewer(taxa_pal, n_display)),
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
    strip.text.x = element_text(size = 10, margin = margin(t = 5))
  )

# Save plot
ggsave(fname, plot = p, width = 6, height = 4, dpi = 300)
