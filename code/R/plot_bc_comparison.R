# The bias corrected data are only comparable for individual taxa changes
# This means you cannot make a combined plot like a stacked bar plot 

source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")
source("./code/R/04_bias_correct.R")

fname_out <- "./figures/BC_comparison.png"
  
# names of significant taxa
all_taxa_ancom <- get_ancom_taxa(ancom_fname, p_threshold, write2excel = FALSE, fname_out = NULL)
high_ab_taxa <- high_ab_genera[high_ab_genera %in% all_taxa_ancom]

# Load Data
table_bc_long <- get_bc_abund(ancom_fname) %>%
  filter(Genus %in% high_ab_taxa) %>%
  mutate(Genus = factor(Genus, levels = high_ab_genera)) 

table_rel_long <- get_rel_genus(ps_ASV) %>%
  filter(Genus %in% high_ab_taxa) %>%
  dplyr::select(Sample, size.mm, size.name, Genus, Abundance) %>%
  mutate(Genus = factor(Genus, levels = high_ab_genera))

# Merge into one long-format dataframe
merged_long <- full_join(table_bc_long, table_rel_long, 
                         by = c("Sample", "Genus", "size.mm", "size.name")) %>%
  pivot_longer(cols = c(bc_abund, Abundance),
               names_to = "Metric", values_to = "Value")


# ---- Plotting -------

merged_scaled <- merged_long %>%
  group_by(Genus) %>%
  mutate(
    # Scale bc_abund to the range of Abundance
    Value_scaled = case_when(
      Metric == "bc_abund" ~ (Value - min(Value[Metric=="bc_abund"], na.rm=TRUE)) /
        (max(Value[Metric=="bc_abund"], na.rm=TRUE) -
           min(Value[Metric=="bc_abund"], na.rm=TRUE)) *
        (max(Value[Metric=="Abundance"], na.rm=TRUE) -
           min(Value[Metric=="Abundance"], na.rm=TRUE)) +
        min(Value[Metric=="Abundance"], na.rm=TRUE),
      TRUE ~ Value  # Abundance stays as-is
    )
  ) %>%
  ungroup()

my_colors <- met.brewer(taxa_pal, 2)

p <- ggplot(merged_scaled, aes(x = size.name, y = Value_scaled, color = Metric)) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_y_continuous(name = "Abundance") +   # single axis
  scale_color_manual(
    values = c("bc_abund" = my_colors[1], "Abundance" = my_colors[2]),
    labels = c("bc_abund" = "Bias-corrected abundance", 
               "Abundance" = "Relative abundance [%]")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top"
  ) +
  facet_wrap(~Genus, scales = "free_y", ncol = 3) +
  labs(x = "Size", 
       color = "Metric")

# Save the plot
ggsave(fname_out, p, width = 9, height = 8, dpi = 300)
