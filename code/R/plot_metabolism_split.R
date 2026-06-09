# NOTE: Using all data, not just top n_show
# This is to verify that major metabolic groups are not being dropped

rm(list = ls())
source("./code/R/01_load_ps.R")
source("./code/R/02_metab_and_DA.R")

write2excel <- 0

# Define relative abundance
rel_size <- get_rel(ps) %>%
  group_by(Genus, size.name) %>%
  summarize(
    mean_ab = mean(Abundance),
    std_ab  = sd(Abundance),
    .groups = "drop"
  ) 

# Load metabolism data
# Input must contain Genus
m <- get_metabolism(rel_size, metab_fname)

# define taxa in each metabolism group
taxa_P <- map(m, ~ rownames(m)[which(.x == "P")])
taxa_V <- map(m, ~ rownames(m)[which(.x == "V")])
taxa_PV <- map(m, ~ rownames(m)[.x %in% c("P", "V")])

# function sums up values when taxa = P or V
summarize_metab <- function(taxa_list, value_col) {
  map_dfr(names(taxa_list), function(nm) {
    rel_size %>%
      filter(Genus %in% taxa_list[[nm]]) %>%
      group_by(size.name) %>%
      summarize(
        sum_mean = sum(mean_ab, na.rm = TRUE),
        sum_sd   = sum(std_ab, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(metab_val = value_col, metab = nm) %>%
      dplyr::select(metab_val, metab, size.name, sum_mean, sum_sd)
  }) 
}

df_P  <- summarize_metab(taxa_P, "P") 
df_V  <- summarize_metab(taxa_V, "V") 
# P + V: only sum metab found in V
df_PV <- summarize_metab(taxa_PV, "P + V") %>%
  filter(metab %in% unique(df_V$metab))


# joins data sets
df <- bind_rows(df_P, df_PV) %>%
  # define panel grouping
  mutate(
    panel = case_when(
      metab %in% c("GAO", "Nitrite reduction") ~ "Nitrite reduction & GAO",
      metab %in% c("PAO", "Filamentous") ~ "Filamentous & PAO",
      metab %in% c("AOB", "NOB") ~ "AOB & NOB"
    ),
    metab_val = recode(metab_val, "P" = "Positive", "P + V" = "Positive + Variable")
  ) %>%
  # reorder
  dplyr::select(panel, metab, metab_val, size.name, sum_mean, sum_sd)

# Convert to factor
df$size.name <- factor(df$size.name, levels = size$name)
df$metab <- factor(
  df$metab, levels = c("GAO", "PAO", "Nitrite reduction", "Filamentous", "AOB", "NOB")
)
df$panel <- factor(
  df$panel, levels = c("Nitrite reduction & GAO", "Filamentous & PAO", "AOB & NOB")
)

# ------------ Plot ------------------

p <- ggplot(data = df, 
            aes(x = size.name, y = sum_mean, linetype = metab_val, group = metab_val)) +
  geom_point(color = "steelblue") +
  geom_line(color = "steelblue") +
  geom_errorbar(
    aes(ymin = sum_mean - sum_sd, ymax = sum_mean + sum_sd),
    color = "steelblue",
    width = 0.2,
    position = position_dodge(width = 0.2)
  ) +
  facet_wrap(~metab, nrow = 3, scales = "free") +
  labs(
    x = "Size",
    y = "Relative Abundance [%]"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.key.spacing.x = unit(0.5, "in"),
    legend.title = element_blank()
  )

# Save plot
fname <- "./figures/metabolism.png"
ggsave(fname, plot = p, width = 6.5, height = 5, dpi = 300)



# ------ Write Data to Excel

if (write2excel == 1) {
  ### Do not exclude data
  # define relative abundance
  rel_wide <- get_rel_wide(ps) %>%
    rownames_to_column(var = "Genus")
  
  new_m <- get_metabolism(rel_wide, metab_fname) %>%
    # true if any metabolic groups in row are defined
    mutate(tf = as.integer(if_any(everything(), ~ !is.na(.x)))) %>%
    rownames_to_column(var = "Genus")
  
  full_df <- left_join(rel_wide, new_m, by = "Genus") %>%
    filter(tf == 1) %>%
    dplyr::select(-tf) %>%
    relocate(where(is.numeric), .after = where(is.character)) %>%
    arrange(Genus)
  
  library(writexl)
  write_xlsx(full_df, path = "./data/rel_ab_metab.xlsx")
}