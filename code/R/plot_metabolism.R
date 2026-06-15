# NOTE: Using all data, not just top n_show
# This is to verify that major metabolic groups are not being dropped

rm(list = ls())
library(patchwork)
source("./code/R/01_load_ps.R")
source("./code/R/02_join_rel_ab_and_function.R")

write2excel <- 0

metab_order <- c("GAO", "Filamentous", "AOB")
n_rows <- length(metab_order)

DA_df <- readRDS("./data/DA_metab_processed.rds") %>%
  mutate(
    metab = factor(metab, levels = metab_order)
    )

rel_ab_df <- join_rel_ab_and_function(ps) %>%
  filter(metab_val == "Positive",
         metab %in% metab_order) %>%
  mutate(
    metab = factor(metab, levels = metab_order)
    )


# ------------ Plot ------------------

min_y = floor(min(DA_df$lfc))
max_y = ceiling(max(DA_df$lfc))

p1 <- ggplot(DA_df, aes(x = size, y = lfc)) +
  geom_col(fill = "steelblue", width = 0.6) +
  facet_wrap(~metab, scales = "fixed", ncol = 1) +
  ylim(min_y, max_y) +
  labs(
    y = "Log fold-change (relative to S)",
    x = "Size"
    ) +
  theme_minimal(base_size = 12) 

p2 <- ggplot(data = rel_ab_df, 
             aes(x = size.name, y = mean_sum)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_sum - sd_sum, ymax = mean_sum + sd_sum),
    width = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~metab, scales = "free_y", ncol = 1) +
  labs(
    y = "Relative Abundance (%)",
    x = "Size"
  ) +
  theme_minimal(base_size = 12) 


p <- p1 + p2


# Save plot
fname <- "./figures/metabolism.png"
ggsave(fname, plot = p, width = 6.5, height = 6, dpi = 300)



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