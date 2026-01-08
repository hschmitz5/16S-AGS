source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")
library(writexl)
library(ComplexHeatmap)

write2excel      <- FALSE
n_display_low    <- 30
effect_threshold <- 1

fname_excel <- "../results/ALDEx2_taxa.xlsx"
fname_high  <- "./figures/DA_aldex_high.png"
fname_low   <- "./figures/DA_aldex_low.png"

# Cell height in inches (adjust as needed)
cell_h <- 0.2
cell_w <- 0.6 # same as cell_h

# Font sizes
row_fontsize <- 10
col_fontsize <- 11

sig_taxa <- get_aldex_taxa(aldex_fname, ps, p_threshold, effect_threshold, rel_ab_cutoff, write2excel, fname_out)

# Read in ALDEx data
output <- readRDS(aldex_fname) 
output$res <- map(output$res, ~ .x %>% 
                    rownames_to_column("OTU") %>%
                    left_join(sig_taxa$name_tbl, join_by(OTU))
                  ) 

process_clr <- function(df, taxa, p_threshold, effect_threshold) {
  map2(df$res, df$comparison, ~ .x %>%
         filter(Species_updated %in% taxa) %>%
         mutate(
           effect_update = if_else(
             wi.eBH < p_threshold & abs(effect) > effect_threshold,
             effect,
             0
           )
         ) %>%
         group_by(Species_updated) %>%
         summarise(
           !!as.character(.y) := mean(effect_update, na.rm = TRUE),
           .groups = "drop"
         )
  ) %>%
    reduce(full_join, by = "Species_updated") %>%
    mutate(
      mean_effect = rowMeans(across(all_of(df$comparison), abs), na.rm = TRUE)
    ) %>%
    arrange(desc(mean_effect))
}


# Apply function to high and low abundance taxa
effect_high <- process_clr(output, sig_taxa$high_ab, p_threshold, effect_threshold)
effect_low  <- process_clr(output, sig_taxa$low_ab, p_threshold, effect_threshold)

### For plotting
fig_high <- effect_high %>%
  dplyr::select(-mean_effect) %>%
  tibble::column_to_rownames("Species_updated") %>%
  as.matrix()

fig_low <- effect_low %>%
  head(n_display_low) %>%
  dplyr::select(-mean_effect) %>%
  tibble::column_to_rownames("Species_updated") %>%
  as.matrix()

# ---- Plotting

create_heatmap <- function(mat, rowname_w = NULL, col_title = NULL) {
  n_cols <- ncol(mat)
  n_rows <- nrow(mat)
  
  args <- list(
    mat,
    name = "effect",
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 0,
    column_names_centered = TRUE,
    column_title = col_title,
    # size
    width  = unit(n_cols * cell_w, "inches"),
    height = unit(n_rows * cell_h, "inches"),
    row_names_gp = gpar(fontsize = row_fontsize),
    column_names_gp = gpar(fontsize = col_fontsize)
  )
  # Only add row_names_max_width if it is not NULL
  if (!is.null(rowname_w)) {
    args$row_names_max_width <- unit(rowname_w, "inches")
  }
  
  do.call(Heatmap, args)
}

ht_high <- create_heatmap(fig_high, NULL, paste0("Abundance > ",rel_ab_cutoff,"%"))

# Define rowname width in first plot
ht_grob <- draw(ht_high, show_heatmap_legend = FALSE)
ht1 <- ht_grob@ht_list[[1]]
fig_props <- ht1@layout$layout_size

rowname_width_in <- convertWidth(
  fig_props$row_names_right_width,
  "inches", valueOnly = TRUE
)

ht_low  <- create_heatmap(fig_low, rowname_width_in, paste0(n_display_low, " largest effect"))

# Save images
common_width <- n_cols * cell_w + rowname_width_in + 2 # in
png(fname_high, 
    width = common_width,
    height = nrow(fig_high) * cell_h + 1,
    units = "in", res = 300)
draw(ht_high)
dev.off()

# png(fname_low, 
#     width = common_width,
#     height = nrow(fig_low) * cell_h + 1,
#     units = "in", res = 300)
# draw(ht_low)
# dev.off()