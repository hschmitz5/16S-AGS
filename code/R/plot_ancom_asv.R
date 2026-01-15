rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")
library(writexl)
library(ComplexHeatmap)

write2excel <- FALSE
n_display_low <- 30

fname_excel <- "./data/ANCOM_ASV.xlsx"

fname_high  <- "./figures/DA_ancom_high.png"
fname_low   <- "./figures/DA_ancom_low.png"

# Cell height in inches (adjust as needed)
cell_h <- 0.2
cell_w <- 0.6 # same as cell_h

# Font sizes
row_fontsize <- 10
col_fontsize <- 11

# Define matrix of log-fold change data per size
process_lfc <- function(df, taxa) {
  df %>%
    filter(OTU %in% taxa) %>%
    # change value to zero if not significant
    mutate(
      M   = ifelse(q_size.nameM   < p_threshold & passed_ss_size.nameM,   lfc_size.nameM, 0),
      L   = ifelse(q_size.nameL   < p_threshold & passed_ss_size.nameL,   lfc_size.nameL, 0),
      XL  = ifelse(q_size.nameXL  < p_threshold & passed_ss_size.nameXL,  lfc_size.nameXL, 0),
      XXL = ifelse(q_size.nameXXL < p_threshold & passed_ss_size.nameXXL, lfc_size.nameXXL, 0)
    ) %>%
    dplyr::select(OTU, M, L, XL, XXL) %>% # Species - new name
    mutate(
      mean_lfc = rowMeans(across(M:XXL, abs), na.rm = TRUE)
    ) %>%
    rename_with(~ paste0(., "-S"), c("M", "L", "XL", "XXL")) %>%
    arrange(desc(mean_lfc))
}

# names of significant taxa
sig_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel, fname_excel)

# Read in data and add new names
output <- readRDS(ancom_fname)
res_prim <- output$res %>% 
  rename(OTU = taxon)

# Apply function to high and low abundance taxa
lfc_high <- process_lfc(res_prim, sig_taxa$high_ab)
lfc_low  <- process_lfc(res_prim, sig_taxa$low_ab)

### For plotting
fig_high <- lfc_high %>%
  dplyr::select(-mean_lfc) %>%
  tibble::column_to_rownames("OTU") %>%
  as.matrix()

fig_low <- lfc_low %>%
  head(n_display_low) %>%
  dplyr::select(-mean_lfc) %>%
  tibble::column_to_rownames("OTU") %>%
  as.matrix()

# ---- Plotting

create_heatmap <- function(mat, rowname_w = NULL, col_title = NULL) {
  n_cols <- ncol(mat)
  n_rows <- nrow(mat)
  
  args <- list(
    mat,
    name = "log fold\nchange",
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

ht_low  <- create_heatmap(fig_low, rowname_width_in, paste0(n_display_low, " largest lfc"))

# Save images
common_width <- ncol(fig_high) * cell_w + rowname_width_in + 2 # in
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