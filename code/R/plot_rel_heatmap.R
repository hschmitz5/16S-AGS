source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")
source("./code/R/04_subset_DA.R")
library(ComplexHeatmap)

fname_rel <- "./figures/rel_ab_heatmap.png"

# Cell height in inches (adjust as needed)
cell_h <- 0.2
cell_w <- 0.2 

# Font sizes
row_fontsize <- 10
col_fontsize <- 11

# Choose DA Taxa
all_taxa_ancom <- get_ancom_taxa(ancom_fname, p_threshold, write2excel = FALSE)
DA_taxa <- high_ab_genera[high_ab_genera %in% all_taxa_ancom]
# Common taxa (ANCOM & ALDEx)
# DA_taxa <- get_common_taxa(ancom_fname, aldex_fname, p_threshold, effect_threshold, high_ab_genera)

# -------- Define groups ---------

metadata <- data.frame(ps@sam_data) %>% .[sam_name, ] # reorder

# Define colors for size categories directly from metadata
size_colors <- setNames(met.brewer(size_pal, n_sizes), size$name)

# Create the annotation directly from metadata
bot_annot <- HeatmapAnnotation(
  Size = factor(metadata$size.name, levels = size$name),
  col = list(Size = size_colors)
)

# ---- Plotting

make_rel_heatmap <- function(df, annot = NULL, legend_name = "Rel Ab [%]") {
  data_mat <- df %>%
    column_to_rownames("Genus") %>%
    as.matrix() %>%
    .[, sam_name]  # reorder columns 
  
  n_cols <- ncol(data_mat)
  n_rows <- nrow(data_mat)
  row_fontface <- ifelse(rownames(data_mat) %in% DA_taxa, "bold", "plain")
  
  Heatmap(
    data_mat,
    name = legend_name,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    column_names_rot = 0,
    column_names_centered = TRUE,
    bottom_annotation = annot,
    col = met.brewer("Hokusai2", type = "continuous"),
    width  = unit(n_cols * cell_w, "inches"),
    height = unit(n_rows * cell_h, "inches"),
    row_names_gp = gpar(fontsize = row_fontsize, fontface = row_fontface),
    column_names_gp = gpar(fontsize = col_fontsize)
  )
}

rel_wide <- get_rel_genus(ps) %>%
  filter(Genus %in% high_ab_genera) %>%
  dplyr::select(Genus, Sample, Abundance) %>%  
  pivot_wider(
    names_from = Sample,
    values_from = Abundance
  ) 

rel_wide_top  <- rel_wide %>%
  filter(Genus %in% high_ab_genera[1:3])

rel_wide_rest <- rel_wide %>%
  filter(Genus %in% high_ab_genera[4:length(high_ab_genera)])

ht_top  <- make_rel_heatmap(rel_wide_top, NULL, legend_name = "Top 3\nRel Ab [%]")
ht_rest <- make_rel_heatmap(rel_wide_rest, bot_annot, legend_name = "Bottom\nRel Ab [%]")
ht_list <- ht_top %v% ht_rest

# Draw combined heatmap
png(fname_rel,
    width = 6,  # width in inches; can adjust
    height = 8, # height in inches; can adjust
    units = "in", res = 300)
draw(ht_list,
     heatmap_legend_side = "right",
     annotation_legend_side = "top")
dev.off()