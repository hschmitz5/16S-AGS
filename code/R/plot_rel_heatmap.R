rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")
library(ComplexHeatmap)
library(circlize)

fname_rel <- "./figures/rel_ab_heatmap.png"

# Cell height in inches (adjust as needed)
cell_h <- 0.2
cell_w <- 0.2 

# Font sizes
row_fontsize <- 10
col_fontsize <- 11

write2excel <- FALSE 

# Choose DA Taxa
ancom_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel = FALSE)
DA_taxa <- ancom_taxa$high_ab 

# Common taxa (ANCOM & ALDEx)
# DA_taxa2 <- get_common_taxa(ancom_fname, aldex_fname, p_threshold, effect_threshold)

# -------- Define groups ---------

high_ab_taxa <- get_rel_ASV(ps) %>%
  filter(Abundance > rel_ab_cutoff) %>%
  distinct(OTU) %>%
  pull(OTU)

rel_wide <- get_rel_ASV(ps) %>%
  filter(OTU %in% high_ab_taxa) %>%
  dplyr::select(OTU, Sample, Abundance) %>%  
  pivot_wider(
    names_from = Sample,
    values_from = Abundance
  ) %>%
  .[, c("OTU", sam_name)] %>%
  mutate(
    sig_taxa = ifelse(OTU %in% DA_taxa, "T", "F")
  )

combine_names <- function(data) {
  data %>%
    mutate(
      name_prefix = if_else(
        grepl("[0-9]$", OTU),       # last character is a number
        sub("-[^-]+$", "", OTU),    # remove last dash and everything after
        OTU                         # last character is a letter → keep full OTU
      )
    ) %>%
    group_by(sig_taxa, name_prefix) %>%
    summarise(
      n = n(),
      OTU_orig = first(OTU),  # only used when n = 1
      across(all_of(sam_name), \(x) sum(x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      OTU = if_else(
        n == 1,
        OTU_orig,
        paste0(name_prefix, "_", sig_taxa)
      )
    ) %>%
    select(OTU, sig_taxa, all_of(sam_name))
}

pseudo <- 1e-6  # choose based on detection limit

data_mat <- combine_names(rel_wide) %>%
  select(-sig_taxa) %>%
  column_to_rownames("OTU") %>%
  as.matrix() %>%
  { log10(. + pseudo) }

DA_taxa_renamed <- combine_names(rel_wide) %>%
  filter(sig_taxa == "T") %>%
  pull(OTU)

# ---- Plotting
n_cols <- ncol(data_mat)
n_rows <- nrow(data_mat)
split = rep(1:n_sizes, each = n_replicates)

row_fontface <- ifelse(rownames(data_mat) %in% DA_taxa_renamed, "bold", "plain")

ht_colors <- met.brewer(taxa_pal, type = "continuous")
ha_colors <- met.brewer(size_pal, n_sizes)

size_annot <- HeatmapAnnotation(sz = anno_block(gp = gpar(fill = ha_colors),
                                                labels = size$name,
                                                labels_gp = gpar(col = "white", fontsize = col_fontsize)))

ht <- Heatmap(
  data_mat,
  column_split = split,
  top_annotation = size_annot,
  #name = "Rel Ab [Log(%)]",
  cluster_columns = FALSE,
  show_heatmap_legend = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_title = NULL,
  col = ht_colors,
  width  = unit(n_cols * cell_w, "inches"),
  height = unit(n_rows * cell_h, "inches"),
  row_names_gp = gpar(fontsize = row_fontsize, fontface = row_fontface),
  column_names_gp = gpar(fontsize = col_fontsize)
)

lgd <- Legend(
  col_fun = colorRamp2(
    seq(min(data_mat), max(data_mat), length.out = length(ht_colors)),
    ht_colors
  ),
  title = "Rel Ab [Log(%)]",
  direction = "horizontal",
  title_position = "lefttop",
  at = pretty(c(min(data_mat), max(data_mat))),
  labels = pretty(c(min(data_mat), max(data_mat)))
)

# Draw combined heatmap
png(fname_rel,
    width = 5.5,  # width in inches; can adjust
    height = 9.5, # height in inches; can adjust
    units = "in", res = 300)
draw(ht)
draw(lgd, x = unit(0.41, "npc"), y = unit(0.99, "npc"), just = c("center", "top"))
dev.off()
