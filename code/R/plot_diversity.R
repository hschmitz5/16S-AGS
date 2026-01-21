source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
library(picante)

fname_div <- "./figures/diversity.png"

# Shannon with plot_richness
div_type = c("Shannon")

# Faith's phylogenetic diversity
metadata <- get_metadata(ps)

comm <- as.data.frame(as.matrix(ps@otu_table)) %>%
  t() 

phy <- ps@phy_tree

pd_results <- pd(comm, phy) %>%
  rownames_to_column(var = "Sample") %>%
  left_join(metadata, join_by(Sample))

# head(pd_results)

p1 <- plot_richness(ps, x = "size.name", measures = div_type) +
  geom_point() + 
  labs(
    x = "Size",
    y = paste(div_type,"index")
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_blank(), # removes default title
    axis.title.x = element_text(margin = margin(t = 10)) # adds space above x axis title
  )

p2 <- pd_results %>%
  ggplot(aes(x = size.name, y = PD)) +
  geom_point() +
  labs(
    x = "Size",
    y = "Faith's PD index",
  ) +
  theme_minimal(base_size = 10) +  # sets text size
  theme(
    axis.title.x = element_text(margin = margin(t = 10)) # same for second plot
  )

diversity_plot <- p1 + p2 
ggsave(fname_div, plot = diversity_plot, width = 8, height = 3, dpi = 300)
