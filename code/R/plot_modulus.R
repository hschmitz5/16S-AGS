rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
library(readxl)
library(tidyverse)
library(patchwork)
library(MetBrewer)

fname_in  <- "./data/Rheometry_Nov_2024.xlsx"
fname_out <- "./figures/moduli.png"

# define sample names
size <- data.frame(
  name = c("S", "M", "L", "XL", "XXL")
)

#### Plot

p1 <- ggplot(modulus, aes(x = freq_rad, y = G_avg, color = size)) +
  geom_point() +
  geom_line(aes(group = size)) +
  geom_errorbar(
    aes(ymin = pmax(G_avg - G_sd, 0), ymax = G_avg + G_sd),
    width = 0.2
  ) +
  scale_color_manual(
    name = "Size", 
    values = met.brewer(size_pal, n_sizes)
  ) +
  labs(
    x = "Frequency [rad/s]",
    y = "Storage Modulus [Pa]",
  ) 

p2 <- ggplot(modulus, aes(x = freq_rad, y = G2_avg, color = size)) +
  geom_point() +
  geom_line(aes(group = size)) +
  geom_errorbar(
    aes(ymin = pmax(G2_avg - G2_sd, 0), ymax = G2_avg + G2_sd),
    width = 0.2
  ) +
  scale_color_manual(
    name = "Size",
    values = met.brewer(size_pal, n_sizes)
  ) +
  labs(
    x = "Frequency [rad/s]",
    y = "Loss Modulus [Pa]",
  ) 


# horizontal
p <- p1 + p2 + 
  plot_layout(guides = "collect") & 
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") &
  guides(
    color = guide_legend(
      title.position = "bottom",
      title.hjust = 0.5 # centers title
      )
  )

ggsave(fname_out, plot = p, width = 8, height = 3, dpi = 600)
