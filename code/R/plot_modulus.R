rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
library(readxl)
library(tidyverse)
library(patchwork)
library(MetBrewer)

fname_in  <- "./data/Bio Granules Nov 2024.xlsx"
fname_out <- "./figures/moduli.png"

# define sample names
size <- data.frame(
  name = c("S", "M", "L", "XL", "XXL")
)

read_loss <- function(sheet) {
  # read + reshape avg
  avg <- read_excel(fname_in, sheet = sheet, range = cell_cols("A:F")) %>%
    rename(freq = `Frequency (rad/s)`) %>%
    pivot_longer(!freq, names_to = "size.name", values_to = "avg")
  
  # read + reshape sd
  sd <- read_excel(fname_in, sheet = sheet, range = cell_cols("H:M")) %>%
    rename(freq = `Frequency (rad/s)`) %>%
    pivot_longer(!freq, names_to = "size.name", values_to = "sd")
  
  # combine + clean up
  avg %>%
    left_join(sd, join_by(freq, size.name)) %>%
    mutate(size.name = factor(size.name, levels = size$name))
}

G1 <- read_loss("G1")
G2 <- read_loss("G2")

#### Plot

p1 <- ggplot(G1, aes(x = freq, y = avg, color = as.factor(size.name))) +
  geom_point() +
  geom_line(aes(group = size.name)) +
  geom_errorbar(
    aes(ymin = pmax(avg - sd, 0), ymax = avg + sd),
    width = 0.2
  ) +
  scale_color_manual(
    name = "Size (q = 0.0833)",
    values = met.brewer(size_pal, n_sizes)
  ) +
  labs(
    x = "Frequency [rad/s]",
    y = "Storage Modulus [Pa]",
  ) 

p2 <- ggplot(G2, aes(x = freq, y = avg, color = as.factor(size.name))) +
  geom_point() +
  geom_line(aes(group = size.name)) +
  geom_errorbar(
    aes(ymin = pmax(avg - sd, 0), ymax = avg + sd),
    width = 0.2
  ) +
  scale_color_manual(
    name = "Size (q = 0.0833)",
    values = met.brewer(size_pal, n_sizes)
  ) +
  labs(
    x = "Frequency [rad/s]",
    y = "Loss Modulus [Pa]",
  ) 


# vertical
# p <- p1 / p2 + plot_layout(guides = "collect") &
#   theme(legend.position = "right")

# horizontal
p <- p1 + p2 + 
  plot_layout(guides = "collect") & 
  theme_minimal(base_size = 12) +
  theme(legend.position = "top") &
  guides(
    color = guide_legend(
      title.position = "top",
      title.hjust = 0.5 # centers title
      )
  )

ggsave(fname_out, plot = p, width = 8, height = 3, dpi = 300)
