rm(list = ls())

library(readxl)
library(tidyverse)
library(MetBrewer)

fname_in  <- "./data/EPS_moduli.xlsx"
fname_out <- "./figures/PNPS.png"

taxa_pal <- "Hiroshige"

# define sample names
size <- data.frame(
  name = c("S", "M", "L", "XL", "XXL")
)

# Other Data
eps <- read_excel(fname_in, range = cell_cols("A:E")) 
eps$size.name = factor(eps$size.name, levels = size$name)

mu <- read_excel(fname_in, range = cell_cols("G:H")) 
mu$size.name = factor(mu$size.name, levels = size$name)

eps_long <- eps %>%
  pivot_longer(
    cols = -size.name,                              # all columns except size.mm
    names_to = c("extract.type", ".value"),        # split names into type + variable
    names_pattern = "(.*)\\.(avg|sd)"              # regex to match "LB.avg", "TB.sd", etc.
  ) %>%
  rename(
    pnps.avg = avg,
    pnps.sd  = sd
  )

ggplot(data = eps_long, aes(x = size.name, y = pnps.avg, color = extract.type)) +
  geom_point(position = position_dodge(width = 0.2), size = 2) +
  geom_line(aes(group = extract.type)) +
  geom_errorbar(
    aes(ymin = pnps.avg - pnps.sd, ymax = pnps.avg + pnps.sd),
    width = 0.2,
    position = position_dodge(width = 0.2)
  ) +
  scale_color_manual(values = met.brewer(taxa_pal, 2)) +
  labs(
    x = "Size",
    y = "Mean PN/PS",
    color = "Extract Type"
  ) +
  theme_minimal(base_size = 10)

ggsave(fname_out, height = 2, width = 4, dpi = 300)
