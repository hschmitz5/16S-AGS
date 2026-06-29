rm(list = ls())
library(tidyverse)
library(patchwork)

# File names for concentration data
fname_pn    <- paste0("./data/EPS/PN_conc.rds")
fname_polys <- paste0("./data/EPS/PS_conc.rds")

# Calculate average and std of replicates
group_data <- function(fname) {
  df <- readRDS(fname) %>%
    filter(size != "XS") %>%
    group_by(size, extract) %>%
    summarize(
      avg = mean(C_VSS),
      sd = sd(C_VSS),
      .groups = "drop"
    )
  }
# Apply function to each assay
PN <- group_data(fname_pn) 
PS <- group_data(fname_polys)

# Calculate PN + PS and PN/PS
df_wide <- left_join(
  PN %>% select(size, extract, PN_avg = avg, PN_sd = sd), 
  PS %>% select(size, extract, PS_avg = avg, PS_sd = sd), 
  by = c("size", "extract")
  ) %>%
  mutate(
    total = PN_avg + PS_avg,
    PNPS = PN_avg/PS_avg,
    sd = NA
    ) 

# Combine into single data frame
df <- bind_rows(
  'Protein (PN)' = PN,
  'Polysaccharide (PS)' = PS,
  'PN + PS' = df_wide %>% select(size, extract, avg = total, sd),
  .id = "assay"
  ) %>%
  mutate(
    extract = recode(extract,"LB" = "Loosely Bound","TB" = "Tightly Bound"),
    assay = factor(assay, levels = c("Polysaccharide (PS)", "Protein (PN)", "PN + PS"))
    ) 

# Calculate PN/PS
PNPS <- df_wide %>% 
  select(size, extract, avg = PNPS) %>%
  mutate(
    extract = recode(extract,"LB" = "Loosely Bound","TB" = "Tightly Bound")
    ) 

# Determine maximum avg + sd
max_y <- ceiling(max(df$avg + df$sd))

# -------------------------------

# Make Plot

p <- ggplot(data = df, aes(x = size, y = avg, fill = assay)) +
  geom_col(position = "dodge", width = 0.8) +
  geom_errorbar(
    aes(ymin = avg - sd, ymax = avg + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  facet_wrap(~extract) + 
  ylim(0, max_y) +
  labs(
    x = "Size",
    y = expression(paste(mu, "g/mgVSS")),
    fill = NULL, # legend titles
    color = NULL 
  ) +
  scale_fill_manual(
    values = c(
      "Polysaccharide (PS)" = "lavenderblush2",
      "Protein (PN)"        = "lightblue",
      "PN + PS"      = "gray"
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(
      colour = NA # facet label outline
    )
  ) 

annot <- ggplot(data = PNPS) +
  geom_tile(aes(x = size, y = 1, fill = round(avg,1))) +
  geom_text(aes(x = size, y = 1, label = round(avg,1))) +
  scale_fill_gradient(low="white", high="lightgray") +
  facet_wrap(~extract) +
  labs(
    y = expression(frac(PN,PS))
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text  = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(angle = 0, hjust = 0.5),
    axis.title.x  = element_blank(),
    axis.text = element_blank(),
    axis.ticks  = element_blank()
  )

p2 <- p / annot +
  plot_layout(heights = c(3, 1))

fname_out <- "./figures/EPS.png"
ggsave(fname_out, plot = p2, width = 6.5, height = 3, dpi = 300)
