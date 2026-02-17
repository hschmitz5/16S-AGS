suppressPackageStartupMessages({
  # ---- Packages ----
  library(readxl)
  library(phyloseq)
  library(tidyverse)
  # formatting figures
  library(patchwork)
  # colors
  library(RColorBrewer)
  library(MetBrewer) # fun color palettes
})
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, cache = TRUE)