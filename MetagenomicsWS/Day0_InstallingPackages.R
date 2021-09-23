install.packages(c('readxl', 'tidyverse', 'dplyr', 'vegan', 'ggpubr'))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c('phyloseq', 'microbiome'))
