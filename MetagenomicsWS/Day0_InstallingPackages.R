install.packages(c('readxl', 'tidyverse', 'dplyr', 'vegan'))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c('phyloseq', 'microbiome'))
