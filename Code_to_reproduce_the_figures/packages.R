
# Packages required 
packages <- 
  c("tidyverse",
    "ggplot2",
    "ggrepel",
    "knitr",
    "ggforce",
    "biomaRt",
    "vroom",
    "stringr",
    "corrplot",
    "ggh4x",
    "openxlsx",
    "htmltools",
    "pheatmap",
    "purrr",
    "ggpubr", 
    "drc", 
    "progress", 
    "data.table", 
    # Bioconductor = BiocManager::install()
    "Biostrings",
    "limma")

packages_to_install <- 
  packages_cran[!packages %in% installed.packages()[, 1]]

print(packages_to_install)

message("Please intall the packages identified missing with either install.packages() or BiocManager::install().")

