# Install CRAN Packages
install.packages("pak")

CRAN_REQUIRED <- c(
  "tidyverse", "rio", "ggplot2", "ggrepel", "ggpubr", "pheatmap", "patchwork",
  "igraph", "ggraph", "scales", "httr", "jsonlite", "Matrix", "metafor",
  "survival", "survminer", "glmnet", "randomForest", "pROC", "caret",
  "VennDiagram", "ggVennDiagram", "msigdbr", "gridExtra",
  "circlize", "future"
)
pak::pkg_install(CRAN_REQUIRED)


# Bioconductor Packages 
BIOC_REQUIRED <- c(
  "DESeq2", "sva", "edgeR", "limma", "RUVSeq", "clusterProfiler", "enrichplot",
  "org.Hs.eg.db", "ReactomePA", "GSVA", "recount3", "SummarizedExperiment",
  "MultiAssayExperiment", "curatedTCGAData", "ComplexHeatmap",
  "genekitr", "TCGAbiolinks", "biomaRt", "GEOquery", "ExperimentHub"
)
pak::pkg_install(BIOC_REQUIRED)

# Github Packages 
# Install MetaVolcanoR
devtools::install_github("iza-mcac/MetaVolcanoR")


# Optional: each enables a richer path, none is required.
BIOC_OPTIONAL <- c(
  "decoupleR",   # 18: TF activity via run_ulm instead of the built-in estimator
  "multiMiR",    # 17: experimentally validated miRNA-target pairs
  "depmap",      # 19: DepMap CRISPR dependencies through ExperimentHub
  "cBioPortalData",
  "SingleR", "celldex",   # 16: reference-based cell-type labels
  "slingshot",            # 16: pseudotime
  "AUCell"
)

pak::pkg_install(BIOC_OPTIONAL)
