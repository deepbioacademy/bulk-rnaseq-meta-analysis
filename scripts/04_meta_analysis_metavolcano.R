# Meta-Analysis of RNA-seq data using MetaVolcanoR
# Colorectal cancer (tumor vs normal) across GEO datasets

# Load required packages
library(MetaVolcanoR)
library(tidyverse)
library(rio)

# Pin dplyr verbs (plyr, if attached, masks these)
mutate <- dplyr::mutate; summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange; rename <- dplyr::rename; count <- dplyr::count
desc <- dplyr::desc; select <- dplyr::select; filter <- dplyr::filter


# Output directories
dir.create("results/tables/meta-analysis",  showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/meta-analysis", showWarnings = FALSE, recursive = TRUE)

# Compatibility patch: MetaVolcanoR 1.0.1 x ggplot2 4.0
# ggplot2 4.0 plot objects are S7; MetaVolcanoR's S4 "MetaVolcano" class types
# its plot slots as "gg" and rejects them at construction (validObject error).
# Relax the plot slots to "ANY" so the object builds; metaresult is unaffected.
local({
  ns <- asNamespace("MetaVolcanoR")
  nm <- ".__C__MetaVolcano"
  if (bindingIsLocked(nm, ns)) unlockBinding(nm, ns)
  suppressWarnings(setClass("MetaVolcano",
                            representation(input = "data.frame", inputnames = "character",
                                           metaresult = "data.frame", MetaVolcano = "ANY", degfreq = "ANY"),
                            where = ns))
})

# Gene annotation lookup, reused from the pre-annotated tables produced by
# 02_annotation.R (results/tables/annotated). No re-annotation / network calls.
anno_files <- list.files("results/tables/annotated", pattern = "[.]csv$", full.names = TRUE)
gene_annotation <- anno_files |>
  lapply(function(f) {
    import(f) |> select(any_of(c("Gene_ID", "Gene_Symbol", "Gene_Description")))
  }) |>
  bind_rows() |>
  mutate(Gene_ID = as.character(Gene_ID)) |>
  distinct(Gene_ID, .keep_all = TRUE)

annotate_genes <- function(df, id_col = "Gene_ID") {
  df |>
    mutate(!!id_col := as.character(.data[[id_col]])) |>
    left_join(gene_annotation, by = id_col)
}

# Data pre-processing
# Auto-detect datasets from DESeq2 output
deg_files <- list.files("results/tables/DESeq2", pattern = "\\.csv$", full.names = TRUE)
geo_ids   <- tools::file_path_sans_ext(basename(deg_files))

# Read all studies into a named list
studies <- lapply(deg_files, import)
names(studies) <- geo_ids

# 04. Meta-Analysis -- Random Effect Model
meta_degs_rem <- rem_mv(
  diffexp       = studies,
  pcriteria     = "padj",
  foldchangecol = "log2FoldChange",
  genenamecol   = "Gene_ID",
  geneidcol     = NULL,
  collaps       = TRUE,
  vcol          = "lfcSE",
  cvar          = FALSE,
  metathr       = 0.01,
  jobname       = "MetaVolcano_REM",
  outputfolder  = "results/figures/meta-analysis/",
  draw          = "PDF",
  ncores        = 1 # adjust as per your PC core! 
)

meta_results <- meta_degs_rem@metaresult

# Annotate + export full REM result
annotated_results <- annotate_genes(meta_results)
export(annotated_results, "results/tables/meta-analysis/random_effect_model.csv")

# Filter DEGs: significant (randomP < 0.05) and reliable effect (|randomSummary| >= 1)
key_genes <- annotated_results |>
  filter(randomP < 0.05, abs(randomSummary) >= 1)
export(key_genes, "results/tables/meta-analysis/filtered_meta_degs.csv")

# Meta-Analysis -- Combining approach (Mean)
meta_degs_comb <- combining_mv(
  diffexp       = studies,
  pcriteria     = "padj",
  foldchangecol = "log2FoldChange",
  genenamecol   = "Gene_ID",
  metafc        = "Mean",
  metathr       = 0.01,
  collaps       = TRUE,
  jobname       = "MetaVolcano_Combining",
  outputfolder  = "results/figures/meta-analysis/",
  draw          = "PDF"
)

combined_meta <- meta_degs_comb@metaresult
export(annotate_genes(combined_meta),
       "results/tables/meta-analysis/meta_combining_mean.csv")

# Cross-reference with Venn intersection (genes DE across all datasets)
venn_file <- "results/tables/venn/common_genes_all_datasets.csv"
if (file.exists(venn_file)) {
  intersect_genes <- import(venn_file)
  
  intersects <- combined_meta |>
    filter(Gene_ID %in% as.character(intersect_genes$Gene_ID))
  
  intersec_results <- annotate_genes(intersects)
  export(intersec_results,
         "results/tables/meta-analysis/meta_combining_mean_intersect_genes.csv")
  message("Intersection genes annotated: ", nrow(intersec_results))
} else {
  message("Venn intersection file not found; skipping intersect cross-reference. Run 03_venn.R first.")
}
