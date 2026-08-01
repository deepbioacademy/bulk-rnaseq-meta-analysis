# Summary of meta-analysis results (colorectal cancer, tumor vs normal)
# Load packages
library(tidyverse)
library(rio)

# Pin dplyr verbs (plyr, if attached, masks these)
mutate <- dplyr::mutate; summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange; rename <- dplyr::rename; count <- dplyr::count
desc <- dplyr::desc; select <- dplyr::select; filter <- dplyr::filter

meta_dir <- "results/tables/meta-analysis"

# Up / down regulation summary from the Random Effect Model result
meta_result <- import(file.path(meta_dir, "random_effect_model.csv")) |>
  select(Gene_Symbol, randomSummary, randomP, Gene_Description) |>
  rename(log2FC = randomSummary, P.Value = randomP)

# Significance call: padj < 0.05 and |log2FC| > 1
meta_result <- meta_result |>
  mutate(Significance = case_when(
    P.Value < 0.05 & log2FC >  1 ~ "Up",
    P.Value < 0.05 & log2FC < -1 ~ "Down",
    TRUE                         ~ "NS"
  ))

# Counts and percentage of up- / down-regulated genes
gene_stats <- meta_result |>
  filter(Significance != "NS") |>
  group_by(Significance) |>
  summarise(
    Count      = n(),
    Percentage = n() / nrow(meta_result) * 100,
    .groups    = "drop"
  )

print(gene_stats)
export(gene_stats, file.path(meta_dir, "meta_degs_summary_stats.csv"))

# Significant DEGs that carry a gene symbol (annotated only)
annotated_genes <- meta_result |>
  mutate(Gene_Symbol = na_if(Gene_Symbol, "")) |>
  filter(!is.na(Gene_Symbol), Significance != "NS")

export(annotated_genes, file.path(meta_dir, "filtered_meta_degs_annotated_only.csv"))

# Regulation table from the filtered meta DEGs
meta_key_results <- import(file.path(meta_dir, "filtered_meta_degs.csv"), na.strings = "") |>
  select(Gene_ID, Gene_Symbol, Gene_Description,
         log2FoldChange = randomSummary, adjusted.P.Value = randomP) |>
  drop_na(Gene_Symbol) |>
  mutate(Regulation = ifelse(log2FoldChange > 0, "UP", "DOWN"))

export(meta_key_results, file.path(meta_dir, "meta_degs_regulation.csv"))

# ------------------------------------------------------------------------------
# 04. Cross-reference Venn intersection genes with meta-analysis effect sizes
venn_common <- "results/tables/venn/common_genes_all_datasets.csv"
comb_file   <- file.path(meta_dir, "meta_combining_mean.csv")

if (file.exists(venn_common) && file.exists(comb_file)) {
  intersect_ids <- import(venn_common) |> pull(Gene_ID) |> as.character()
  
  combined_meta <- import(comb_file) |>
    mutate(Gene_ID = as.character(Gene_ID))
  
  intersects <- combined_meta |> filter(Gene_ID %in% intersect_ids)
  
  export(intersects, "results/tables/venn/intersect_genes_meta_analysis.csv")
  message("Venn-intersect genes with meta stats: ", nrow(intersects))
} else {
  message("Missing Venn or combining-approach file; skipping cross-reference. Run 03_venn.R and 04_metavolcano.R first.")
}
