# Common differentially expressed genes across colorectal cancer datasets
# Load required packages
library(tidyverse)
library(VennDiagram)
library(ggVennDiagram)
library(rio)

# Pin dplyr verbs (plyr, if attached, masks these)
mutate <- dplyr::mutate; summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange; rename <- dplyr::rename; count <- dplyr::count
desc <- dplyr::desc; select <- dplyr::select; filter <- dplyr::filter


# Create output directories
dir.create("results/figures/venn", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables/venn", showWarnings = FALSE, recursive = TRUE)

# Locate DEG tables
# Prefer annotated tables (contain Gene_Symbol); fall back to raw DESeq2 output.
input_dir <- if (length(list.files("results/tables/annotated", pattern = "\\.csv$")) > 0) {
  "results/tables/annotated"
} else {
  "results/tables/DESeq2"
}

deg_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
geo_ids   <- tools::file_path_sans_ext(basename(deg_files))
names(deg_files) <- geo_ids

# Load, filter for significant DEGs, and collect gene sets
sig_cutoff <- 0.05
lfc_cutoff <- 1

gene_sets <- map(deg_files, function(f) {
  import(f) |>
    filter(!is.na(padj), padj < sig_cutoff, abs(log2FoldChange) > lfc_cutoff) |>
    pull(Gene_ID) |>
    unique() |>
    as.character()
})

# Drop empty sets so the plotting/intersection code stays valid
gene_sets <- gene_sets[map_int(gene_sets, length) > 0]
n_sets <- length(gene_sets)

# Visualise overlaps
# A classic Venn is only legible up to ~5 sets. With more, use an UpSet plot.
if (n_sets >= 2 && n_sets <= 5) {
  
  fill_cols <- c("#99d8c9", "#addd8e", "#bcbddc", "#fec44f", "#fc9272")[seq_len(n_sets)]
  cat_cols  <- c("#2ca25f", "#31a354", "#756bb1", "#e6550d", "#de2d26")[seq_len(n_sets)]
  
  png("results/figures/venn/venn_diagram.png",
      units = "in", width = 7, height = 7, res = 600)
  
  venn.plot <- venn.diagram(
    x               = gene_sets,
    filename        = NULL,
    fill            = fill_cols,
    alpha           = 0.6,
    cex             = 0.8,
    cat.cex         = 1,
    cat.col         = cat_cols,
    margin          = 0.12,
    lwd             = 2,
    disable.logging = TRUE
  )
  grid::grid.draw(venn.plot)
  dev.off()
  
} else if (n_sets > 5) {
  
  # UpSet-style plot for many sets
  upset <- ggVennDiagram::ggVennDiagram(gene_sets, force_upset = TRUE)
  ggsave("results/figures/venn/upset_plot.png", upset,
         width = 12, height = 7, units = "in", bg = "white", dpi = 600)
}

# Compute every intersection region (works for any number of sets)
regions <- ggVennDiagram::process_region_data(ggVennDiagram::Venn(gene_sets))

# Full intersection table: region id, member datasets, and gene count
regions_out <- regions |>
  select(id, name, count) |>
  arrange(desc(count))
export(as.data.frame(regions_out),
       "results/tables/venn/intersection_summary.csv")

# Genes shared by ALL datasets (the region present in every set)
all_id <- paste(seq_len(n_sets), collapse = "/")
common_genes <- regions |>
  filter(id == all_id) |>
  pull(item) |>
  unlist()

common_df <- data.frame(Gene_ID = as.character(common_genes))
export(common_df, "results/tables/venn/common_genes_all_datasets.csv")

# DEG frequency table: how many datasets each gene is DE in (useful when
# a full N-way overlap is small). Ranked most-recurrent first.
deg_frequency <- gene_sets |>
  unlist() |>
  table() |>
  as.data.frame() |>
  setNames(c("Gene_ID", "n_datasets")) |>
  arrange(desc(n_datasets))

export(deg_frequency, "results/tables/venn/deg_frequency.csv")
