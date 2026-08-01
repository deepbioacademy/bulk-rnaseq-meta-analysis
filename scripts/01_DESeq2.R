# Load packages
library(tidyverse)
library(DESeq2)
library(sva)
library(rio)

# Pin dplyr verbs (plyr, if attached, masks these)
mutate <- dplyr::mutate; summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange; rename <- dplyr::rename; count <- dplyr::count
desc <- dplyr::desc; select <- dplyr::select; filter <- dplyr::filter

# Create output directories
dir.create("results/tables/DESeq2", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/PCA", showWarnings = FALSE, recursive = TRUE)

# Define all GEO dataset IDs
geo_ids <- c(
  "GSE89393", "GSE136630"
)

# Define color palette for PCA plots
custom_colors <- c("Cancer" = "#de2d26", "Normal" = "#2171b5",
                   "cancer" = "#de2d26", "normal" = "#2171b5",
                   "Tumor" = "#de2d26", "tumor" = "#de2d26")

# Process each dataset
for (geo_id in geo_ids) {
  
  # Step 1: Load count data
  count_file <- paste0("data/raw_counts/", geo_id, "_raw_counts.tsv")
  if (!file.exists(count_file)) {
    count_file <- paste0("data/raw_counts/", geo_id, "_raw_counts.csv")
  }
  count_data_raw <- import(count_file)
  colnames(count_data_raw)[1] <- "GeneID"
  
  # Step 2: Load metadata
  metadata <- import(paste0("data/metadata/", geo_id, "_metadata.csv"))
  
  # Step 3: Standardize column names to lowercase
  colnames(metadata) <- tolower(colnames(metadata))
  
  # Step 4: Prepare count matrix
  count_data <- count_data_raw |>
    filter(!duplicated(GeneID)) |>
    column_to_rownames("GeneID") |>
    as.matrix()
  mode(count_data) <- "integer"
  
  # Step 5: Match metadata with count data
  metadata <- metadata |>
    filter(.data$sample %in% colnames(count_data)) |>
    arrange(match(.data$sample, colnames(count_data)))
  
  # Step 6: Prepare sample information
  colData <- data.frame(
    condition = relevel(factor(ifelse(grepl("normal|control|adjacent|healthy|non.?tumou?r",
                                            tolower(as.character(metadata$condition))),
                                      "Normal", "Tumor")), ref = "Normal"),
    row.names = colnames(count_data)
  )
  
  # Step 7: Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = colData,
    design = ~ condition
  )
  
  # Step 8: Filter low count genes
  keep <- rowSums(counts(dds) >= 10) >= round(nrow(metadata) / 2)
  dds <- dds[keep, ]
  
  # Step 9: Estimate surrogate variables for batch correction
  mod <- model.matrix(~ condition, data = colData(dds))
  mod0 <- model.matrix(~ 1, data = colData(dds))
  svobj <- svaseq(counts(dds), mod, mod0)
  n.sv <- svobj$n.sv
  
  # Step 10: Add surrogate variables to colData
  for (i in seq_len(n.sv)) {
    colData(dds)[[paste0("SV", i)]] <- svobj$sv[, i]
  }
  
  # Step 11: Update design formula with SVs
  sv_terms <- paste0("SV", seq_len(n.sv), collapse = " + ")
  design(dds) <- as.formula(paste("~", sv_terms, "+ condition"))
  
  # Step 12: Run DESeq2 analysis
  dds_sva <- DESeq(dds)
  
  # Step 13: Extract results
  res <- results(dds_sva)
  
  # Step 14: Convert to dataframe and add Gene_ID
  res_df <- res |>
    as.data.frame() |>
    rownames_to_column("Gene_ID")
  
  # Step 15: Export differential expression results
  output_file <- paste0("results/tables/DESeq2/", geo_id, ".csv")
  export(res_df, output_file)
  
  # Step 16: Variance stabilizing transformation for PCA
  dds_original <- dds
  design(dds_original) <- ~ condition
  vsd_pre <- vst(dds_original, blind = TRUE)
  vsd_post <- vst(dds_sva, blind = FALSE)
  
  # Step 17: Generate PCA data
  pca_data_pre <- plotPCA(vsd_pre, intgroup = "condition", returnData = TRUE)
  pca_data_post <- plotPCA(vsd_post, intgroup = "condition", returnData = TRUE)
  
  # Step 18: Calculate variance explained
  percent_pre <- round(100 * attr(pca_data_pre, "percentVar"))
  percent_post <- round(100 * attr(pca_data_post, "percentVar"))
  
  # Step 19: Create pre-correction PCA plot
  pca_pre <- ggplot(pca_data_pre, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percent_pre[1], "% variance")) +
    ylab(paste0("PC2: ", percent_pre[2], "% variance")) +
    ggtitle(paste0(geo_id, " (Pre-Correction)")) +
    scale_color_manual(values = custom_colors) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.title = element_text(face = "bold", size = 10),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      legend.position = "right",
      aspect.ratio = 1
    )
  
  # Step 20: Create post-correction PCA plot
  pca_post <- ggplot(pca_data_post, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percent_post[1], "% variance")) +
    ylab(paste0("PC2: ", percent_post[2], "% variance")) +
    ggtitle(paste0(geo_id, " (Post-Correction)")) +
    scale_color_manual(values = custom_colors) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.title = element_text(face = "bold", size = 10),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      legend.position = "right",
      aspect.ratio = 1
    )
  
  # Step 21: Save PCA plots
  pca_dir <- paste0("results/figures/PCA/", geo_id)
  dir.create(pca_dir, showWarnings = FALSE, recursive = TRUE)
  
  ggsave(
    filename = paste0(pca_dir, "/", geo_id, "_pre.png"),
    plot = pca_pre,
    width = 5,
    height = 5,
    units = "in",
    bg = "white",
    dpi = 600
  )
  
  ggsave(
    filename = paste0(pca_dir, "/", geo_id, "_post.png"),
    plot = pca_post,
    width = 5,
    height = 5,
    units = "in",
    bg = "white",
    dpi = 600
  )
}
