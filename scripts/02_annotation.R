library(genekitr)
library(tidyverse)
library(rio)

# Pin dplyr verbs (plyr, if attached, masks these)
mutate <- dplyr::mutate; summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange; rename <- dplyr::rename; count <- dplyr::count
desc <- dplyr::desc; select <- dplyr::select; filter <- dplyr::filter

# Create output directory
dir.create("results/tables/annotated", showWarnings = FALSE, recursive = TRUE)

# Define all GEO dataset IDs
geo_ids <- c(
  "GSE89393", "GSE136630"
)

# Process each dataset
for (geo_id in geo_ids) {
  
  message("Processing ", geo_id, " ...")
  
  # Step 1: Load DESeq2 results
  input_file <- paste0("results/tables/DESeq2/", geo_id, ".csv")
  if (!file.exists(input_file)) {
    message("  Skipping ", geo_id, ": DESeq2 output not found")
    next
  }
  deseq_results <- import(input_file)
  
  # Step 2: Convert Gene_ID to character
  deseq_results <- deseq_results |>
    mutate(Gene_ID = as.character(Gene_ID))
  
  # Step 3: Get gene information from genekitr
  gene_info <- tryCatch(
    genInfo(id = deseq_results$Gene_ID, org = "hs", unique = TRUE, keepNA = FALSE),
    error = function(e) {
      message("  genInfo failed for ", geo_id, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(gene_info)) next
  
  # Step 4: Join gene information with DESeq2 results
  annotated <- deseq_results |>
    left_join(gene_info, by = c("Gene_ID" = "input_id"))
  
  # Step 5: Rename columns if present
  colnames(annotated)[colnames(annotated) == "symbol"]     <- "Gene_Symbol"
  colnames(annotated)[colnames(annotated) == "gene_name"]  <- "Gene_Description"
  colnames(annotated)[colnames(annotated) == "gene_biotype"] <- "Gene_Biotype"
  
  # Step 6: Keep only relevant columns
  annotated <- annotated |>
    select(any_of(c("Gene_ID", "Gene_Symbol", "Gene_Description", "Gene_Biotype",
                    "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")))
  
  # Step 7: Filter for protein-coding genes only
  if ("Gene_Biotype" %in% colnames(annotated)) {
    annotated <- annotated |>
      filter(Gene_Biotype == "protein_coding")
  }
  
  # Step 8: Export annotated results
  output_file <- paste0("results/tables/annotated/", geo_id, ".csv")
  export(annotated, output_file)
}
