# Functional enrichment (ORA + GSEA) of colorectal cancer meta-analysis DEGs
# Load packages
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)
library(tidyverse)
library(rio)

# org.Hs.eg.db (AnnotationDbi) masks dplyr verbs -- pin them back
select <- dplyr::select
filter <- dplyr::filter

# Pin dplyr verbs (plyr, if attached, masks these)
mutate <- dplyr::mutate; summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange; rename <- dplyr::rename; count <- dplyr::count
desc <- dplyr::desc; select <- dplyr::select; filter <- dplyr::filter

# Thresholds
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1

# Input / output
META_FILE <- "results/tables/meta-analysis/random_effect_model.csv"
OUT_CSV   <- "results/tables/enrichment"
OUT_FIG   <- "results/figures/enrichment"
PREFIX    <- "colorectal"
dir.create(OUT_CSV, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# Helper functions
# Ranked list for GSEA. Names = Entrez (Gene_ID), score = sign(log2FC)*-log10(p).
build_ranked_list <- function(meta_df) {
  df <- meta_df |>
    filter(!is.na(randomP), !is.na(randomSummary), randomP > 0, !is.na(Gene_ID)) |>
    mutate(rank_score = sign(randomSummary) * -log10(randomP),
           Gene_ID    = as.character(Gene_ID)) |>
    filter(!duplicated(Gene_ID)) |>
    arrange(desc(rank_score))
  setNames(df$rank_score, df$Gene_ID)
}

# ORA: GO (BP/MF/CC) + KEGG + Reactome
run_ora <- function(entrez_ids, universe_entrez) {
  if (length(entrez_ids) < 5) {
    message("    Too few genes for ORA (n=", length(entrez_ids), "), skipping.")
    return(NULL)
  }
  results <- list()
  
  for (ont in c("BP", "MF", "CC")) {
    results[[paste0("GO_", ont)]] <- tryCatch(
      enrichGO(gene = entrez_ids, universe = universe_entrez, OrgDb = org.Hs.eg.db,
               ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05,
               qvalueCutoff = 0.2, readable = TRUE),
      error = function(e) { message("    GO-", ont, " error: ", e$message); NULL }
    )
  }
  
  results[["KEGG"]] <- tryCatch(
    enrichKEGG(gene = entrez_ids, universe = universe_entrez, organism = "hsa",
               pAdjustMethod = "BH", pvalueCutoff = 0.05),
    error = function(e) { message("    KEGG error: ", e$message); NULL }
  )
  
  results[["Reactome"]] <- tryCatch(
    enrichPathway(gene = entrez_ids, universe = universe_entrez, organism = "human",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE),
    error = function(e) { message("    Reactome error: ", e$message); NULL }
  )
  
  results
}

# GSEA: GO-BP + KEGG + Hallmark
run_gsea <- function(ranked_list) {
  if (length(ranked_list) < 100) {
    message("    Too few ranked genes (n=", length(ranked_list), "), skipping GSEA.")
    return(NULL)
  }
  results <- list()
  
  results[["GSEA_GO_BP"]] <- tryCatch(
    gseGO(geneList = ranked_list, OrgDb = org.Hs.eg.db, ont = "BP",
          minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE),
    error = function(e) { message("    GSEA GO error: ", e$message); NULL }
  )
  
  results[["GSEA_KEGG"]] <- tryCatch(
    gseKEGG(geneList = ranked_list, organism = "hsa",
            minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE),
    error = function(e) { message("    GSEA KEGG error: ", e$message); NULL }
  )
  
  # msigdbr >= 10: use collection= (category deprecated); Entrez col is ncbi_gene
  hallmark <- msigdbr(species = "Homo sapiens", collection = "H") |>
    select(gs_name, ncbi_gene) |>
    mutate(ncbi_gene = as.character(ncbi_gene))
  
  results[["GSEA_Hallmark"]] <- tryCatch(
    GSEA(geneList = ranked_list, TERM2GENE = hallmark,
         minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE),
    error = function(e) { message("    GSEA Hallmark error: ", e$message); NULL }
  )
  
  results
}

# Export enrichment results to CSV
export_enrichment <- function(result_list, out_dir, prefix) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (name in names(result_list)) {
    res <- result_list[[name]]
    if (is.null(res)) next
    df <- tryCatch(as.data.frame(res), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    export(df, file.path(out_dir, paste0(prefix, "_", name, ".csv")))
  }
}

# Tidy long MSigDB / GO set names for plotting.
.pretty_set_label <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^GOBP_", "", x); x <- sub("^KEGG_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  acronyms <- c("DNA","RNA","TNFA","NFKB","IL1","IL2","IL6","JAK","STAT1","STAT3",
                "STAT5","KRAS","MTORC1","TGF","E2F","UV","WNT","MYC")
  vapply(strsplit(x, " ", fixed = TRUE), function(words) {
    hit <- toupper(words) %in% acronyms
    words[hit] <- toupper(words[hit])
    paste(words, collapse = " ")
  }, character(1))
}

# Dotplot. GSEA drawn on NES axis; ORA uses clusterProfiler dotplot.
save_dotplot <- function(res, title, path, width = 10, height = 8, show = 20) {
  if (is.null(res)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(res), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  
  is_gsea <- all(c("NES", "setSize") %in% colnames(df))
  if (is_gsea) {
    df <- df[order(df$p.adjust), , drop = FALSE]
    df <- head(df, show)
    df$Label <- stringr::str_wrap(.pretty_set_label(df$Description), 34)
    p <- ggplot(df, aes(x = NES, y = reorder(Label, NES))) +
      geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
      geom_segment(aes(x = 0, xend = NES, yend = reorder(Label, NES)),
                   color = "grey80", linewidth = 0.5) +
      geom_point(aes(size = setSize, fill = p.adjust), shape = 21, color = "grey30") +
      scale_fill_gradient(low = "#b2182b", high = "#2166ac", trans = "log10",
                          name = "p.adjust",
                          labels = function(v) format(v, scientific = TRUE, digits = 1)) +
      scale_size_continuous(range = c(3, 10), name = "Set size") +
      labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            axis.text.y = element_text(color = "black"),
            panel.grid.minor = element_blank())
  } else {
    p <- dotplot(res, showCategory = show) + ggtitle(title)
  }
  ggsave(path, plot = p, width = width, height = height, dpi = 600, bg = "white")
  invisible(p)
}

# GSEA enrichment plot (top n pathways)
save_gsea_plot <- function(res, title, path, n = 3) {
  if (is.null(res)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(res), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  ids <- head(df$ID, n)
  p <- gseaplot2(res, geneSetID = ids, title = title)
  ggsave(path, plot = p, width = 12, height = 6 * ceiling(n / 2), dpi = 300)
  invisible(p)
}

# ------------------------------------------------------------------------------
# Driver (single cancer: colorectal)
# ------------------------------------------------------------------------------
message("\n=== ", toupper(PREFIX), " ENRICHMENT ===")

meta <- import(META_FILE) |>
  mutate(Gene_ID = as.character(Gene_ID))
message("  Genes in meta: ", nrow(meta))

# REM output has no adjusted p-value; derive BH-adjusted padj from randomP.
meta <- meta |>
  mutate(rem_padj = p.adjust(randomP, method = "BH"))

# Universe: all tested genes with an Entrez ID
universe_entrez <- meta |> filter(!is.na(Gene_ID)) |> pull(Gene_ID) |> unique()

# Significant DEGs for ORA
sig <- meta |> filter(rem_padj < PADJ_CUTOFF, abs(randomSummary) >= LFC_CUTOFF)
message("  Significant DEGs (padj<", PADJ_CUTOFF, ", |LFC|>=", LFC_CUTOFF, "): ", nrow(sig))

entrez_up   <- sig |> filter(randomSummary > 0) |> pull(Gene_ID) |> unique()
entrez_down <- sig |> filter(randomSummary < 0) |> pull(Gene_ID) |> unique()
entrez_all  <- sig |> pull(Gene_ID) |> unique()
message("  Up: ", length(entrez_up), "  Down: ", length(entrez_down))

# ----- ORA -----
message("  Running ORA...")
ora_up   <- run_ora(entrez_up,   universe_entrez)
ora_down <- run_ora(entrez_down, universe_entrez)
ora_all  <- run_ora(entrez_all,  universe_entrez)

export_enrichment(ora_up,   OUT_CSV, paste0(PREFIX, "_ora_up"))
export_enrichment(ora_down, OUT_CSV, paste0(PREFIX, "_ora_down"))
export_enrichment(ora_all,  OUT_CSV, paste0(PREFIX, "_ora_all"))

for (name in names(ora_all)) {
  save_dotplot(ora_all[[name]], paste(PREFIX, name),
               file.path(OUT_FIG, paste0(PREFIX, "_ora_all_", name, "_dotplot.png")))
}

# ----- GSEA -----
message("  Running GSEA...")
ranked <- build_ranked_list(meta)
message("  Ranked genes for GSEA: ", length(ranked))

gsea_res <- run_gsea(ranked)
export_enrichment(gsea_res, OUT_CSV, paste0(PREFIX, "_gsea"))

for (name in names(gsea_res)) {
  save_dotplot(gsea_res[[name]], paste(PREFIX, name),
               file.path(OUT_FIG, paste0(PREFIX, "_", name, "_dotplot.png")))
  save_gsea_plot(gsea_res[[name]], paste(PREFIX, name, "- Top Pathways"),
                 file.path(OUT_FIG, paste0(PREFIX, "_", name, "_enrichplot.png")))
}
