# PPI network analysis of colorectal cancer meta-analysis DEGs (STRING + igraph)
# Load packages
library(httr)
library(igraph)
library(ggraph)
library(ggplot2)
library(tidyverse)
library(rio)
library(scales)

# NOTE: the STRINGdb package downloads the full proteome (hundreds of MB), which
# is impractically slow here. Instead we query the STRING REST API for only our
# DEG set -- it returns just the edges among the submitted genes (small payload).

select <- dplyr::select
filter <- dplyr::filter

# Pin dplyr verbs (plyr, if attached, masks these)
mutate <- dplyr::mutate; summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange; rename <- dplyr::rename; count <- dplyr::count
desc <- dplyr::desc; select <- dplyr::select; filter <- dplyr::filter

# Thresholds
PADJ_CUTOFF     <- 0.05
LFC_CUTOFF      <- 1
STRING_SCORE     <- 700    # combined score threshold (0-1000)
STRING_SPECIES   <- 9606   # Homo sapiens
STRING_MAX_NODES <- 2000   # STRING API rejects networks > 2000 nodes
HUB_QUANTILE     <- 0.90   # top 10% degree = hub gene
MIN_MODULE_SIZE <- 5

# Input / output
META_FILE <- "results/tables/meta-analysis/random_effect_model.csv"
OUT_CSV   <- "results/tables/network"
OUT_FIG   <- "results/figures/network"
PREFIX    <- "colorectal"
dir.create(OUT_CSV, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# Helper functions
# Query STRING REST API for interactions among a set of gene symbols.
# Returns a data frame with preferredName_A / preferredName_B / score (0-1).
string_api_network <- function(symbols, required_score = STRING_SCORE,
                               species = STRING_SPECIES) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  if (length(symbols) < 2) return(NULL)
  
  resp <- tryCatch(
    httr::POST(
      "https://string-db.org/api/tsv/network",
      body = list(
        identifiers     = paste(symbols, collapse = "\n"),
        species         = species,
        required_score  = required_score,   # API scale 0-1000
        caller_identity = "crc_meta_pipeline"
      ),
      encode = "form",
      httr::timeout(120)
    ),
    error = function(e) { message("  STRING API error: ", e$message); NULL }
  )
  if (is.null(resp) || httr::http_error(resp)) {
    message("  STRING API request failed."); return(NULL)
  }
  
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  df <- tryCatch(utils::read.delim(text = txt, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || !all(c("preferredName_A", "preferredName_B", "score") %in% names(df))) {
    message("  STRING API returned no usable network."); return(NULL)
  }
  df
}

# Build igraph PPI network from STRING API interactions.
# lfc_map: named vector (Symbol -> pooled_log2FC) for node attributes/colouring.
build_ppi_graph <- function(net_df, lfc_map) {
  if (is.null(net_df) || nrow(net_df) == 0) return(NULL)
  
  edges <- net_df |>
    transmute(gene_from = preferredName_A, gene_to = preferredName_B,
              combined_score = score) |>
    filter(gene_from %in% names(lfc_map), gene_to %in% names(lfc_map))
  if (nrow(edges) == 0) return(NULL)
  
  nodes <- unique(c(edges$gene_from, edges$gene_to))
  vertex_df <- tibble(Symbol = nodes, pooled_log2FC = unname(lfc_map[nodes]))
  
  g <- graph_from_data_frame(d = edges, vertices = vertex_df, directed = FALSE)
  g <- simplify(g)
  
  # Drop isolated nodes (no interaction at this score threshold)
  iso <- which(degree(g) == 0)
  if (length(iso) > 0) g <- delete_vertices(g, iso)
  g
}

# Network centrality metrics
compute_metrics <- function(g) {
  tibble(
    Symbol      = V(g)$name,
    Degree      = degree(g),
    Betweenness = betweenness(g, normalized = TRUE),
    Closeness   = closeness(g, normalized = TRUE),
    Eigenvector = eigen_centrality(g)$vector,
    Log2FC      = V(g)$pooled_log2FC
  ) |> arrange(desc(Degree))
}

# Hub genes (top N% by degree)
get_hub_genes <- function(metrics_df, quantile_cut = HUB_QUANTILE) {
  thresh <- quantile(metrics_df$Degree, quantile_cut)
  metrics_df |> filter(Degree >= thresh)
}

# Community detection (Louvain)
detect_communities <- function(g) {
  comm <- cluster_louvain(g)
  V(g)$community <- membership(comm)
  list(graph = g, communities = comm)
}

# Visual attributes; label only top-N hubs so labels stay legible
set_visual_attrs <- function(g, metrics_df, label_top = 15) {
  deg <- degree(g)
  top_hubs <- names(sort(deg, decreasing = TRUE))[seq_len(min(label_top, length(deg)))]
  V(g)$size  <- rescale(deg, to = c(3, 16))
  V(g)$color <- ifelse(V(g)$pooled_log2FC > 0, "#de2d26", "#2171b5")
  V(g)$label <- ifelse(V(g)$name %in% top_hubs, V(g)$name, NA)
  g
}

# Full PPI network plot (base igraph)
save_network_plot <- function(g, title, path, seed = 42) {
  if (is.null(g) || vcount(g) < 2) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  png(path, width = 12, height = 12, units = "in", res = 300)
  set.seed(seed)
  plot(g, vertex.size = V(g)$size, vertex.color = V(g)$color,
       vertex.label = V(g)$label, vertex.label.cex = 0.65,
       vertex.label.color = "black", edge.width = 0.4, edge.color = "gray75",
       layout = layout_with_fr(g), main = title)
  legend("topright", legend = c("Upregulated", "Downregulated"),
         col = c("#de2d26", "#2171b5"), pch = 19, bty = "n", pt.cex = 1.5)
  dev.off()
  invisible(NULL)
}

# Publication-quality ggraph plot
save_ggraph_plot <- function(g, metrics_df, title, path, seed = 42) {
  if (is.null(g) || vcount(g) < 2) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  
  # Attributes must be keyed to V(g)$name (metrics_df is degree-sorted)
  deg <- degree(g)
  top_hubs <- names(sort(deg, decreasing = TRUE))[seq_len(min(15, length(deg)))]
  V(g)$deg       <- deg
  V(g)$direction <- ifelse(V(g)$pooled_log2FC > 0, "Up", "Down")
  V(g)$lab       <- ifelse(V(g)$name %in% top_hubs, V(g)$name, NA_character_)
  
  set.seed(seed)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.25, color = "gray60", width = 0.3) +
    geom_node_point(aes(size = deg, color = direction), alpha = 0.85) +
    geom_node_text(aes(label = lab), repel = TRUE, size = 3,
                   max.overlaps = Inf, na.rm = TRUE, fontface = "italic") +
    scale_color_manual(values = c(Up = "#de2d26", Down = "#2171b5"), name = "Direction") +
    scale_size_continuous(range = c(1.5, 9), name = "Degree") +
    labs(title = title) +
    theme_graph(base_family = "sans") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          legend.position = "right")
  
  ggsave(path, plot = p, width = 14, height = 12, dpi = 300)
  invisible(p)
}

# Hub gene barplot
save_hub_barplot <- function(hub_df, title, path, n = 25) {
  if (is.null(hub_df) || nrow(hub_df) == 0) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  df <- head(hub_df, n)
  p <- ggplot(df, aes(x = reorder(Symbol, Degree), y = Degree, fill = Log2FC)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(low = "#2171b5", mid = "white", high = "#de2d26",
                         midpoint = 0, name = "Log2FC") +
    coord_flip() +
    labs(title = title, x = NULL, y = "Degree Centrality") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(path, plot = p, width = 8, height = 7, dpi = 300)
  invisible(p)
}

# Module subgraph plots (top N modules)
save_module_plots <- function(g, fig_dir, prefix, n_modules = 5, seed = 42) {
  if (is.null(g) || is.null(V(g)$community)) return(invisible(NULL))
  module_sizes <- sort(table(V(g)$community), decreasing = TRUE)
  top_ids      <- names(head(module_sizes, n_modules))
  
  for (mid in top_ids) {
    vids <- which(V(g)$community == mid)
    subg <- induced_subgraph(g, vids = vids)
    if (vcount(subg) < MIN_MODULE_SIZE) next
    
    path <- file.path(fig_dir, paste0(prefix, "_module_", mid, "_network.png"))
    png(path, width = 8, height = 8, units = "in", res = 300)
    set.seed(seed)
    sub_deg <- degree(subg)
    plot(subg, vertex.size = rescale(sub_deg, to = c(5, 18)),
         vertex.color = V(subg)$color, vertex.label = V(subg)$name,
         vertex.label.cex = 0.6, vertex.label.color = "black",
         edge.width = 1, edge.color = "gray70", layout = layout_with_fr(subg),
         main = paste0(prefix, " - Module ", mid, " (", vcount(subg), " genes)"))
    dev.off()
  }
}

# Export metrics + hub genes + modules to CSV
export_network <- function(metrics_df, hub_df, modules_df, out_dir, prefix) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  export(metrics_df, file.path(out_dir, paste0(prefix, "_network_metrics.csv")))
  export(hub_df,     file.path(out_dir, paste0(prefix, "_hub_genes.csv")))
  export(modules_df, file.path(out_dir, paste0(prefix, "_network_modules.csv")))
}

# Driver (single cancer: colorectal)
meta <- import(META_FILE)
# REM output has no adjusted p-value; derive BH-adjusted padj from randomP.
meta <- meta |> mutate(rem_padj = p.adjust(randomP, method = "BH"))

# Significant DEGs (one row per symbol), ranked most-significant first so that,
# if we have to cap at the STRING 2000-node limit, the strongest signal is kept.
sig <- meta |>
  filter(rem_padj < PADJ_CUTOFF, abs(randomSummary) >= LFC_CUTOFF,
         !is.na(Gene_Symbol), Gene_Symbol != "") |>
  arrange(randomP, desc(abs(randomSummary))) |>
  distinct(Gene_Symbol, .keep_all = TRUE) |>
  transmute(Symbol = Gene_Symbol, pooled_log2FC = randomSummary)
message("  Significant DEGs: ", nrow(sig))

if (nrow(sig) < 5) stop("Too few significant DEGs for network analysis.")

if (nrow(sig) > STRING_MAX_NODES) {
  message("  Capping to top ", STRING_MAX_NODES,
          " DEGs by significance (STRING API 2000-node limit).")
  sig <- head(sig, STRING_MAX_NODES)
}

lfc_map <- setNames(sig$pooled_log2FC, sig$Symbol)

message("  Querying STRING API...")
net_df <- string_api_network(sig$Symbol)
if (is.null(net_df)) stop("STRING API returned no network.")
message("  Interactions returned: ", nrow(net_df))

message("  Building PPI network...")
g <- build_ppi_graph(net_df, lfc_map)
if (is.null(g)) stop("No STRING interactions found at this score threshold.")
message("  Nodes: ", vcount(g), "  Edges: ", ecount(g))

# Centrality + hubs
metrics   <- compute_metrics(g)
hub_genes <- get_hub_genes(metrics)
message("  Hub genes (top ", round((1 - HUB_QUANTILE) * 100), "%): ", nrow(hub_genes))

# Communities
comm_result <- detect_communities(g)
g <- comm_result$graph
modules_df <- tibble(Symbol = V(g)$name, Module = V(g)$community,
                     Degree = degree(g), Log2FC = V(g)$pooled_log2FC) |>
  arrange(Module, desc(Degree))

# Export CSVs
export_network(metrics, hub_genes, modules_df, OUT_CSV, PREFIX)

# Figures
g <- set_visual_attrs(g, metrics)
message("  Saving figures...")
save_network_plot(g, paste(tools::toTitleCase(PREFIX), "PPI Network"),
                  file.path(OUT_FIG, paste0(PREFIX, "_ppi_network.png")))
save_ggraph_plot(g, metrics, paste(tools::toTitleCase(PREFIX), "PPI Network"),
                 file.path(OUT_FIG, paste0(PREFIX, "_ppi_ggraph.png")))
save_hub_barplot(hub_genes, paste(tools::toTitleCase(PREFIX), "- Top Hub Genes"),
                 file.path(OUT_FIG, paste0(PREFIX, "_hub_genes_barplot.png")))
save_module_plots(g, OUT_FIG, PREFIX)
