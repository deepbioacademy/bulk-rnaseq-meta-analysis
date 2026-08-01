# Volcano plot visualizing meta-analysis results (colorectal cancer)
# Load packages
library(tidyverse)
library(ggrepel)
library(rio)

# Pin dplyr verbs (plyr, if attached, masks these)
mutate <- dplyr::mutate; summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange; rename <- dplyr::rename; count <- dplyr::count
desc <- dplyr::desc; select <- dplyr::select; filter <- dplyr::filter

dir.create("results/figures/meta-analysis", showWarnings = FALSE, recursive = TRUE)

# Load Random Effect Model result
meta_result <- import("results/tables/meta-analysis/random_effect_model.csv") |>
  select(Gene_Symbol, randomSummary, randomP) |>
  rename(log2FC = randomSummary, P.Value = randomP)

# Significance thresholds: padj < 0.05 and |log2FC| > 1
meta_result <- meta_result |>
  mutate(Significance = case_when(
    P.Value < 0.05 & log2FC >  1 ~ "Up",
    P.Value < 0.05 & log2FC < -1 ~ "Down",
    TRUE                         ~ "NoSignificant"
  ))

# Top genes to label (annotated, significant, largest |log2FC|)
top_genes <- meta_result |>
  mutate(Gene_Symbol = na_if(Gene_Symbol, "")) |>
  filter(!is.na(Gene_Symbol), P.Value < 0.05, abs(log2FC) > 1) |>
  slice_max(order_by = abs(log2FC), n = 500)

# Symmetric x-axis limit from the data (avoids clipping large fold changes)
x_lim   <- ceiling(max(abs(meta_result$log2FC), na.rm = TRUE))
x_break <- if (x_lim > 10) 2 else 1

# Volcano plot
volcano <- ggplot(meta_result, aes(x = log2FC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.8, size = 2.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Down" = "#2c7fb8",
                                "NoSignificant" = "#636363",
                                "Up" = "#e34a33")) +
  theme_minimal() +
  labs(title = "", x = "log2 (Fold Change)", y = "-log10 (adjusted P-value)") +
  theme(
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.position = "right",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(colour = "black", hjust = 1, size = 12),
    axis.text.y  = element_text(colour = "black", hjust = 1, size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  geom_text_repel(
    data = top_genes, aes(label = Gene_Symbol),
    colour = "black", size = 3, max.overlaps = 10,
    direction = "both", max.time = 5,
    force = 5, force_pull = 5, point.padding = 0.5, seed = 40
  ) +
  scale_x_continuous(limits = c(-x_lim, x_lim),
                     breaks = seq(-x_lim, x_lim, by = x_break))

# Save
ggsave("results/figures/meta-analysis/volcano_plot.png",
       plot = volcano, width = 12, height = 12, dpi = 600, bg = "white")
