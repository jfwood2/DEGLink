# Load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)
library(ggplot2)
library(ggridges)
library(cowplot)

# === SHARED PATHWAY HEATMAP === #

top_n <- 200
top_deg1 <- gsea_deg1@result %>% arrange(p.adjust) %>% head(top_n)
top_deg2 <- gsea_deg2@result %>% arrange(p.adjust) %>% head(top_n)

# Truncate descriptions
top_deg1$Description <- substr(top_deg1$Description, 1, 50)
top_deg2$Description <- substr(top_deg2$Description, 1, 50)

# Find shared pathway descriptions
shared_pathways <- intersect(top_deg1$Description, top_deg2$Description)

# Filter GSEA results for shared pathways
shared_deg1 <- gsea_deg1@result %>%
  filter(Description %in% shared_pathways) %>%
  select(Description, NES) %>%
  mutate(Timepoint = "2-Day")

shared_deg2 <- gsea_deg2@result %>%
  filter(Description %in% shared_pathways) %>%
  select(Description, NES) %>%
  mutate(Timepoint = "5-Day")

# Combine and pivot to wide format
combined_nes <- bind_rows(shared_deg1, shared_deg2) %>%
  pivot_wider(names_from = Timepoint, values_from = NES) %>%
  column_to_rownames("Description")

# Optional: sort by 2-Day NES for order
combined_nes <- combined_nes[order(combined_nes$`2-Day`, decreasing = TRUE), ]

# Panel A: Heatmap
heatmap_plot <- pheatmap(
  combined_nes,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "A: 2 vs. 5 dpi GSEA Shared Pathway Comparison Heatmap",
  angle_col = 45,
  fontsize_row = 5,
  fontsize_col = 12,
  border_color = NA
)

# === C: Ridge Plots of Log2FC per Pathway === #

# Load DEG files (ensure the columns match as mentioned before)
deg_2d <- read.csv("2DayCVB_3Col_hgnc.csv")
deg_5d <- read.csv("5DayCVB_3Col_hgnc.csv")

# Join DEG results on hgnc_id
deg_df <- full_join(deg_2d, deg_5d, by = "hgnc_id", suffix = c("_2Day", "_5Day"))

# Function to extract logFC data for specific pathways
extract_gene_logfc <- function(df, logfc_column) {
  df %>%
    select(hgnc_id, logFC = !!sym(logfc_column)) %>%
    left_join(gsea_deg1@result %>% select(Description, core_enrichment), by = c("hgnc_id" = "core_enrichment")) %>%
    filter(!is.na(Description))
}

# Extract LogFC data for 2-Day and 5-Day
ridge_data_2d <- extract_gene_logfc(deg_df, "logFC_2Day")
ridge_data_5d <- extract_gene_logfc(deg_df, "logFC_5Day")

# Combine the data and plot ridge plots
combined_ridge_data <- bind_rows(
  ridge_data_2d %>% mutate(Timepoint = "2-Day"),
  ridge_data_5d %>% mutate(Timepoint = "5-Day")
)

# Plot ridge plots (Panel C)
ridge_plot <- ggplot(combined_ridge_data, aes(x = logFC, y = Description, fill = Timepoint)) +
  geom_density_ridges(scale = 0.9) +
  theme_ridges() +
  xlab("Log2 Fold Change") +
  ylab("Pathway") +
  ggtitle("C: Log2 Fold Change Distribution per Pathway") +
  scale_fill_manual(values = c("2-Day" = "skyblue", "5-Day" = "lightcoral"))

# === FINAL COMBINED FIGURE === #

# Combine Panel A (Heatmap) and Panel C (Ridge Plot)
plot_grid(heatmap_plot[[4]], ridge_plot, labels = c("A", "C"), ncol = 2)
