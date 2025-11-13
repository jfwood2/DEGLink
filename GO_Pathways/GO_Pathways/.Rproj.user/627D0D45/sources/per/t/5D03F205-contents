# =========================================
# Multi-pathway heatmap from user-defined pathways
# =========================================

library(dplyr)
library(pheatmap)
library(RColorBrewer)

# -------------------------------
# 1. Load DEG data
# -------------------------------
deg <- read.csv("Nadine_GBS_3Col.csv")  # columns: symbol, log2FC, pvalue
colnames(deg) <- tolower(colnames(deg))
colnames(deg) <- c("symbol", "log2FC", "pvalue")

# Create a named vector for fast lookup
deg_lookup <- setNames(deg$log2FC, deg$symbol)

# -------------------------------
# 2. Load pathways from CSV
# -------------------------------
# Assuming your CSV has 2 columns: Pathway_name, Genes
pathway_df <- read.csv("heatmap.csv", stringsAsFactors = FALSE, header = FALSE)
colnames(pathway_df) <- c("Pathway", "Genes")

# Convert comma-separated gene strings into vectors
pathways <- setNames(
  lapply(pathway_df$Genes, function(x) strsplit(x, ",")[[1]]),
  pathway_df$Pathway
)

# -------------------------------
# 3. Build a heatmap matrix
# -------------------------------
all_genes <- unique(unlist(pathways))  # all genes across pathways
heat_matrix <- matrix(NA, nrow = length(pathways), ncol = length(all_genes))
rownames(heat_matrix) <- names(pathways)
colnames(heat_matrix) <- all_genes

# Fill in log2FC values
for(pw in names(pathways)){
  genes <- pathways[[pw]]
  heat_matrix[pw, genes] <- deg_lookup[genes]
}

# Optional: replace NA with 0 if you want missing genes to appear white
heat_matrix_filled <- heat_matrix
heat_matrix_filled[is.na(heat_matrix_filled)] <- 0

# -------------------------------
# 4. Create diverging color palette centered at 0
# -------------------------------
max_val <- max(abs(heat_matrix_filled), na.rm = TRUE)
breaks <- seq(-max_val, max_val, length.out = 51)
colors <- colorRampPalette(c("blue", "white", "red"))(50)

# -------------------------------
# 5. Plot heatmap
# -------------------------------
pheatmap(
  heat_matrix_filled,
  color = colors,
  breaks = breaks,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  main = "Multi-pathway Heatmap (log2FC)",
  fontsize_row = 4,
  fontsize_col = 4
)