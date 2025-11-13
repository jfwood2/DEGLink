# ============================================
# 2D UMAP of KEGG Pathways Based on DEG GSEA (gseKEGG)
# ============================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(umap)
library(tibble)
library(enrichplot)
library(openxlsx)


# --- User input: your DEG CSV file ---
deg <- read.csv("Nadine_GBS_3Col.csv")

# Standardize column names
colnames(deg) <- tolower(colnames(deg))
colnames(deg)[colnames(deg) == "gene_name"] <- "gene"
colnames(deg)[colnames(deg) == "log2foldchange"] <- "log2fc"

# --- Convert gene symbols to Entrez IDs (for KEGG) ---
deg$entrez <- mapIds(org.Hs.eg.db,
                     keys = deg$gene,
                     keytype = "SYMBOL",
                     column = "ENTREZID",
                     multiVals = "first")

deg <- deg %>% filter(!is.na(entrez))

# --- Prepare ranked gene list for GSEA ---
gene_list <- deg$log2fc
names(gene_list) <- deg$entrez
gene_list <- sort(gene_list, decreasing = TRUE)

# --- Run GSEA KEGG ---
gsea_kegg <- gseKEGG(
  geneList     = gene_list,
  organism     = 'hsa',
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# --- Extract results ---
gsea_res <- as.data.frame(gsea_kegg@result)
cat("Number of significant KEGG pathways:", nrow(gsea_res), "\n")

# --- Build pathway × gene score matrix ---
# Each pathway has a gene set -> extract from core enrichment genes
kegg2gene <- strsplit(gsea_res$core_enrichment, "/")
names(kegg2gene) <- gsea_res$Description

genes <- unique(unlist(kegg2gene))
pathway_matrix <- matrix(0, nrow = length(kegg2gene), ncol = length(genes),
                         dimnames = list(names(kegg2gene), genes))

for (i in seq_along(kegg2gene)) {
  pathway_matrix[i, kegg2gene[[i]]] <- 1
}

# Map log2FC values to Entrez IDs
fc_map <- deg$log2fc
names(fc_map) <- deg$entrez

# Compute weighted pathway signatures
common_genes <- intersect(colnames(pathway_matrix), names(fc_map))
pathway_matrix <- pathway_matrix[, common_genes, drop = FALSE]
weighted_mat <- sweep(pathway_matrix, 2, fc_map[common_genes], "*")

# --- Use NES (Normalized Enrichment Score) for coloring ---
pathway_scores <- gsea_res$NES
names(pathway_scores) <- gsea_res$Description

# --- Sanity check before UMAP ---
cat("Number of KEGG pathways:", nrow(weighted_mat), "\n")
cat("Number of overlapping genes:", length(common_genes), "\n")

if (nrow(weighted_mat) == 0 | ncol(weighted_mat) == 0) {
  stop("Error: No overlapping genes/pathways found. Check gene ID formats (Entrez vs SYMBOL).")
}

# --- Run UMAP ---
set.seed(123)
umap_res <- umap(weighted_mat)

umap_df <- data.frame(umap_res$layout) %>%
  rownames_to_column("Pathway") %>%
  mutate(NES = pathway_scores[Pathway],
         Direction = ifelse(NES > 0, "Upregulated", "Downregulated"))

# --- Plot 2D UMAP ---
ggplot(umap_df, aes(x = X1, y = X2, color = NES)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 14) +
  labs(
    title = "2D UMAP of KEGG Pathways Based on GSEA (gseKEGG)",
    subtitle = "Pathways positioned by similarity of DEG profiles; NES color indicates direction",
    color = "Normalized Enrichment Score (NES)"
  )

# Create an output data frame for export
umap_export <- umap_df %>%
  dplyr::select(Pathway, X1, X2, NES, Direction) %>%
  arrange(desc(NES))

# Save to Excel
output_file <- "KEGG_UMAP_Pathway_Coordinates.xlsx"
write.xlsx(umap_export, file = output_file, row.names = FALSE)

cat("✅ UMAP pathway coordinates exported to:", output_file, "\n")