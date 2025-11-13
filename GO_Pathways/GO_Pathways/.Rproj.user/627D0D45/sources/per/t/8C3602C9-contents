library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(dplyr)
library(KEGGREST)

# 2. Load your DEG data
deg <- read.csv("Nadine_GBS_3Col.csv")  # columns: symbol, log2FC, pvalue

# Standardize column names
colnames(deg) <- tolower(colnames(deg))
colnames(deg) <- c("symbol", "log2FC", "pvalue")

# Convert gene symbols to Entrez IDs (required by KEGG)
deg_entrez <- bitr(deg$symbol,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# Merge back to main DEG data
deg_merged <- merge(deg, deg_entrez, by.x = "symbol", by.y = "SYMBOL")

# 3. Prepare input vector for Pathview
gene_data <- deg_merged$log2FC
names(gene_data) <- deg_merged$ENTREZID

# 4. Define pathway info
pathway_id <- "hsa04010"
pathway_info <- keggGet(pathway_id)[[1]]
pathway_name <- gsub(" ", "_", pathway_info$NAME)  # replace spaces with underscores

# 5. Create folder for pathway outputs
if(!dir.exists(pathway_name)) dir.create(pathway_name)

# 6. Run Pathview for the full pathway, output goes into folder
pathview(
  gene.data = gene_data,
  pathway.id = pathway_id,
  species = "hsa",
  kegg.native = TRUE,
  same.layer = TRUE,
  limit = list(gene = c(-3.5, 3.5)),
  low = list(gene = "blue"),
  mid = list(gene = "white"),
  high = list(gene = "red"),
  out.suffix = "Response_to_estrogen",
  out.dir = pathway_name  # save PNG/PDF in folder
)

# 7. Export KEGG gene mapping with DEG info
kegg_genes <- pathway_info$GENE
kegg_df <- data.frame(
  entrez_id = kegg_genes[seq(1, length(kegg_genes), 2)],
  info = kegg_genes[seq(2, length(kegg_genes), 2)],
  stringsAsFactors = FALSE
)
kegg_df$symbol <- sapply(strsplit(kegg_df$info, ";"), `[`, 1)
kegg_df$desc <- sapply(strsplit(kegg_df$info, "; "), `[`, 2)

deg_kegg <- left_join(kegg_df, deg_merged, by = c("symbol" = "symbol"))
deg_kegg_present <- deg_kegg[!is.na(deg_kegg$log2FC), ]

write.table(deg_kegg_present,
            file = file.path(pathway_name, "KEGG_pathway_gene_mapping.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

cat("KEGG pathway gene mapping exported to", file.path(pathway_name, "KEGG_pathway_gene_mapping.txt"), "\n")

# 8. Optional: restrict to only subset genes
pathway_genes <- c("GATA3","SERPINB9","PELP1","KMT2D","RARA","WBP2","ESR2",
                   "EP300","ARID5A","BRCA1","ESR1","ASH2L","KRT19","BCAS3")

deg_subset <- deg_merged[deg_merged$symbol %in% pathway_genes, ]
gene_data_subset <- deg_subset$log2FC
names(gene_data_subset) <- deg_subset$ENTREZID

pathview(
  gene.data = gene_data_subset,
  pathway.id = pathway_id,
  species = "hsa",
  kegg.native = TRUE,
  limit = list(gene = c(-3.5, 3.5)),
  low = list(gene = "blue"),
  mid = list(gene = "white"),
  high = list(gene = "red"),
  out.suffix = "Response_to_estrogen_subset",
  out.dir = pathway_name  # save PNG/PDF in folder
)