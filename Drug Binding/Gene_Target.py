import pandas as pd
import os
import itertools

# === INPUT FILES ===
files = [
    '/Users/jakewood/Downloads/UTCOMLS/UTCOMLS Research/Dr. Robert Smith/Code Projects/Drug Binding/Erlotinib.xlsx',
    '/Users/jakewood/Downloads/UTCOMLS/UTCOMLS Research/Dr. Robert Smith/Code Projects/Drug Binding/Neratinib.xlsx',
    '/Users/jakewood/Downloads/UTCOMLS/UTCOMLS Research/Dr. Robert Smith/Code Projects/Drug Binding/Pelitinib.xlsx',
    '/Users/jakewood/Downloads/UTCOMLS/UTCOMLS Research/Dr. Robert Smith/Code Projects/Drug Binding/Tyrphostin AG-1478.xlsx'
]
sheet_name = 0

# === SET OUTPUT PATHS ===
script_dir = os.path.dirname(os.path.abspath(__file__))
summary_path = os.path.join(script_dir, "gene_overlap_summary.txt")
csv_path = os.path.join(script_dir, "shared_gene_targets.csv")

# === LOAD FILES ===
dfs = []
for f in files:
    df = pd.read_excel(f, sheet_name=sheet_name)
    df.columns = [str(col).strip() for col in df.columns]
    dfs.append(df)

# === EXTRACT GENE TARGETS ===
gene_sets = [set(df.iloc[:, 0]) for df in dfs]
shared_genes = set.intersection(*gene_sets)

# === WRITE SUMMARY ===
summary_lines = []

# Unique gene targets per file
for i, genes in enumerate(gene_sets):
    unique_genes = genes - shared_genes
    summary_lines.append(f"{os.path.basename(files[i])}: {len(unique_genes)} unique gene targets")

summary_lines.append(f"\nShared gene targets across all 4 files: {len(shared_genes)}\n")

# === PAIRWISE SHARED GENE COMPARISON ===
drug_names = [os.path.splitext(os.path.basename(f))[0] for f in files]
pairwise_matrix = pd.DataFrame(index=drug_names, columns=drug_names)

for i in range(len(dfs)):
    for j in range(len(dfs)):
        if i == j:
            pairwise_matrix.iloc[i, j] = "—"
        else:
            shared = len(gene_sets[i].intersection(gene_sets[j]))
            pairwise_matrix.iloc[i, j] = str(shared)

summary_lines.append("Shared gene targets between drugs (pairwise):\n")
summary_lines.append(pairwise_matrix.to_string())

# === TRIPLET SHARED GENE COMPARISON ===
summary_lines.append("\nShared gene targets between all combinations of 3 drugs:\n")
triplets = list(itertools.combinations(range(len(dfs)), 3))
for triplet in triplets:
    drugs_in_triplet = [drug_names[i] for i in triplet]
    shared_triplet = gene_sets[triplet[0]].intersection(gene_sets[triplet[1]], gene_sets[triplet[2]])
    summary_lines.append(f"{', '.join(drugs_in_triplet)}: {len(shared_triplet)}")

# === SAVE SUMMARY FILE ===
with open(summary_path, "w") as f:
    f.write('\n'.join(summary_lines))

print(f"📄 Summary written to: {summary_path}")

# === COLLECT PROBABILITIES AND COMMON NAMES ===
probability_dict = {}
common_name_map = {}

for i, df in enumerate(dfs):
    drug_name = os.path.splitext(os.path.basename(files[i]))[0]
    df_filtered = df[df.iloc[:, 0].isin(shared_genes)]
    
    # Build gene → probability map
    gene_to_prob = dict(zip(df_filtered.iloc[:, 0], df_filtered.iloc[:, 5]))
    probability_dict[drug_name] = gene_to_prob

    # Build gene → common name map, only once
    for idx, row in df_filtered.iterrows():
        gene = row.iloc[0]
        common = row.iloc[1]
        if gene not in common_name_map:
            common_name_map[gene] = common

# === BUILD FINAL TABLE ===
final_df = pd.DataFrame(probability_dict)
final_df.index.name = 'Gene Target'
final_df['Common Name'] = final_df.index.map(common_name_map)
final_df = final_df.reset_index()
final_df = final_df[['Gene Target', 'Common Name'] + [col for col in final_df.columns if col not in ['Gene Target', 'Common Name']]]

# === EXPORT FINAL PIVOTED TABLE ===
final_df.to_csv(csv_path, index=False)
print(f"✅ Shared gene binding probability matrix saved to: {csv_path}")
