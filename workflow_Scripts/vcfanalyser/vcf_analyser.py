#!/usr/bin/env python3
# Usage: python plot_vep_variants.py <input_table> <output_folder>

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid")

table_file = sys.argv[1]
out_dir = sys.argv[2]
os.makedirs(out_dir, exist_ok=True)

# Load table
df = pd.read_csv(table_file, sep='\t')

# Extract CSQ first transcript
df['CSQ_FIRST'] = df['CSQ'].str.split(',').str[0]
csq_fields = df['CSQ_FIRST'].str.split('|', expand=True)

# Assign columns from VEP header
csq_cols = [
    "Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE",
    "EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position",
    "Amino_acids","Codons","Existing_variation","REF_ALLELE","UPLOADED_ALLELE","DISTANCE",
    "STRAND","FLAGS","SYMBOL_SOURCE","HGNC_ID","MANE","MANE_SELECT","MANE_PLUS_CLINICAL",
    "TSL","APPRIS","ENSP","REFSEQ_MATCH","REFSEQ_OFFSET","GIVEN_REF","USED_REF","BAM_EDIT",
    "SOURCE","SIFT","PolyPhen","HGVS_OFFSET","AF","CLIN_SIG","SOMATIC","PHENO","PUBMED",
    "MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS",
    "Pathway_KEGG_full","clinvar_clnsig","clinvar_hgvs","GO","ClinVar_SV","ClinVar_SV_CLNSIG",
    "ClinVar_SV_CLNACC"
]

csq_fields.columns = csq_cols[:csq_fields.shape[1]]  # match number of columns

# Add key columns to df
df['GENE'] = csq_fields['SYMBOL']
df['CONSEQUENCE'] = csq_fields['Consequence']
df['CLIN_SIG'] = csq_fields['CLIN_SIG']
df['AF'] = pd.to_numeric(csq_fields['AF'], errors='coerce')

print(df.columns.tolist())

# --------------------------
# 1️⃣ Top mutated genes
# --------------------------
# --------------------------
# Top 25 mutated genes
# --------------------------
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

# Ensure clean gene names
df['GENE'] = df['GENE'].astype(str).str.strip()
invalid_genes = ['', '.', '-', 'nan', 'None', 'NA', 'Unknown', ' ', 'null']
valid_genes = df[~df['GENE'].isin(invalid_genes)]

# Top 25 genes
top25_genes = valid_genes['GENE'].value_counts().head(25)

# --- Plot aesthetics ---
plt.figure(figsize=(10, 6))
sns.set(style="whitegrid")
palette = sns.color_palette("Spectral", n_colors=len(top25_genes))

sns.barplot(
    x=top25_genes.values,
    y=top25_genes.index,
    palette=palette,
    edgecolor='black'
)

# Add text labels
for i, v in enumerate(top25_genes.values):
    plt.text(v + 0.1, i, str(v), color='black', va='center', fontsize=11)

plt.xlabel("Number of Variants", fontsize=13)
plt.ylabel("Gene", fontsize=13)
plt.title("Top 25 Mutated Genes", fontsize=15, weight='bold')
plt.xlim(0, max(top25_genes.values) + 1.5)

plt.tight_layout()
plt.savefig(f"{out_dir}/top25_mutated_genes.png", dpi=400)
plt.close()

# Export variants for top 25 genes
top25_df = df[df['GENE'].isin(top25_genes.index)]

# Export to .table
table_file = f"{out_dir}/top25_gene_variants.table"
top25_df.to_csv(table_file, sep='\t', index=False)

# Export to .xlsx
xlsx_file = f"{out_dir}/top25_gene_variants.xlsx"
top25_df.to_excel(xlsx_file, index=False)

print(f"✅ Exported {len(top25_df)} variants from top 25 genes to:")
print(f"   • {table_file}")
print(f"   • {xlsx_file}")

#----------
# 2️⃣ Top consequences
top_cons = df['CONSEQUENCE'].value_counts().head(15)
plt.figure(figsize=(12,8))
sns.barplot(x=top_cons.values, y=top_cons.index, palette='magma')
plt.xlabel("Number of Variants", fontsize=14)
plt.ylabel("Consequence", fontsize=14)
plt.title("Top 15 Consequences", fontsize=16)
plt.tight_layout()
plt.savefig(f"{out_dir}/top_consequences.png", dpi=400)
plt.close()

# --------------------------
# 3️⃣ Variants per chromosome 
chrom_counts = df['CHROM'].value_counts().sort_index()
plt.figure(figsize=(22,8))
sns.barplot(x=chrom_counts.index, y=chrom_counts.values, palette='cubehelix')
plt.xlabel("Chromosome", fontsize=16)
plt.ylabel("Number of Variants", fontsize=16)
plt.title("Variants per Chromosome", fontsize=18)
plt.xticks(rotation=90, ha='center', fontsize=12)
plt.tight_layout()
plt.savefig(f"{out_dir}/variants_per_chromosome.png", dpi=400)
plt.close()

# --------------------------
# 4️⃣ Ti/Tv ratio
# --------------------------
# Define transition and transversion pairs
transitions = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
transversions = [
    ('A', 'C'), ('C', 'A'), ('A', 'T'), ('T', 'A'),
    ('G', 'C'), ('C', 'G'), ('G', 'T'), ('T', 'G')
]

# Count transitions and transversions
ti_count = df.apply(lambda row: (row['REF'], row['ALT']) in transitions, axis=1).sum()
tv_count = df.apply(lambda row: (row['REF'], row['ALT']) in transversions, axis=1).sum()

ti_tv_ratio = ti_count / tv_count if tv_count != 0 else float('nan')

print(f"Ti/Tv ratio: {ti_tv_ratio:.2f} (Ti={ti_count}, Tv={tv_count})")

# Optional: Plot Ti vs Tv counts
import matplotlib.pyplot as plt

plt.figure(figsize=(6,6))
plt.bar(['Transitions', 'Transversions'], [ti_count, tv_count], color=['skyblue','salmon'])
plt.ylabel("Number of Variants", fontsize=14)
plt.title(f"Ti/Tv Ratio = {ti_tv_ratio:.2f}", fontsize=16)
plt.tight_layout()
plt.savefig(f"{out_dir}/TiTv_ratio.png", dpi=400)
plt.close()

