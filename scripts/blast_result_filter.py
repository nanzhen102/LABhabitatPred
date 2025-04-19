"""
# Apply ProkAtlas-style filters:
# - percent identity >= 97
# - alignment length >= 150 bp
"""
import pandas as pd
# import os

# base_path = "/Users/nanzhen/Documents/1_phd_2020/5_doctoral_program/0.11_16S_predict_202504/step4"
# sample = "AJ306297"
# blast_file = os.path.join(base_path, f"{sample}_blast_results.tsv")
# filtered_output = os.path.join(base_path, f"{sample}_blast_results_filtered.tsv")

blast_file = snakemake.input[0]
filtered_output = snakemake.output [0]

# ==== Define BLAST outfmt 6 column names ====
columns = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

# ==== Load and filter the BLAST output ====
blast_df = pd.read_csv(blast_file, sep="\t", names=columns)

# Apply ProkAtlas-style filters:
# - percent identity >= 97
# - alignment length >= 150 bp
filtered_df = blast_df[
    (blast_df["pident"] >= 97) &
    (blast_df["length"] >= 150)
].copy()

# ==== Save filtered output ====
filtered_df.to_csv(filtered_output, sep='\t', index=False)

print(f"✅ Total input hits: {blast_df.shape[0]}")
print(f"✅ Filtered hits (≥97% identity, ≥150 bp): {filtered_df.shape[0]}")
print(f"📄 Saved to: {filtered_output}")