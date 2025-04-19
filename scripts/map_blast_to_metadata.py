import pandas as pd
import re
# import os

# base_path = "/Users/nanzhen/Documents/1_phd_2020/5_doctoral_program/0.11_16S_predict_202504/step4"

# project = "AJ306297"

# blast_file = os.path.join(base_path, f"{project}_blast_results_filtered.tsv")
# metadata_file = os.path.join(base_path, "ssu_all_r220_Lactobacillaceae_deduplicated_1200bp_noN_matched_metadata.csv")
# output_file = os.path.join(base_path, f"{project}_blast_results_filtered_metadata.csv")
# unmatched_file = os.path.join(base_path, f"{project}_blast_results_filtered_unmatched.csv")

blast_file = snakemake.input[0]
metadata_file = snakemake.input[1]
output_file = snakemake.output [0]
unmatched_file = snakemake.output [1]

# ====== Load input files ======
blast_df = pd.read_csv(blast_file, sep="\t")
metadata_df = pd.read_csv(metadata_file)

# ====== Extract accession from sseqid ======
# e.g., RS_GCF_947381685.1~NZ_CANCWV010000002.1 → RS_GCF_947381685.1
blast_df["full_accession"] = blast_df["sseqid"].str.extract(r"^([A-Z]+_GCA_\d+\.\d+|[A-Z]+_GCF_\d+\.\d+)")


# Merge with metadata on full accession
merged_df = blast_df.merge(metadata_df, how="left", left_on="full_accession", right_on="accession")

# ====== Save full matched metadata ======
output_columns = [
    "qseqid", "sseqid", "pident", "length", "accession",
    "ncbi_isolation_source", "ncbi_country", "gtdb_taxonomy"
]
merged_df[output_columns].to_csv(output_file, index=False)

# ====== Print and save unmatched sseqids ======
unmatched = merged_df[merged_df["ncbi_isolation_source"].isna()]["sseqid"].unique()
pd.DataFrame({"unmatched_sseqid": unmatched}).to_csv(unmatched_file, index=False)

print(f"✅ Mapped metadata saved to: {output_file}")
print(f"⚠️  Unmatched sseqids saved to: {unmatched_file} ({len(unmatched)} entries)")