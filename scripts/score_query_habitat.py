#!/usr/bin/env python3
# import os
import pandas as pd
import numpy as np

# project = "AJ306297"

# base_path = "/Users/nanzhen/Documents/1_phd_2020/5_doctoral_program/0.11_16S_predict_202504/step4"
# blast_results  = os.path.join(base_path, f"{project}_blast_results_filtered_metadata.csv")
# source_weight  = os.path.join(base_path, "source_class_weight.csv")
# output_profile = os.path.join(base_path, f"{project}_habitat_profile.csv")

blast_results = snakemake.input[0]
source_weight = snakemake.input[1]
output_profile = snakemake.output[0]

# === 1. Load BLAST results (already has ncbi_isolation_source) ===
df_hits = pd.read_csv(blast_results, dtype=str)
df_hits["ncbi_isolation_source"] = (
    df_hits["ncbi_isolation_source"]
    .str.strip()
    .str.lower()
)

# === 2. Count how many hits per source, call it “hits” not “count” ===
hit_counts = (
    df_hits["ncbi_isolation_source"]
    .value_counts()
    .rename_axis("ncbi_isolation_source")
    .reset_index(name="hits")
)

# === 3. Load the source → category + weight table ===
df_lookup = pd.read_csv(source_weight, dtype=str)
df_lookup["ncbi_isolation_source"] = (
    df_lookup["ncbi_isolation_source"]
    .str.strip()
    .str.lower()
)
# drop the old “count” from lookup so it won’t conflict
if "count" in df_lookup.columns:
    df_lookup = df_lookup.drop(columns="count")
df_lookup["weight"] = df_lookup["weight"].astype(float)

# === 4. Merge hits with lookup ===
df_merged = hit_counts.merge(
    df_lookup,
    on="ncbi_isolation_source",
    how="left"
)

# === 5. Compute weighted hits and build the habitat profile ===
df_merged["weighted_hits"] = df_merged["hits"] * df_merged["weight"]
habitat_scores = (
    df_merged
    .groupby("Classification_v1", dropna=False)["weighted_hits"]
    .sum()
    .reset_index(name="score")
)
total = habitat_scores["score"].sum()
habitat_scores["proportion"] = habitat_scores["score"] / total
habitat_scores["proportion"] = habitat_scores["proportion"] * 100
habitat_scores["proportion"] = habitat_scores["proportion"].map(lambda x: f"{x:.2f}%")

# === 6. Save out ===
habitat_scores.to_csv(output_profile, index=False)
print(f"✅ Habitat profile written to: {output_profile}")