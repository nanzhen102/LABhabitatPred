import glob
import os

SAMPLES = [os.path.basename(f).replace(".fasta", "") for f in glob.glob("data/*.fasta")]

rule all:
    input:
        expand("results/habitat_profiles/{sample}_habitat_profile.csv", sample=SAMPLES)

# Step 1: Run blastn
rule blastn:
    input:
        query="data/{sample}.fasta",
        db="16S_database/LAB_16S_db.nsq"
    output:
        "results/blast_raw/{sample}_blast_results.tsv"
    log:
        "logs/{sample}_blast.log"
    shell:
        """
        blastn -query {input.query} \
               -db 16S_database/LAB_16S_db \
               -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
               -evalue 1e-5 \
               -out {output} \
               2> {log}
        """

# Step 2: Filter blast results
rule filter_blast:
    input:
        "results/blast_raw/{sample}_blast_results.tsv"
    output:
        "results/blast_filtered/{sample}_blast_results_filtered.tsv"
    log:
        "logs/{sample}_filter.log"
    script:
        "scripts/blast_result_filter.py"

# Step 3: Map filtered results to metadata
rule map_metadata:
    input:
        filtered="results/blast_filtered/{sample}_blast_results_filtered.tsv",
        metadata="16S_database/ssu_all_r220_Lactobacillaceae_deduplicated_1200bp_noN_matched_metadata.csv"
    output:
        matched="results/mapped_metadata/{sample}_blast_results_filtered_metadata.csv",
        unmatched="results/mapped_metadata/{sample}_blast_results_filtered_unmatched.csv"
    log:
        "logs/{sample}_map_metadata.log"
    script:
        "scripts/map_blast_to_metadata.py"

# Step 4: Score habitat preference
rule habitat_profile:
    input:
        mapped="results/mapped_metadata/{sample}_blast_results_filtered_metadata.csv",
        source_weight="16S_database/source_class_weight.csv"
    output:
        "results/habitat_profiles/{sample}_habitat_profile.csv"
    log:
        "logs/{sample}_habitat_profile.log"
    script:
        "scripts/score_query_habitat.py"