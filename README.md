 ![Pipeline DAG](./docs/images/logo.jpg)
LABhabitPred is a reproducible and scalable Snakemake pipeline for predicting environmental habitat preferences of 16S rRNA sequences, specifically targeting the *Lactobacillaceae* family.

## Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Input](#input)
- [Output](#output)
- [Options](#options)
- [Dependencies](#dependencies)
- [Example Data](#example-data)

## Introduction

The `LABhabitPred` Pipeline is designed for researchers who need to infer the environmental habitat preferences of *Lactobacillaceae* sequences using 16S rRNA gene data. The pipeline:

- Runs BLASTn against a curated LAB-specific 16S reference database.
- Filters BLAST results by percent identity (≥97%) and alignment length (≥150 bp).
- Maps sequences to environmental metadata (isolation sources, taxonomy).
- Scores habitat preferences using category-specific weighted scoring adapted from the ProkAtlas method.

## Installation

### Prerequisites

- **Conda/Mamba**: [Miniconda](https://docs.conda.io/en/latest/miniconda.html) is required. [Mamba](https://mamba.readthedocs.io/en/latest/) is recommended for faster dependency resolution.
- **Git**: To clone the repository.

### Steps

1. **Clone the repository:**

```bash
git clone https://github.com/nanzhen102/LABhabitatPre.git
cd LABhabitatPre
```

2. **Create Conda Environment:**

Install dependencies manually:

```bash
conda create -n lab_habitat -c bioconda snakemake pandas numpy blast biopython
conda activate lab_habitat
```

## Usage

From the main project directory, execute the pipeline with:

```bash
snakemake --cores 8
```

This command will:

- Process each sample from the `data/` directory.
- Run BLAST and filtering steps.
- Map metadata and generate habitat profiles.
- Log execution in the `logs/` directory.

## Input

**Data Files:**

- Place your query 16S rRNA FASTA sequences in the `data/` directory.
- Filename format: `sample1.fasta`, `sample2.fasta`, etc.

**Database Files:**

- LAB-specific BLAST database files should be placed in `16S_database/`.
- Ensure `ssu_all_r220_Lactobacillaceae_deduplicated_1200bp_noN_matched_metadata.csv` and `source_class_weight.csv` are available in this directory.

## Output

After pipeline completion, the `results/` directory will contain:

- **Raw BLAST results:** `results/blast_raw/sample_blast_results.tsv`
- **Filtered BLAST results:** `results/blast_filtered/sample_blast_results_filtered.tsv`
- **Metadata-mapped results:** `results/mapped_metadata/sample_blast_results_filtered_metadata.csv`
- **Habitat profiles:** `results/habitat_profiles/sample_habitat_profile.csv`

Logs for each rule are stored in the `logs/` directory.

## Options

- **Conda environment:**  
  Use `--use-conda` flag to activate Conda environments.
- **Cores:**  
  `--cores <N>` specifies the number of cores.

## Dependencies

- Snakemake
- Conda/Mamba
- BLAST+
- Python ≥3.8
- Pandas
- NumPy
- Biopython

## Example Data

Available in the `data/` directory.
