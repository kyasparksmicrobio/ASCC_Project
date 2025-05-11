#!/bin/bash

# === Create README.md for 16S ===
mkdir -p 16S
cat << 'EOF' > 16S/README.md
# 16S

This directory contains 16S rRNA amplicon analysis materials, including scripts, metadata, and output files.

## Key Subfolders

- `16S_Rscripts/`: R scripts for all stages of 16S analysis (alpha, beta, core taxa, MaAsLin2, taxonomy).
- `16S_metadata/`: Metadata tables used in QIIME2, MaAsLin2, and visualization steps.
- `16S_coremetrics/`: QIIME2 core diversity metrics outputs.
- `16S_alpha_diversity/`, `16S_beta_diversity/`: R-based diversity analysis outputs.
- `16S_Maaslin/`: R scripts and outputs from MaAsLin2 differential abundance analysis.
- `16S_MAGs_coniferous/`: Matched ASVs to MAGs via taxonomy.
- `16S_CoreTaxaAnalysis/`: Scripts for identifying and visualizing the core microbiome.
- `16S_taxonomy/`: Taxonomic bar plots and taxonomic summaries.
- `16S_slurm_scripts/`: SLURM scripts used for HPC processing.
- `16S_slurm_outputs/`: Job logs from HPC execution.

EOF

# === Create README.md for ITS ===
mkdir -p ITS
cat << 'EOF' > ITS/README.md
# ITS

This directory contains internal transcribed spacer (ITS) fungal community analysis data.

## Key Subfolders

- `ITS_Rscripts/`: All R scripts used for alpha, beta, and core taxa analyses.
- `ITS_metadata/`: Sample metadata and experimental design tables.
- `ITS_coremetrics/`: QIIME2 diversity outputs.
- `ITS_FUNguild/`: Functional annotations of fungi using FUNGuild.
- `ITS_maaslin/`: Differential abundance testing using MaAsLin2.
- `ITS_reads_OG/`: Original raw reads and SampleSheet from sequencing runs.
- `ITS_EMP_paired_end_sequences/`: EMP formatted raw read directories for QIIME2 import.
- `ITS_taxonomy/`: Taxonomy assignment plots and summaries.
- `ITS_Plots/`: Finalized visualizations of key results.
- `ITS_slurm_scripts/` and `ITS_slurm_outputs/`: Scripts and logs for SLURM job execution.
- `ITS_NetworkAnalysis/`: In-progress or archived co-occurrence network analysis.

EOF

# === Create README.md for MAGs ===
mkdir -p MAGs
cat << 'EOF' > MAGs/README.md
# MAGs

This folder contains outputs from metagenome-assembled genomes (MAGs) processing and annotation.

## Subfolders & Files

- `outputs_versionDRAM2/`: Functional annotation results from DRAM v2.
- `outputs_versionDRAM1.4.4/`: Older DRAM outputs.
- `ASCC_16S_allbins_coremaslin.xlsx`: OTUs matched to MAGs for linking with MaAsLin2 results.
- `.svg`, `.ai`, `.pptx`: Concept figures related to carbon/nitrogen metabolism or MAG summaries.

EOF

# === Create README.md for Soil Chemistry ===
mkdir -p Soil_Chemistry_TimFegelRMRS
cat << 'EOF' > Soil_Chemistry_TimFegelRMRS/README.md
# Soil Chemistry (Tim Fegel RMRS)

Contains soil chemistry datasets and analyses provided by collaborators.

## Contents

- Water-extractable chemistry tables.
- Soil texture and bulk density measurements.
- PERMANOVA and envfit results.
- Associated RMarkdown analysis scripts and figures.

EOF

# === Create README.md for maaslin_files_probwrong ===
mkdir -p maaslin_files_probwrong
cat << 'EOF' > maaslin_files_probwrong/README.md
# maaslin_files_probwrong

Miscellaneous or problematic MaAsLin2 output files for debugging or validation.

This folder contains:
- Possibly mislabeled or duplicate MaAsLin2 result files
- Cross-checks of OTU presence/absence across result sheets

EOF

# === Create README.md for Manuscript folder ===
mkdir -p Manuscript_2025
cat << 'EOF' > Manuscript_2025/README.md
# Manuscript_2025

In-development figures, supplementary tables, and notes for the ASCC manuscript.

This includes:
- Final plots exported from R and Qiime2
- Supporting Excel files for reference
- Draft text or comments for figures

EOF

# Stage and commit all new README files
git add */README.md
git commit -m "Auto-generate subdirectory README files"
echo "README.md files created and staged for commit."
