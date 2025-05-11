# ASCC Project: Soil Microbiome Analysis

This repository contains data processing, analysis scripts, and results for the **Adaptive Silviculture for Climate Change (ASCC)** soil microbiome study. The project includes 16S rRNA and ITS amplicon data, soil chemical metadata, and code used for bioinformatics and statistical analysis.

---

## 📁 Project Structure

```
ASCC_Project/
├── 16S/                     # 16S rRNA analysis
│   ├── 16S_Rscripts/        # Alpha, beta diversity, MaAsLin2, taxonomy
│   ├── 16S_metadata/        # Feature tables and metadata
│   └── 16S_slurm_scripts/   # SLURM scripts for cluster jobs
│
├── ITS/                     # ITS fungal analysis
│   ├── ITS_Rscripts/        # Scripts for diversity, guilds, core taxa
│   ├── ITS_metadata/        # ITS-specific metadata and feature tables
│   └── ITS_slurm_scripts/   # SLURM scripts for ITS
│
├── MAGs/                    # DRAM outputs and genome annotations
├── Soil_Chemistry_TimFegelRMRS/   # Soil chemical data and analyses
├── maaslin_files_probwrong/       # Debugging/test MaAsLin2 files
├── Manuscript_2025/        # Manuscript drafts and supporting figures
├── .gitattributes          # Git LFS configuration
├── .gitignore              # Ignore rules for untracked files
└── README.md               # Project overview
```

---

## 🔬 Summary

This project examines bacterial and fungal community responses across different depths and sites in forest soils under climate-adapted silvicultural treatments. It includes:

- **16S**: bacterial/archaeal community profiling  
- **ITS**: fungal community profiling  
- **MAGs**: metagenome-assembled genomes & DRAM annotation  
- **Soil chemistry**: pH, nutrients, and other abiotic drivers

---

## 💻 Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/kyasparksmicrobio/ASCC_Project.git
cd ASCC_Project
```

### 2. Install Git LFS (for large file support)

```bash
brew install git-lfs   # macOS
git lfs install
```

---

## 📊 Key Tools and Pipelines

| Component     | Tools Used                              |
|---------------|------------------------------------------|
| Amplicon QC   | Qiime2, DADA2                            |
| Diversity     | `vegan`, `phyloseq`, custom `ggplot2`    |
| Differential Abundance | MaAsLin2                        |
| Functional Analysis | DRAM (MAGs)                        |
| Guild Analysis | FUNGuild, custom R scripts              |

---

## 📁 Data Notes

- Large files (e.g., `.csv`, `.xlsx`, `.tsv`) are tracked with **Git LFS**.
- Intermediate `.qza`, `.fastq`, and `.pptx` files are excluded via `.gitignore`.
- Sensitive raw sequencing data is stored on [external drive/cluster path].

---

## 🧬 Git LFS Warnings

GitHub may warn about large files (>50MB), but these are safely managed with Git LFS. If you clone this repo, make sure you run:

```bash
git lfs install
git lfs pull
```

---

## ✍️ Contact

Maintained by **Kya Sparks** (kyasparks@colostate.edu)  
Colorado State University, Soil & Crop Sciences / Forest Microbiomes Lab

---

## 📄 License

[MIT License](LICENSE) — open for academic use and collaboration.
