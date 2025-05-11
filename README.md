# ASCC\_Project: Adaptive Silviculture for Climate Change (ASCC) Soil Microbiome Analyses

This repository contains all scripts, metadata, outputs, and visualizations for 16S and ITS rRNA gene sequencing analyses associated with the **Adaptive Silviculture for Climate Change (ASCC)** soil microbial project. The project spans coniferous forest sites in Colorado and involves bioinformatic and statistical workflows for understanding microbial community structure and function in response to climate-adaptive forest management.

---

## Repository Structure

### `16S/`

Contains all analyses related to bacterial and archaeal 16S rRNA gene sequencing.

Subdirectories:

* `16S_alpha_diversity/`: Alpha diversity metrics, rarefaction, and plots.
* `16S_beta_diversity/`: NMDS, PCoA, PERMANOVA, and environmental vector fitting.
* `16S_coremetrics/`: Qiime2-generated metrics.
* `16S_CoreTaxaAnalysis/`: Core taxa identification and ternary plots.
* `16S_MAGs_coniferous/`: Mapping ASVs to MAGs and associated results.
* `16S_Maaslin/`: MaAsLin2 differential abundance analyses, heatmaps, point-range plots.
* `16S_metadata/`: Cleaned and formatted metadata, barcode files.
* `16S_Rscripts/`: R scripts used for analysis and plotting.
* `16S_slurm_outputs/`: Cluster output logs.
* `16S_slurm_scripts/`: All bash scripts used for SLURM-based Qiime2 pipeline.
* `16S_taxonomy/`: Taxonomic classification and barcharts.
* `16S_barcodes/`: Barcode mapping files.

### `ITS/`

Contains analyses of fungal ITS rRNA gene sequencing data.

Subdirectories:

* `ITS_alpha_diversity/`: Rarefaction, diversity comparisons.
* `ITS_beta_diversity/`: NMDS and PCoA with depth/site comparisons.
* `ITS_coremetrics/`: Core metrics and EM fungal detection.
* `ITS_CoreTaxaAnalysis/`: Dominant/consistent taxa analysis.
* `ITS_FUNguild/`: Functional guild parsing and results.
* `ITS_maaslin/`: Differential abundance outputs.
* `ITS_metadata/`: ITS-specific metadata files.
* `ITS_reads_OG/`: Raw read sample sheet.
* `ITS_Plots/`: Barcharts, volcano plots, and exploratory graphics.
* `ITS_Rscripts/`: R scripts for diversity, guild, and plotting.
* `ITS_slurm_outputs/`: SLURM output files.
* `ITS_slurm_scripts/`: ITS pipeline bash scripts for Qiime2.
* `ITS_taxonomy/`: Taxonomy classification and barcharts.

### `MAGs/`

Contains metagenome-assembled genome (MAG) analysis outputs.

* `outputs_versionDRAM2/` and `outputs_versionDRAM1.4.4/`: Annotated metabolism summaries, genome stats, product tables.
* All `.tsv`, `.csv`, `.html`, and `.xlsx` files from DRAM are tracked.

### `Soil_Chemistry_TimFegelRMRS/`

Processed soil chemistry data provided by Tim Fegel (USFS RMRS). Includes:

* Water extractable chemistry
* Mineral/organic soil fractions
* Soil texture and bulk density
* PERMANOVA results and ordinations

### `Manuscript_2025/`

Draft text and figures supporting the ASCC 2025 manuscript.

### `maaslin_files_probwrong/`

A scratch folder for potentially incorrect or early versions of MaAsLin2 outputs.

---

## Notebooks and Reports

* `.Rmd`, `.R`, and `.md` files throughout the repo document analytical decisions and steps.
* All heatmaps, stacked barplots, and point-range visualizations for MaAsLin2 are included in `geom_point_range` and `stacked_bar_maaslin` folders.

---

## Git LFS Usage

Some files exceed GitHub's default 100 MB limit. [Git LFS](https://git-lfs.github.com/) is used to track large:

* `.csv` output tables from ITS-FUNGuild
* DRAM2 `.tsv` annotation tables
* Feature tables for taxonomy and abundance

Make sure Git LFS is installed before cloning:

```bash
brew install git-lfs
# or
sudo apt install git-lfs

# Then
git lfs install
```

---

## .gitignore Summary

This repository tracks only the necessary text files and figures. Images, PowerPoint slides, and intermediate outputs are ignored globally **except** where explicitly preserved in `MAGs/` and `Maaslin/` folders.

---

## Citation

If you use any part of this repository in your work, please cite:

> Sparks, K. et al. (2025). Soil microbial responses to adaptive silviculture for climate change in western coniferous forests. *In prep.*

---

## Contact

Kya Sparks — Soil Microbiology PhD Student, Colorado State University
[kya.sparks@colostate.edu](mailto:kya.sparks@colostate.edu)
GitHub: [kyasparksmicrobio](https://github.com/kyasparksmicrobio)

---

Let me know if you'd like this broken into sub-README files by folder.

