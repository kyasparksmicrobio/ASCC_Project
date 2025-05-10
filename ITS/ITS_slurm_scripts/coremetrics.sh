#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=128gb
#SBATCH --mail-user=kysparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low

source /opt/Miniconda2/miniconda2/bin/activate qiime2-2021.2

 qiime diversity core-metrics-phylogenetic \
  --i-phylogeny /home/projects-wilkins-2/ASCC/summer_2024/ITS/tree/rooted-tree.qza \
  --i-table /home/projects-wilkins-2/ASCC/summer_2024/ITS/denoise_files/table-ITS.qza \
  --p-sampling-depth 6000 \
  --m-metadata-file /home/projects-wilkins-2/ASCC/summer_2024/ITS/ITS_metadata/ASCC_ITS_barcodes3.txt \
  --output-dir core-metrics-results-6000
