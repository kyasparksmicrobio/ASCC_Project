#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=128gb
#SBATCH --mail-user=kysparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low

source /opt/Miniconda2/miniconda2/bin/activate qiime2-2021.2

qiime diversity alpha-rarefaction \
  --i-table /home/projects-wilkins-2/ASCC/summer_2024/ITS/taxonomy/taxonomy_unite_ITS.qza 
  --i-phylogeny /home/projects-wilkins-2/ASCC/summer_2024/ITS/tree/rooted-tree.qza \
  --p-max-depth 6000 \
  --m-metadata-file /home/projects-wilkins-2/ASCC/summer_2024/ITS/ITS_metadata/ASCC_ITS_barcodes3.txt \
  --o-visualization alpha-rarefaction.qzv
