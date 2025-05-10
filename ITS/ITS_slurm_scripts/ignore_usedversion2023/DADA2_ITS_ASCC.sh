#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=200gb
#SBATCH --mail-user=kya.sparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low

source /home/opt/Miniconda3/miniconda3/bin/activate qiime2-2021.2


qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /home/projects-wilkins-2/ASCC/summer_2024/ITS/demultiplex_files/demux-ITS.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 250 \
  --o-table table-ITS.qza \
  --o-representative-sequences rep-seqs-ITS.qza \
  --o-denoising-stats denoising-stats-ITS.qza
