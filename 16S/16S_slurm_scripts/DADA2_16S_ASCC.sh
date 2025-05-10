#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=20gb
#SBATCH --mail-user=kysparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low

source /home/opt/Miniconda3/miniconda3/bin/activate qiime2-2023.9

qiime dada2 denoise-paired \
--i-demultiplexed-seqs /home/projects-wilkins-2/ASCC/summer_2024/16S/demultiplex_files/demux-16S.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 250 \
--p-trunc-len-r 250 \
--o-table table-16S.qza \
--o-representative-sequences rep-seqs-16S.qza \
--o-denoising-stats denoising-stats-16S.qza
