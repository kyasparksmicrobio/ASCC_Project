#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=128gb
#SBATCH --mail-user=kysparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low

source /home/opt/Miniconda3/miniconda3/bin/activate qiime2-2021.2

qiime feature-classifier classify-sklearn \
--i-classifier /home/Database/qiime2/classifiers/qiime2-2021.2/classifier_unite_ver9.qza \
--i-reads /home/projects-wilkins-2/ASCC/summer_2024/ITS/denoise_files/rep-seqs-ITS.qza \
--o-classification taxonomy_unite_ITS.qza
