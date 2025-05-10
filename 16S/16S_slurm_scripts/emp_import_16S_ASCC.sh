#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=128gb
#SBATCH --mail-user=kysparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low

source /home/opt/Miniconda3/miniconda3/bin/activate qiime2-2023.9

qiime tools import \
--type EMPPairedEndSequences \
--input-path /home/projects-wilkins-2/ASCC/summer_2024/16S/reads_renamed \
--output-path /home/projects-wilkins-2/ASCC/summer_2024/16S/EMP_paired_end_sequences/emp-paired-end-sequences.qza
