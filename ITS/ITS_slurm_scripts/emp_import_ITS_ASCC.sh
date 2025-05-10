#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=128gb
#SBATCH --mail-user=kysparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low

source /opt/Miniconda2/miniconda2/bin/activate qiime2-2021.2

qiime tools import \
--type EMPPairedEndSequences \
--input-path /home/projects-wilkins-2/ASCC/summer_2024/ITS/reads_renamed \
--output-path emp-paired-end-sequences.qza
