#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --mem=170gb
#SBATCH --time=28-00:00:00
#SBATCH --job-name=dram1.4.4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kysparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low
#SBATCH --nodelist=marmot

source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4

DRAM.py annotate -i '/home/projects-wilkins-2/ASCC/summer_2024/link_16S_to_MAGs/MAGs/*.fa' -o /home/projects-wilkins-2/ASCC/summer_2024/link_16S_to_MAGs/MAGs/dram1.4.4_ANNOTATIONS --threads 15 --use_fegenie --use_sulfur

DRAM.py distill -i /home/projects-wilkins-2/ASCC/summer_2024/link_16S_to_MAGs/MAGs/dram1.4.4_ANNOTATIONS/annotations.tsv -o /home/projects-wilkins-2/ASCC/summer_2024/link_16S_to_MAGs/MAGs/dram1.4.4_ANNOTATIONS/summarize --rrna_path /home/projects-wilkins-2/ASCC/summer_2024/link_16S_to_MAGs/MAGs/dram1.4.4_ANNOTATIONS/rrnas.tsv --trna_path /home/projects-wilkins-2/ASCC/summer_2024/link_16S_to_MAGs/MAGs/dram1.4.4_ANNOTATIONS/trnas.tsv


# 1. DRAM.py : you can output separate vs all together, single quotes because wildcard with python