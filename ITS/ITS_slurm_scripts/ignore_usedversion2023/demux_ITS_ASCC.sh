#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=128gb
#SBATCH --mail-user=kysparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low

source /home/opt/Miniconda3/miniconda3/bin/activate qiime2-2023.9

qiime demux emp-paired \
--m-barcodes-file /home/projects-wilkins-2/ASCC/summer_2024/ITS/ITS_metadata/ASCC_ITS_barcodes.txt \
--m-barcodes-column Golay-Barcode \
--p-rev-comp-mapping-barcodes \
--i-seqs /home/projects-wilkins-2/ASCC/summer_2024/ITS/EMP_paired_end_sequences/emp-paired-end-sequences.qza \
--o-per-sample-sequences demux-ITS.qza \
--o-error-correction-details demux-details-ITS.qza 
