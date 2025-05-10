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
  --m-barcodes-file  /home/projects-wilkins-2/ASCC/summer_2024/16S/16S_metadata/ASCC_16S_barcodes.txt  \
  --m-barcodes-column Golay-Barcode \
  --p-no-rev-comp-mapping-barcodes \
  --p-no-golay-error-correction\
  --i-seqs /home/projects-wilkins-2/ASCC/summer_2024/16S/EMP_paired_end_sequences/emp-paired-end-sequences.qza \
  --o-per-sample-sequences demux-16S.qza \
  --o-error-correction-details demux-16S-details.qza 

