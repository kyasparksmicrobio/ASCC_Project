** July 27th 2024-- correct script for ASCC 16S processing **

# Use zsh / shell script /UTF-8

# working directory as soon as you enter the terminal
# /Users/kyasparks

# log in to the Wilkins lab server

ssh ksparks@zenith.aggie.colostate.edu

# will prompt you for a password 

password = in4kyams

# make sure ASCC folder is cleaned up and organized before starting processing
# Now naviagte to the project folder (ASCC)
# how I get to the correct directory to find the ASCC folders:

cd ..

pwd

# this brings me to /home
# Navigate to ASCC folder
# since we are starting from scratch, make sure we make a new directory

cd projects-wilkins-2/ASCC/summer_2024

mkdir 16S

cd projects-wilkins-2/ASCC/summer_2024/16S

ls

# check folder to be sure its there
# enter new directory

cd 16S

# Source into qiime2! Use version 2023.9!

source /home/opt/Miniconda3/miniconda3/bin/activate qiime2-2023.9

### The above gives you the latest version of qiime2. Sometimes you can't run files from 
    ### a previous version of qiime2 in a new version. If you want to work with older files 
    ### you can activate the previous version that is still on the server (below).
    
    ### source /opt/Miniconda2/miniconda2/bin/activate qiime2-2021.2


# Import sequencing files into qiime2
# make new directory and navigate to that directory

mkdir reads_OG

cd reads_OG

# use filezilla to copy original reads into new directory!!
# when finished, check to make sure files copied, then go back a directory
# go to 16S

ls

cd ..


# Make new directory for reads that are renamed! 
    ### make in 16S
# go back a directory to 16S (two directories back)

mkdir reads_renamed

cd ../../

# Enter reads_OG folder
# Copy files from reads_OG to new reads_renamed directory

cd reads_OG

# copy all files ending in .gz to new folder, go back a directory

cp *.gz /home/projects-wilkins-2/ASCC/summer_2024/16S/reads_renamed

cd ../..

# navigate to the new reads_renamed directory

cd reads_renamed

# Rename forward, reverse, and barcode fastq.gz files, then navigate back a directory

mv Undetermined_S0_L001_I1_001.fastq.gz barcodes.fastq.gz
mv Undetermined_S0_L001_R1_001.fastq.gz forward.fastq.gz
mv Undetermined_S0_L001_R2_001.fastq.gz reverse.fastq.gz

cd..

# Make new directory for slurm scripts and slurm outputs
# make sure to move output to slurm_output folders after finished running

mkdir slurm_scripts

mkdir slurm_outputs

# Import Sequence files using slurm 
    ### advice is to use slurm for basically everything 
# Create slurm script!!

cd slurm_scripts

nano emp_import_16S_ASCC.sh

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
--output-path emp-paired-end-sequences.qza

# Exit script and save
    ### will sometimes need to deactivate qiime before running script
    ### go back a directory to 16S_amplicon_analyses_May2024_KMS

cd ..

# Make another directory for the EMP_paired_end_sequences and navigate to it

mkdir EMP_paired_end_sequences

cd EMP_paired_end_sequences

# Run script!
    ### check email to make sure it doesnt fail
    ### if it fails use tail slurm-<#>.out - dont use <>

sbatch /home/projects-wilkins-2/ASCC/summer_2024/16S/slurm_scripts/emp_import_16S_ASCC.sh

# check to see if sequences imported, then go back a directory 

ls

cd ..

# Demultiplex the sequences!
    ### create slurm script!

cd slurm_scripts

nano demux_16S_ASCC.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=128gb
#SBATCH --mail-user=your.email@colostate.edu
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


# Exit script and save
    ### will sometimes need to deactivate qiime before running script
    ### go back a directory to 16S_amplicon_analyses_May2024_KMS

cd ..

# before we run the script lets make another directory for metadata

mkdir 16S_metadata

# move metadata file into new directory

mv *.csv 16S_metadata

# Make another directory for the demultiplex files and navigate to it

mkdir demultiplex_files

cd demultiplex_files

# Run script!
    ### check email to make sure it doesnt fail

sbatch /home/projects-wilkins-2/ASCC/summer_2024/16S/slurm_scripts/demux_16S_ASCC.sh

# Visualize sequence reads and sequence quality no need to run in slurm, download .qzv,
# go back a directory  

qiime demux summarize --i-data demux-16S.qza --o-visualization demux-16S.qzv

cd ..

# Run DADA2 to denoise and merge sequence reads 
# Create slurm script!!

cd slurm_scripts

nano DADA2_16S_ASCC.sh

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

# Exit script and save
    ### will sometimes need to deactivate qiime before running script
    ### go back a directory to 16S_amplicon_analyses_May2024_KMS

cd ..

# Make another directory for the EMP_paired_end_sequences and navigate to it

mkdir denoise_files

cd denoise_files

# Run script!
    ### check email to make sure it doesnt fail
    ### if it fails use tail slurm-<#>.out - dont use <>

sbatch /home/projects-wilkins-2/ASCC/summer_2024/16S/slurm_scripts/DADA2_16S_ASCC.sh

# check to see if sequences imported, then go back a directory 

ls

cd ..


qiime metadata tabulate \
--m-input-file denoising-stats-16S.qza \
--o-visualization denoising-stats-16S.qzv

# Visualize feature table!

qiime feature-table summarize \
 --i-table table-16S.qza  \
 --o-visualization table-16S.qzv \
 --m-sample-metadata-file /home/projects-wilkins-2/ASCC/summer_2024/16S/16S_metadata/ASCC_16S_barcodes.txt

# Visualize representative sequences!

qiime feature-table tabulate-seqs \
--i-data rep-seqs-16S.qza \
--o-visualization rep-seqs-16S.qzv

cd ...

cd slurm_scripts

nano classifier.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=128gb
#SBATCH --mail-user=kysparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low

source /home/opt/Miniconda3/miniconda3/bin/activate qiime2-2023.9

qiime feature-classifier classify-sklearn \
--i-classifier /home/Database/qiime2/classifiers/qiime2-2023.9/silva-138-99-515-806-nb-classifier.qza \
--i-reads /home/projects-wilkins-2/ASCC/summer_2024/16S/denoise_files/rep-seqs-16S.qza \
--o-classification taxonomy_silva138_16S.qza

# Exit!

cd ..

mkdir taxonomy 

cd taxonomomy

sbatch /home/projects-wilkins-2/ASCC/summer_2024/16S/slurm_scripts/classifier.sh

qiime metadata tabulate \
--m-input-file taxonomy_silva138_16S.qza \
--o-visualization taxonomy_silva138_16S.qzv

qiime taxa filter-table \
--i-table /home/projects-wilkins-2/ASCC/summer_2024/16S/denoise_files/table-16S.qza  \
--i-taxonomy /home/projects-wilkins-2/ASCC/summer_2024/16S/taxonomy/taxonomy_silva138_16S.qza \
--p-exclude mitochondria,chloroplast,Unassigned,Eukaryota \
--o-filtered-table ttaxonomy_silva138_16S_no_chloro.qza

unzip taxonomy_silva138_16S_no_chloro.qza

biom convert -i feature-table.biom -o table.from_biom.txt --to-tsv

qiime taxa barplot \
--i-table /home/projects-wilkins-2/ASCC/summer_2024/16S/denoise_files/table-16S.qza \
--i-taxonomy taxonomy_silva138_16S.qza \
--m-metadata-file /home/projects-wilkins-2/ASCC/summer_2024/16S/16S_metadata/ASCC_16S_metadata.txt \
--o-visualization barplots-16S.qzv

