** July 27th 2024-- correct script for ASCC ITS processing **

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

mkdir ITS

cd projects-wilkins-2/ASCC/summer_2024/ITS

ls

# check folder to be sure its there
# enter new directory

cd ITS

# Source into qiime2! Use version 2023.9!

source /opt/Miniconda2/miniconda2/bin/activate qiime2-2021.2

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
# go to ITS

ls

cd ..


# Make new directory for reads that are renamed! 
    ### make in ITS
# go back a directory to ITS (two directories back)

mkdir reads_renamed

cd ../../

# Enter reads_OG folder
# Copy files from reads_OG to new reads_renamed directory

cd reads_OG

# copy all files ending in .gz to new folder, go back a directory

cp *.gz /home/projects-wilkins-2/ASCC/summer_2024/1TS/reads_renamed

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

nano emp_import_ITS_ASCC.sh

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

# Exit script and save
    ### will sometimes need to deactivate qiime before running script
    ### go back a directory 

cd ..

# Make another directory for the EMP_paired_end_sequences and navigate to it

mkdir EMP_paired_end_sequences

cd EMP_paired_end_sequences

# Run script!
    ### check email to make sure it doesnt fail
    ### if it fails use tail slurm-<#>.out - dont use <>

sbatch /home/projects-wilkins-2/ASCC/summer_2024/ITS/slurm_scripts/emp_import_ITS_ASCC.sh

# check to see if sequences imported, then go back a directory 

ls

cd ..

# Demultiplex the sequences!
    ### create slurm script!

cd slurm_scripts

nano demux_ITS_ASCC.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00sque
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=128gb
#SBATCH --mail-user=kysparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low


source /opt/Miniconda2/miniconda2/bin/activate qiime2-2021.2

qiime demux emp-paired \
  --m-barcodes-file /home/projects-wilkins-2/ASCC/summer_2024/ITS/ITS_metadata/ASCC_ITS_barcodes.txt \
  --m-barcodes-column Golay-Barcode \
  --p-rev-comp-mapping-barcodes \
  --i-seqs /home/projects-wilkins-2/ASCC/summer_2024/ITS/EMP_paired_end_sequences/emp-paired-end-sequences.qza \ 
  --o-per-sample-sequences demux-ITS.qza \
  --o-error-correction-details demux-details-ITS.qza 


# Exit script and save
    ### will sometimes need to deactivate qiime before running script
    ### go back a directory

cd ..

# before we run the script lets make another directory for metadata

mkdir ITS_metadata

# move metadata file into new directory

mv *.csv metadata

# Make another directory for the demultiplex files and navigate to it

mkdir demultiplex_files

cd demultiplex_files

# Run script!
    ### check email to make sure it doesnt fail

sbatch /home/projects-wilkins-2/ASCC/summer_2024/ITS/slurm_scripts/demux_ITS_ASCC.sh

# Visualize sequence reads and sequence quality no need to run in slurm, download .qzv,
# go back a directory  

qiime demux summarize --i-data demux-ITS.qza --o-visualization demux-ITS.qzv

cd ..

cd slurm_scripts

nano DADA2_ITS_ASCC.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=200gb
#SBATCH --mail-user=kysparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-lo

source /opt/Miniconda2/miniconda2/bin/activate qiime2-2021.2

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /home/projects-wilkins-2/ASCC/summer_2024/ITS/demultiplex_files/demux-ITS.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 250 \
  --o-table table-ITS.qza \
  --o-representative-sequences rep-seqs-ITS.qza \
  --o-denoising-stats denoising-stats-ITS.qza

# Exit

cd ..

cd denoise_files

# run script

sbatch /home/projects-wilkins-2/ASCC/summer_2024/ITS/slurm_scripts/DADA2_ITS_ASCC.sh

qiime metadata tabulate \
--m-input-file denoising-stats-ITS.qza \
--o-visualization denoising-stats-ITS.qzv

# Visualize feature table!

qiime feature-table summarize \
 --i-table table-ITS.qza  \
 --o-visualization table-ITS.qzv \
 --m-sample-metadata-file /home/projects-wilkins-2/ASCC/summer_2024/ITS/ITS_metadata/ASCC_ITS_barcodes.txt

# Visualize representative sequences!

qiime feature-table tabulate-seqs \
--i-data rep-seqs-ITS.qza \
--o-visualization rep-seqs-ITS.qzv

cd ...

cd slurm_scripts

nano classifier.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=128gb
#SBATCH --mail-user=kya.sparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low

source /opt/Miniconda2/miniconda2/bin/activate qiime2-2021.2

qiime feature-classifier classify-sklearn \
--i-classifier /home/Database/qiime2/classifiers/qiime2-2021.2/classifier_unite_ver9.qza\
--i-reads /home/projects-wilkins-2/ASCC/summer_2024/ITS/denoise_files/rep-seqs-ITS.qza \
--o-classification taxonomy_unite_ITS.qza

# Exit!


cd ..

mkdir taxonomy 

cd taxonomomy

sbatch /home/projects-wilkins-2/ASCC/summer_2024/ITS/slurm_scripts/classifier.sh

qiime metadata tabulate \
--m-input-file taxonomy_unite_ITS.qza \
--o-visualization taxonomy_unite_ITS.qzv

unzip taxonomy_unite_ITS.qza

biom convert -i feature-table.biom -o table.from_biom.txt --to-tsv

qiime taxa barplot \
--i-table /home/projects-wilkins-2/ASCC/summer_2024/ITS/denoise_files/table-ITS.qza \
--i-taxonomy taxonomy_unite_ITS.qza \
--m-metadata-file /home/projects-wilkins-2/ASCC/summer_2024/ITS/ITS_metadata/ASCC_ITS_barcodes.txt \
--o-visualization barplots-ITS.qzv

qiime diversity alpha-rarefaction \
--i-table table-ITS.qza \
--m-metadata-file /home/projects-wilkins-2/ASCC/summer_2024/ITS/ITS_metadata/ASCC_ITS_barcodes.txt \
--o-visualization alpha_rarefaction_curves.qzv \
--p-min-depth 10 \
--p-max-depth 15000

** TREE **

mkdir tree

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=128gb
#SBATCH --mail-user=kysparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low

source /opt/Miniconda2/miniconda2/bin/activate qiime2-2021.2

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /home/projects-wilkins-2/ASCC/summer_2024/ITS/denoise_files/rep-seqs-ITS.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

cd tree

sbatch /home/projects-wilkins-2/ASCC/summer_2024/ITS/slurm_scripts/tree.sh


** RARIFY **

mkdir alpha_rarefaction

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
  --i-table /home/projects-wilkins-2/ASCC/summer_2024/ITS/taxonomy/taxonomy_unite_ITS.qza \
  --i-phylogeny /home/projects-wilkins-2/ASCC/summer_2024/ITS/tree/rooted-tree.qza \
  --p-max-depth 6000 \
  --m-metadata-file /home/projects-wilkins-2/ASCC/summer_2024/ITS/ITS_metadata/ASCC_ITS_barcodes3.txt \
  --o-visualization alpha-rarefaction.qzv

cd alpha_rarefaction

sbatch /home/projects-wilkins-2/ASCC/summer_2024/ITS/slurm_scripts/rarefaction.sh

nano coremetrics.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=128gb
#SBATCH --mail-user=kysparks@colostate.edu
#SBATCH --partition=wilkins-hi,wilkins-low

source /opt/Miniconda2/miniconda2/bin/activate qiime2-2021.2

 qiime diversity core-metrics-phylogenetic \
  --i-phylogeny /home/projects-wilkins-2/ASCC/summer_2024/ITS/tree/rooted-tree.qza \
  --i-table /home/projects-wilkins-2/ASCC/summer_2024/ITS/denoise_files/table-ITS.qza \
  --p-sampling-depth 6000 \
  --m-metadata-file /home/projects-wilkins-2/ASCC/summer_2024/ITS/ITS_metadata/ASCC_ITS_barcodes3.txt \
  --output-dir core-metrics-results-6000

cd coremets

sbatch /home/projects-wilkins-2/ASCC/summer_2024/ITS/slurm_scripts/coremetrics.sh

unzip core-metrics-results-6000/rarefied_table.qza

biom convert -i feature-table.biom -o table.from_biom.txt --to-tsv
