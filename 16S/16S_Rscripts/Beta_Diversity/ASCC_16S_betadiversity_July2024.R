# *** July 2024 ***
# *** ALL DEPTHS ****
# *** REMOVED MITOCHONDRIA AND CHLOROPLAST ***
# *** Redoing ASCC 16S beta diversity processing! ***
# *** be sure to edit metadata file to add location and depth ***
# *** i have removed blanks!***
# *** code for rarefied at bottom ***

# Load packages - install if not already installed!
# install.packages("phyloseq")
# # BiocManager::install("phyloseq")
# # if (!require("BiocManager", quietly = TRUE))
# # install.packages("BiocManager",force=TRUE)
# install.packages("reshape2")
# install.packages("vegan")
# install.packages("ggplot2")
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis"

#install.packages("paletteer")
library(phyloseq)  # Includes microbiome data analysis tools, handles OTU tables
library(vegan)     # For NMDS and other ordination methods
library(ggplot2)   # For customizable data visualization
library(dplyr)     # For data manipulation (tidyverse core)
library(paletteer) # For custom color palettes (optional, if using specific palettes)
library(pairwiseadonis)

# #### this is to check for missing samples between files.
# ## IGNORE UNLESS NEEDED
# sample_names_meta <- metadata %>% pull(Sample_ID)

# column_names_otus <- colnames(otus)

# missing_in_otus <- setdiff(sample_names_meta, column_names_otus)


# unexpected_in_otus <- setdiff(column_names_otus, sample_names_meta)

# list(
#   missing_in_otus = missing_in_otus,
#   unexpected_in_otus = unexpected_in_otus
# )

# Import your metadata

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_metadata_rem.txt",row.names=1) 

# Import feature table
otus <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_feature_table_rem.txt", row.names=1,header=TRUE)

# check that the names match - they must be in the same order! If all works the result of this line will read 'TRUE' 
all.equal(names(otus),row.names(metadata)) #check to make sure row names and column names equal (aka the sample names), if this reads 'TRUE' you are good to move on!

# now we need to re-arrange the files, read in taxa, and make sure everything aligns before downstream steps, and make the phyloseq object
otumat<-as.matrix(otus)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
taxa<-read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_taxonomy/filtered/ASCC_16S_taxonomy_forprocessing_copy.txt",row.names=1) #this is a list of the ASV ids and the corresponding taxa strings
taxmat<-as.matrix(taxa)
all.equal(row.names(taxmat),row.names(otumat)) # again - check to make sure they are in the same order!
TAX = tax_table(taxmat)
physeq<-phyloseq(OTU,TAX) #make your phyloseq object, which is basically just a R object that has all of your data for the analyses stored into it
all.equal(row.names(metadata),sample_names(physeq)) # Final check to make sure everything lines up
sampledata<-sample_data(metadata)
mgd<-merge_phyloseq(physeq,sampledata) #make the final phyloseq object

#Calculate relative abundance from ASV count data - dont need to run this step if you have already done so but it shouldnt really effect! 
mgd_relabund<-transform_sample_counts(mgd,function(x)x/sum(x)) #Calculate relative abundance from ASV count data - this is what you can save and use for other taxonomy analyses. If you're interested in doing these, let me know and i can help out!

# Calculate Bray-Curtis dissimilary (distance) 
mgd_relabund.bray<-distance(mgd_relabund,"bray") 

# Final line to generate NMDS ordination 
mgd_relabund.bray.nmds<-ordinate(mgd_relabund,"NMDS",mgd_relabund.bray) #WAHOO DO NMDS!!

# Use this to read out the NMDS stress (you generally want the stress to be below 0.2)
mgd_relabund.bray.nmds$stress 
# stress =  0.1441404 july2024

#Create plain data frame 'map' of sample metadata
mgd_relabund_map=as(sample_data(mgd_relabund),"data.frame") 
sample_tab<-mgd_relabund_map

# Add NMDS coordinates (scores) to the sample table
sample_tab$NMDS1<-scores(mgd_relabund.bray.nmds$points)[,1] 
sample_tab$NMDS2<-scores(mgd_relabund.bray.nmds$points)[,2]

#Plot using ggplot! 
# Of course, you can change the variables (what you color, size, etc by) as long as they were in your original metadata file

# This plot is points labeled by treatment and points labeled with names
# Need to remove outliers to see plot
# reordering legends, refactoring
sample_tab$Depth <- factor(sample_tab$Depth, levels=c('OM','0_5cm','5_15cm'))
sample_tab$Site <- factor(sample_tab$Site, levels=c('SF','TP','SJ'))


ggplot(sample_tab,aes(x=NMDS1,y=NMDS2,color=Site, shape=Depth))+
  geom_point(size=2, show.legend = TRUE)+
  #scale_shape_manual(values = c(17, 16)) +
  theme(text=element_text(size = 24))+ #change legend text font size)+  # Change axis labels
  labs(title = "NMDS of Microbial Community Composition Across Sites and Depths")+
  scale_colour_paletteer_d("lisa::BridgetRiley")+
  scale_fill_paletteer_d("lisa::BridgetRiley")+
  xlim(-1.4,1.7)+
  ylim(-1.51,1.67)+
  #geom_text_repel(max.overlaps = 20, aes(label=row.names(map_file))) + #use this line if you would like to have the samples labled in your ordination 
  theme_classic() 

## RUN STATS

mrpp(mgd_relabund.bray, sample_tab$Site, permutations=999, distance="bray")

anosim(mgd_relabund.bray, sample_tab2$Site, permutations=999, distance="bray")

# PERMANOVA

adonis2(mgd_relabund.bray ~ metadata$Site*metadata$Depth, dist="bray",perm=999)
adonis2(mgd_relabund.bray ~ sample_tab2$Site*sample_tab$Depth, dist="bray",perm=999)

pairwise.adonis2(mgd_relabund.bray ~ metadata$Site*metadata$Depth, data = metadata, sim.method = "bray", p.adjust.m = "BH", perm = 999)
#not working ?

write.csv(sample_tab,"/Users/kyasparks/Desktop/Desktop - MacBook Pro (2)/Kyas_PhD/CSU_2023_2024/Courses/ESS523A/Final-Project-Explore/shinywebsite/www/sampletab.csv")







