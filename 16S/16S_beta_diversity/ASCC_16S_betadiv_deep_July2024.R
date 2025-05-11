###### DEPTH 5-15cm only
#*** July 2024 ***
#*** Redoing ASCC ITS beta diversity processing! ***
#*** be sure to edit metadata file to add location and depth ***
#*** i have removed blanks!***

# Load packages - install if not already installed!
install.packages("ggrepel")
install.packages("phyloseq")
BiocManager::install("phyloseq")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager",force=TRUE)
install.packages("grid")
install.packages("dplyr")
install.packages("reshape2")
install.packages("vegan")
install.packages("ggplot2")
install.packages("Hmisc")
install.packages("httpgd")

library(phyloseq)
library(ggplot2)
library(grid)
library(dplyr)
library(vegan)
library(Hmisc)
library(reshape2)
library(httpgd)
library(languageserver)

hgd()
hgd_view()
hgd_browse()

# Import your metadata
metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_metadata_rem.txt",row.names=1) 

# Import feature table
otus <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_feature_table_rem.txt",row.names=1,header=TRUE)

# check that the names match - they must be in the same order! If all works the result of this line will read 'TRUE' 
all.equal(names(otus),row.names(metadata)) #check to make sure row names and column names equal (aka the sample names), if this reads 'TRUE' you are good to move on!

otus_deep <- otus %>% 
  select(contains("_1C_15"))

meta_deep <- metadata %>% 
  filter(Depth == "5_15cm")

all.equal(names(otus_shallow),row.names(meta_shallow)) #check to make sure row names and column names equal (aka the sample names), if this reads 'TRUE' you are good to move on!

#### this is to check for missing samples between files.
## IGNORE UNLESS NEEDED
# sample_names_meta <- metadata %>% pull(Sample_ID)

# column_names_otus <- colnames(otus)

# missing_in_otus <- setdiff(sample_names_meta, column_names_otus)


# unexpected_in_otus <- setdiff(column_names_otus, sample_names_meta)

# list(
#   missing_in_otus = missing_in_otus,
#   unexpected_in_otus = unexpected_in_otus
# )
# C

# now we need to re-arrange the files, read in taxa, and make sure everything aligns before downstream steps, and make the phyloseq object
otumat<-as.matrix(otus_deep)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
taxa<-read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_taxonomy/ASCC_16S_taxonomy_forprocessing.txt",row.names=1) #this is a list of the ASV ids and the corresponding taxa strings
taxmat<-as.matrix(taxa)
all.equal(row.names(taxmat),row.names(otumat)) # again - check to make sure they are in the same order!
TAX = tax_table(taxmat)
physeq<-phyloseq(OTU,TAX) #make your phyloseq object, which is basically just a R object that has all of your data for the analyses stored into it
all.equal(row.names(meta_deep),sample_names(physeq)) # Final check to make sure everything lines up
sampledata<-sample_data(meta_deep)
mgd<-merge_phyloseq(physeq,sampledata) #make the final phyloseq object

#Calculate relative abundance from ASV count data - dont need to run this step if you have already done so but it shouldnt really effect! 
mgd_relabund<-transform_sample_counts(mgd,function(x)x/sum(x)) #Calculate relative abundance from ASV count data - this is what you can save and use for other taxonomy analyses. If you're interested in doing these, let me know and i can help out!

# Calculate Bray-Curtis dissimilary (distance) 
mgd_relabund.bray<-distance(mgd_relabund,"bray") 

# Final line to generate NMDS ordination 
mgd_relabund.bray.nmds<-ordinate(mgd_relabund,"NMDS",mgd_relabund.bray) #WAHOO DO NMDS!!

# Use this to read out the NMDS stress (you generally want the stress to be below 0.2)
mgd_relabund.bray.nmds$stress 
# stress = 0.1250581 July2024

#Create plain data frame 'map' of sample metadata
mgd_relabund_map=as(sample_data(mgd_relabund),"data.frame") 
sample_tab<-mgd_relabund_map

# Add NMDS coordinates (scores) to the sample table
sample_tab$NMDS1<-scores(mgd_relabund.bray.nmds$points)[,1] 
sample_tab$NMDS2<-scores(mgd_relabund.bray.nmds$points)[,2]

#Plot using ggplot! 
# Of course, you can change the variavles (what you color, size, etc by) as long as they were in your original metadata file

# This plot is points labeled by treatment and points labeled with names
# Need to remove outliers to see plot

hgd_view()

ggplot(sample_tab,aes(x=NMDS1,y=NMDS2,color=Site))+
  geom_point(size=2, show.legend = TRUE)+
  #scale_shape_manual(values = c(17, 16)) +
  theme(text=element_text(size = 24))+  # Change axis labels
  xlim(-1.4,1.7)+
  ylim(-1.51,1.67)+
  #geom_text_repel(max.overlaps = 20, aes(label=row.names(map_file))) + #use this line if you would like to have the samples labled in your ordination 
  theme_classic() 
  #scale_colour_manual(values = c("#ETRUE#scale_colour_manual(values = c("#EE2A7B", "#27AAE1","#CCCCCB","#2E3192", "#F15A29","#009444", "#AC8A13", "#00A79D", "#7F8080", "#9E1F63", "#877DBB")) # feel free to add custom colors if you'd like! or remove this line


ggplot(sample_tab,aes(x=NMDS1,y=NMDS2,color=Site,shape=Depth))+
  geom_point(size=2, show.legend = TRUE)+
  #scale_shape_manual(values = c(17, 16)) +
  theme(text=element_text(size = 24))+  # Change axis labels
  xlim(-1.4,1.7)+
  ylim(-1.51,1.67)+ 
  #geom_text_repel(max.overlaps = 20, aes(label=row.names(map_file))) + #use this line if you would like to have the samples labled in your ordination 
  theme_classic() 
  #scale_colour_manual(values = c("#ETRUE#scale_colour_manual(values = c("#EE2A7B", "#27AAE1","#CCCCCB","#2E3192", "#F15A29","#009444", "#AC8A13", "#00A79D", "#7F8080", "#9E1F63", "#877DBB")) # feel free to add custom colors if you'd like! or remove this line


mrpp(mgd_relabund.bray, sample_tab$Site, permutations=999, distance="bray")

anosim(mgd_relabund.bray, sample_tab$Site, permutations=999, distance="bray")

adonis2(mgd_relabund.bray ~ meta_deep$Site, dist="bray",perm=999)

pairwise.adonis2(mgd_relabund.bray ~ metadata$Site*metadata$Depth, data = metadata, sim.method = "bray", p.adjust.m = "BH", perm = 999)
#not working ?