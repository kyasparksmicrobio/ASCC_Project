# *** July 2024 ***
# *** BOTH DEPTHS ****
# *** REMOVED MITOCHONDRIA AND CHLOROPLAST ***
# *** Redoing ASCC ITS beta diversity processing! ***
# *** be sure to edit metadata file to add location and depth ***
# *** i have removed blanks!***
# *** code for rarefied at bottom ***
# 
# ** removiing ASCC 77,112 (KEEPING REDO),122 (KEEP REDO),303,73 redo,217,317 (keeping redo),299 (keeping redo),184,327 -- 0 reads, part of ~10% fail **
# ** REMOVING ascc 128 (keeping redo),ascc 101 (keeping redo), 114 (keeping redo), 123 redo, 128 (keeping redo), 13 (keeping redo), 110 dupe, 154 redo, 
#     171 redo, 194 (keeping redo), 220 dupe, 194 dupe2, 239 redo, 247 (keeping redo), 258 dupe, 282 (keeping redo), 129 redo, 292 dupe, 299 (keeping redo),
#     309 redo, 31 (keeping redo), 317 (keeping redo), 326 redo, 332 dupe, 333 dupe, 42 redo, 16 dupe, 50 (keeping redo), 6 dupe,248 REDO, 266 DUPE  **
# ** changed the names to no longer to have "dupe" or "redo"

# Load packages - install if not already installed!
# install.packages("ggrepel")
# install.packages("phyloseq")
# BiocManager::install("phyloseq")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager",force=TRUE)
# install.packages("grid")
# install.packages("dplyr")
# install.packages("reshape2")
# install.packages("vegan")
# install.packages("ggplot2")
# install.packages("Hmisc")
# install.packages("httpgd")

library(phyloseq)
library(ggplot2)
library(grid)
library(dplyr)
library(vegan)
library(Hmisc)
library(reshape2)
# library(httpgd)
# library(languageserver)
# 
# hgd()
# hgd_view()
# hgd_browse()

#### this is to check for missing samples between files.
## IGNORE UNLESS NEEDED
sample_names_meta <- metadata %>% 
pull(samples)

column_names_otus <- colnames(otus)
missing_in_otus <- setdiff(sample_names_meta, column_names_otus)
unexpected_in_otus <- setdiff(column_names_otus, sample_names_meta)
list(
  missing_in_otus = missing_in_otus,
  unexpected_in_otus = unexpected_in_otus
)
#Import your metadata

metadata <- read.delim("/Users/kyasparks/Library/Mobile Documents/com~apple~CloudDocs/Desktop - MacBook Pro (2)/Kyas_PhD/CSU_2023_2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_metadata/ASCC_ITS_metadata_rem.txt",row.names=1) 

# Import feature table
otus <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2028/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG_rem.txt",row.names=1,header=TRUE)

# check that the names match - they must be in the same order! If all works the result of this line will read 'TRUE' 
all.equal(names(otus),row.names(metadata)) #check to make sure row names and column names equal (aka the sample names), if this reads 'TRUE' you are good to move on!

# now we need to re-arrange the files, read in taxa, and make sure everything aligns before downstream steps, and make the phyloseq object
otumat<-as.matrix(otus)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
taxa<-read.delim("/Users/kyasparks/Desktop/CSU 2023-2028/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_taxonomy/taxonomy_ITS_ASCC.txt",row.names=1) #this is a list of the ASV ids and the corresponding taxa strings
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
# stress =  0.2406392 August2024

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
# hgd_view()

ggplot(sample_tab,aes(x=NMDS1,y=NMDS2,color=Site,shape=Depth))+
  geom_point(size=1.5, show.legend = TRUE)+
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
anosim(mgd_relabund.bray, sample_tab$Depth, permutations=999, distance="bray")


# PERMANOVA

adonis2(mgd_relabund.bray ~ metadata$Site*metadata$Depth, dist="bray",perm=999)
adonis2(mgd_relabund.bray ~ metadata$Site*metadata$Depth, dist="bray",perm=999)



pairwise.adonis2(mgd_relabund.bray ~ metadata$Site*metadata$Depth, data = metadata, sim.method = "bray", p.adjust.m = "BH", perm = 999)
#not working ?








#### 0-5cm and 5-15cm ONLY

#Import your metadata

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_metadata/ASCC_ITS_metadata_rem.txt",row.names=1) 

meta_0515 <- metadata %>% 
  filter(Depth == "0_5cm" | Depth == "5_15cm")

# Import feature table
otus <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG_rem.txt", row.names=1,header=TRUE)

otus_0515 <- otus %>%
  select(matches("_1C_(5|15)"))

# check that the names match - they must be in the same order! If all works the result of this line will read 'TRUE' 
all.equal(names(otus_0515),row.names(meta_0515)) #check to make sure row names and column names equal (aka the sample names), if this reads 'TRUE' you are good to move on!

# now we need to re-arrange the files, read in taxa, and make sure everything aligns before downstream steps, and make the phyloseq object
otumat<-as.matrix(otus_0515)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
taxa<-read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_taxonomy/taxonomy_ITS_ASCC.txt",row.names=1) #this is a list of the ASV ids and the corresponding taxa strings
taxmat<-as.matrix(taxa)
all.equal(row.names(taxmat),row.names(otumat)) # again - check to make sure they are in the same order!
TAX = tax_table(taxmat)
physeq<-phyloseq(OTU,TAX) #make your phyloseq object, which is basically just a R object that has all of your data for the analyses stored into it
all.equal(row.names(meta_0515),sample_names(physeq)) # Final check to make sure everything lines up
sampledata<-sample_data(meta_0515)
mgd<-merge_phyloseq(physeq,sampledata) #make the final phyloseq object

#Calculate relative abundance from ASV count data - dont need to run this step if you have already done so but it shouldnt really effect! 
mgd_relabund<-transform_sample_counts(mgd,function(x)x/sum(x)) #Calculate relative abundance from ASV count data - this is what you can save and use for other taxonomy analyses. If you're interested in doing these, let me know and i can help out!

# Calculate Bray-Curtis dissimilary (distance) 
mgd_relabund.bray<-distance(mgd_relabund,"bray") 

# Final line to generate NMDS ordination 
mgd_relabund.bray.nmds<-ordinate(mgd_relabund,"NMDS",mgd_relabund.bray) #WAHOO DO NMDS!!

# Use this to read out the NMDS stress (you generally want the stress to be below 0.2)
mgd_relabund.bray.nmds$stress 
# stress = 0.1244026 August2024

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
# hgd_view()

ggplot(sample_tab,aes(x=NMDS1,y=NMDS2,color=Site,shape=Depth))+
  geom_point(size=2, show.legend = TRUE)+
  #scale_shape_manual(values = c(17, 16)) +
  theme(text=element_text(size = 24))+  # Change axis labels
  xlim(-1.4,1.7)+
  ylim(-1.51,1.67)+
  #geom_text_repel(max.overlaps = 20, aes(label=row.names(map_file))) + #use this line if you would like to have the samples labled in your ordination 
  theme_classic() 
#scale_colour_manual(values = c("#ETRUE#scale_colour_maDepth#scale_colour_manual(values = c("#ETRUE#scale_colour_manual(values = c("#EE2A7B", "#27AAE1","#CCCCCB","#2E3192", "#F15A29","#009444", "#AC8A13", "#00A79D", "#7F8080", "#9E1F63", "#877DBB")) # feel free to add custom colors if you'd like! or remove this line

mrpp(mgd_relabund.bray, sample_tab$Site, permutations=999, distance="bray")
mrpp(mgd_relabund.bray, sample_tab$Depth, permutations=999, distance="bray")

anosim(mgd_relabund.bray, sample_tab$Site, permutations=999, distance="bray")
anosim(mgd_relabund.bray, sample_tab$Depth, permutations=999, distance="bray")

# PERMANOVA
# these both do the same thing
# adonis2(mgd_relabund.bray ~ meta_0515$Site*meta_0515$Depth, dist="bray",perm=999)
adonis2(mgd_relabund.bray ~ sample_tab$Site*sample_tab$Depth, dist="bray",perm=999)
adonis2(mgd_relabund.bray ~ sample_tab$Site, dist="bray",perm=999)
adonis2(mgd_relabund.bray ~ sample_tab$Depth, dist="bray",perm=999)

# not working
pairwise.adonis2(mgd_relabund.bray ~ sample_tab$Site*sample_tab$Depth, data = meta_0515, sim.method = "bray", p.adjust.m = "BH", perm = 999)






# beta disp:

# Assume 'group_factor' is a vector or column in your metadata that groups samples (e.g., Site or Treatment)
# Make sure it's a factor
# Ensure grouping matches row order of community_matrix
group_vector <- metadata$Ecosite

beta_disp <- betadisper(mgd_relabund.bray, group_vector)



# Significance test
permutest(beta_disp)

# Extract distances for plotting
disp_df <- data.frame(
  DistanceToCentroid = beta_disp$distances,
  Ecosite = group_vector
)

# Plot
library(ggplot2)
# reordering legends, refactoring
disp_df$Ecosite <- factor(disp_df$Ecosite, levels=c('SF_OM','SF_5','SF_15',
                                                     'TP_OM','TP_5','TP_15',
                                                     'SJ_OM','SJ_5','SJ_15'))

ggplot(disp_df, aes(x = Ecosite, y = DistanceToCentroid, fill = Ecosite)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5) +
  theme_minimal() +
  labs(title = "Beta Dispersion by Ecosite",
       y = "Distance to Centroid",
       x = "Ecosite") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()

ggplot(disp_df, aes(x = Ecosite, y = DistanceToCentroid, fill = Ecosite)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  facet_wrap(~ Ecosite, scales = "free_x", ncol = 3) +
  theme_minimal() +
  labs(title = "Beta Dispersion by Ecosite",
       x = NULL,
       y = "Distance to Centroid") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold"))+
  theme_bw()


anova(beta_disp)       # F-test for differences in dispersion
permutest(beta_disp)
TukeyHSD(beta_disp)   # Permutation-based test (more robust)

  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = distances ~ group, data = df)

$group
                     diff          lwr           upr     p adj
SF_5-SF_15   0.0021900156 -0.021781218  0.0261612494 0.9999987
SF_OM-SF_15  0.0072761950 -0.017567003  0.0321193933 0.9920002
SJ_15-SF_15 -0.0014757834 -0.023385425  0.0204338581 0.9999999
SJ_5-SF_15  -0.0041288149 -0.026219681  0.0179620508 0.9996724
SJ_OM-SF_15  0.0192801664 -0.002307730  0.0408680629 0.1224587
TP_15-SF_15 -0.0290111424 -0.052621669 -0.0054006157 0.0047011
TP_5-SF_15  -0.0243057337 -0.048091285 -0.0005201822 0.0409089
TP_OM-SF_15  0.0031561290 -0.021222670  0.0275349275 0.9999799
SF_OM-SF_5   0.0050861795 -0.018885054  0.0290574133 0.9991660
SJ_15-SF_5  -0.0036657989 -0.024581533  0.0172499354 0.9997974
SJ_5-SF_5   -0.0063188305 -0.027424325  0.0147866641 0.9907405
SJ_OM-SF_5   0.0170901509 -0.003488305  0.0376686063 0.1930397
TP_15-SF_5  -0.0312011580 -0.053892401 -0.0085099147 0.0007835
TP_5-SF_5   -0.0264957492 -0.049369053 -0.0036224458 0.0103317
TP_OM-SF_5   0.0009661134 -0.022523488  0.0244557147 1.0000000
SJ_15-SF_OM -0.0087519784 -0.030661620  0.0131576631 0.9448555
SJ_5-SF_OM  -0.0114050100 -0.033495876  0.0106858558 0.7969526
SJ_OM-SF_OM  0.0120039714 -0.009583925  0.0335918679 0.7227146
TP_15-SF_OM -0.0362873374 -0.059897864 -0.0126768108 0.0000863
TP_5-SF_OM  -0.0315819287 -0.055367480 -0.0077963772 0.0014286
TP_OM-SF_OM -0.0041200661 -0.028498865  0.0202587324 0.9998460
SJ_5-SJ_15  -0.0026530316 -0.021384135  0.0160780719 0.9999599
SJ_OM-SJ_15  0.0207559498  0.002620758  0.0388911411 0.0119921
TP_15-SJ_15 -0.0275353590 -0.048036697 -0.0070340211 0.0011784
TP_5-SJ_15  -0.0228299503 -0.043532616 -0.0021272851 0.0185575
TP_OM-SJ_15  0.0046319124 -0.016749708  0.0260135331 0.9990300
SJ_OM-SJ_5   0.0234089814  0.005055259  0.0417627037 0.0027126
TP_15-SJ_5  -0.0248823275 -0.045577226 -0.0041874291 0.0063541
TP_5-SJ_5   -0.0201769188 -0.041071279  0.0007174419 0.0677224
TP_OM-SJ_5   0.0072849439 -0.014282338  0.0288522262 0.9798155
TP_15-SJ_OM -0.0482913088 -0.068448435 -0.0281341824 0.0000000
TP_5-SJ_OM  -0.0435859001 -0.063947757 -0.0232240428 0.0000000
TP_OM-SJ_OM -0.0161240374 -0.037175845  0.0049277699 0.2914051
TP_5-TP_15   0.0047054087 -0.017789589  0.0272004065 0.9992483
TP_OM-TP_15  0.0321672714  0.009045890  0.0552886529 0.0006347
TP_OM-TP_5   0.0274618627  0.004161782  0.0507619437 0.0082802