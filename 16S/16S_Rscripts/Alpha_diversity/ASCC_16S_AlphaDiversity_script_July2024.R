###### BOTH Depths
#*** July 2024 ***
#*** Redoing ASCC ITS alpha diversity processing! ***
#*** be sure to edit metadata file to add location and depth ***
#*** i have removed blanks!***


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
install.packages("ggpubr")

# Load packages 
library(plyr)
library(tidyverse)
library(vegan) # For general ecology functions
library(GUniFrac) # For UniFrac
library(picante) # For bNTI and Faith's PD
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)
library(Hmisc)
library(reshape2)
library(httpgd)
library(languageserver)
hgd()
hgd_view()
hgd_browse()

# Read in feature table
otus <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_feature_table_rem.txt") %>% 
  dplyr::rename(OTU_ID = "OTUID") %>%  
  as.data.frame()

row.names(otus) <- otus$OTU_ID

otus <- otus %>%
  dplyr::select(-OTU_ID)

t_otus <- t(otus)

length(row.names(t_otus))

alpha = as.data.frame(matrix(data = otus, nrow = 303, ncol = 4))

row.names(alpha) = row.names(t_otus)

alpha[,1] = diversity(t_otus, index = "shannon")
alpha[,2] = diversity(t_otus, index = "simpson")
alpha[,3] = diversity(t_otus, index = "shannon")/log(specnumber(t_otus)) 
alpha[,4] = specnumber(t_otus)

colnames(alpha) = c("Shannon", "Simpson", "Pielou", "Species_Richness")

write.csv(alpha,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/alpha_diversity/AlphaDiversity_ASCC_16S_May2024.csv") 

data <- read.csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/alpha_diversity/AlphaDiversity_ASCC_16S_May2024.csv", header = TRUE,row.names = 1)

my_comparisons<-list(c('SF','TP'),c('SF','SJ'),c('TP','SJ'))

ggplot(data, aes(x=Site, y=Shannon)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Depth, y=Shannon)) + 
  geom_boxplot(aes(fill=Depth), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=location2, y=Species_Richness)) + 
  geom_boxplot(aes(fill=location2), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=location, y=Species_Richness)) + 
  geom_boxplot(aes(fill=location), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=location2, y=Species_Richness)) + 
  geom_boxplot(aes(fill=location2), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~depth, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=location2, y=Shannon)) + 
  geom_boxplot(aes(fill=location2), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~depth, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

shallow <- subset(data, depth == "0_5cm")

ggplot(shallow, aes(x=location2, y=Shannon)) + 
  geom_boxplot(aes(fill=location2), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()


ggplot(shallow, aes(x=location2, y=Species_Richness)) + 
  geom_boxplot(aes(fill=location2), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

deep <- subset(data, depth == "5_15cm")

ggplot(deep, aes(x=location2, y=Shannon)) + 
  geom_boxplot(aes(fill=location2), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(deep, aes(x=location2, y=Species_Richness)) + 
  geom_boxplot(aes(fill=location2), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

OM <- subset(data, depth == "OM")

ggplot(OM, aes(x=location2, y=Shannon)) + 
  geom_boxplot(aes(fill=location2), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(OM, aes(x=location2, y=Species_Richness)) + 
  geom_boxplot(aes(fill=location2), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

