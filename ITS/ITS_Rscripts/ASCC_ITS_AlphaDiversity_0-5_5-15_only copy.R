*** BOTH Depths ***
*** July 2024 ***
*** Redoing ASCC ITS alpha diversity processing! ***
*** be sure to edit metadata file to add location and depth ***
*** i have removed blanks!***
*** REMOVED MITOCHONDRIA AND CHLOROPLAST ***
*** code for relative abundance at the bottom ***

# install.packages("ggrepel")
# install.packages("phyloseq")
# install.packages("tidyverse")

# BiocManager::install("phyloseq")
#   if (!require("BiocManager", quietly = TRUE,))

# install.packages("BiocManager",force=TRUE)
# install.packages("grid")
# install.packages("dplyr")
# install.packages("reshape2")
# install.packages("vegan")
# install.packages("ggplot2")
# install.packages("Hmisc")
# install.packages("ggpubr")
# update.packages("ggplot2")
# update.packages("ggsignif")


# Load packages 
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
library(paletteer) # For custom color palettes (optional, if using specific palettes)


# Read in feature table
otus <- read.delim("/Users/kyasparks/Library/Mobile Documents/com~apple~CloudDocs/Desktop - MacBook Pro (2)/Kyas_PhD/CSU_2023_2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG_rem.txt") %>% 
  dplyr::rename(OTU_ID = "OTUID") %>%  
  as.data.frame()

row.names(otus) <- otus$OTU_ID

otus <- otus %>%
  dplyr::select(-OTU_ID)

t_otus <- t(otus)


length(row.names(t_otus))

alpha = as.data.frame(matrix(data = otus, nrow = 304, ncol = 4))

row.names(alpha) = row.names(t_otus)

alpha[,1] = diversity(t_otus, index = "shannon")
alpha[,2] = diversity(t_otus, index = "simpson")
alpha[,3] = diversity(t_otus, index = "shannon")/log(specnumber(t_otus)) 
alpha[,4] = specnumber(t_otus)

colnames(alpha) = c("Shannon", "Simpson", "Pielou", "Species_Richness")

write.csv(alpha,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_alpha_diversity/AlphaDiversity_ASCC_ITS_0-5_5-15_August2024_meta_added.csv") 

data <- read.csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_alpha_diversity/AlphaDiversity_ASCC_ITS_0-5_5-15_August2024_meta_added.csv", header = TRUE,row.names = 1)

data_0515 <- data %>%
  filter(grepl("1C_(5|15)$", rownames(data)))

my_comparisons<-list(c('SF','TP'),c('SF','SJ'),c('TP','SJ'))

# data$Depth <- factor(data$Depth, levels=c('OM','0_5cm','5_15cm'))
# data$Site <- factor(data$Site, levels=c('SF','TP','SJ'))

*** LOOKING AT SITE ***
 
ggplot(data_0515, aes(x=Site, y=Shannon)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data_0515, aes(x=Site, y=Shannon)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Depth, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()



ggplot(data_0515, aes(x=Site, y=Simpson)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data_0515, aes(x=Site, y=Simpson)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Depth, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data_0515, aes(x=Site, y=Pielou)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data_0515, aes(x=Site, y=Pielou)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Depth, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data_0515, aes(x=Site, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data_0515, aes(x=Site, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Depth, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

*** LOOKING AT DEPTH ***

ggplot(data, aes(x=Depth, y=Shannon)) + 
  geom_boxplot(aes(fill=Depth), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~Site, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons2, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Depth, y=Simpson)) + 
  geom_boxplot(aes(fill=Depth), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Site, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons2, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Depth, y=Pielou)) + 
  geom_boxplot(aes(fill=Depth), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Site, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons2, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Depth, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Depth), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Site, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons2, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()


***** SHALLOW ONLY -> SITE *****

shallow <- subset(data, Depth == "0_5cm")

ggplot(shallow, aes(x=Site, y=Shannon)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(shallow, aes(x=Site, y=Simpson)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(shallow, aes(x=Site, y=Pielou)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(shallow, aes(x=Site, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

***** DEEP ONLY -> SITE *****

deep <- subset(data, Depth == "5_15cm")

ggplot(deep, aes(x=Site, y=Shannon)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(deep, aes(x=Site, y=Simpson)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(deep, aes(x=Site, y=Pielou)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(deep, aes(x=Site, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

***** OM ONLY -> SITE *****

OM <- subset(data, Depth == "OM")

ggplot(OM, aes(x=Site, y=Shannon)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(OM, aes(x=Site, y=Simpson)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(OM, aes(x=Site, y=Pielou)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(OM, aes(x=Site, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()





*** Same code using relative abundance ***

# Read in feature table
otus <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_rem_RELABUND.txt") %>% 
  dplyr::rename(OTU_ID = "OTUID") %>%  
  as.data.frame()

row.names(otus) <- otus$OTU_ID

otus <- otus %>%
  dplyr::select(-OTU_ID)

t_otus <- t(otus)

length(row.names(t_otus))

alpha = as.data.frame(matrix(data = otus, nrow = 304, ncol = 4))

row.names(alpha) = row.names(t_otus)

alpha[,1] = diversity(t_otus, index = "shannon")
alpha[,2] = diversity(t_otus, index = "simpson")
alpha[,3] = diversity(t_otus, index = "shannon")/log(specnumber(t_otus)) 
alpha[,4] = specnumber(t_otus)

colnames(alpha) = c("Shannon", "Simpson", "Pielou", "Species_Richness")

write.csv(alpha,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_alpha_diversity/AlphaDiversity_ASCC_ITS_bothdepths_August2024_RELABUND.csv") 

data <- read.csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_alpha_diversity/AlphaDiversity_ASCC_ITS_bothdepths_August2024_RELABUND_meta_added.csv", header = TRUE,row.names = 1)

my_comparisons<-list(c('SF','TP'),c('SF','SJ'),c('TP','SJ'))

my_comparisons2<-list(c('0_5cm','5_15cm'),c('0_5cm','OM'),c('OM','5_15cm'))

*** LOOKING AT SITE ***
 
ggplot(data, aes(x=Site, y=Shannon)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Site, y=Shannon)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Depth, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Site, y=Simpson)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Site, y=Simpson)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Depth, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Site, y=Pielou)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Site, y=Pielou)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Depth, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Site, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Site, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Depth, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

*** LOOKING AT DEPTH ***

ggplot(data, aes(x=Depth, y=Shannon)) + 
  geom_boxplot(aes(fill=Depth), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Site, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons2, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Depth, y=Simpson)) + 
  geom_boxplot(aes(fill=Depth), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Site, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons2, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Depth, y=Pielou)) + 
  geom_boxplot(aes(fill=Depth), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Site, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons2, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(data, aes(x=Depth, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Depth), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Site, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons2, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

***** SHALLOW ONLY -> SITE *****

shallow <- subset(data, Depth == "0_5cm")

ggplot(shallow, aes(x=Site, y=Shannon)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(shallow, aes(x=Site, y=Simpson)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(shallow, aes(x=Site, y=Pielou)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(shallow, aes(x=Site, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

***** DEEP ONLY -> SITE *****

deep <- subset(data, Depth == "5_15cm")

ggplot(deep, aes(x=Site, y=Shannon)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(deep, aes(x=Site, y=Simpson)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(deep, aes(x=Site, y=Pielou)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(deep, aes(x=Site, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

***** OM ONLY -> SITE *****

OM <- subset(data, Depth == "OM")

ggplot(OM, aes(x=Site, y=Shannon)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(OM, aes(x=Site, y=Simpson)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(OM, aes(x=Site, y=Pielou)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()

ggplot(OM, aes(x=Site, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  #facet_wrap(~location2, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  theme_bw()
