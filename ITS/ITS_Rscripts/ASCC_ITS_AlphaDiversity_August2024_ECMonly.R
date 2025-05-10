*** BOTH Depths ***
*** July 2024 ***
*** Redoing ASCC ITS alpha diversity processing! ***
*** be sure to edit metadata file to add location and depth ***
*** i have removed blanks!***
*** REMOVED MITOCHONDRIA AND CHLOROPLAST ***
*** code for relative abundance at the bottom ***


# Load packages 
library(plyr)
library(tidyverse)
library(vegan) # For general ecology functionsor UniFrac
library(picante) # For bNTI and Faith's PD
library(phyloseq)
library(ggplot2)
library(dplyr)
library(GUniFrac) # For UniFrac
library(picante)
library(ggpubr) # For bNTI and Faith's PD


# Read in feature table
otus <- read_csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_FUNguild/ECMonly_featuretable.csv") %>% 
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

write.csv(alpha,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_alpha_diversity/AlphaDiversity_ASCC_ITS_bothdepths_August2024_ECM.csv") 

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_metadata/ASCC_ITS_metadata_rem.txt") 

data <- read.csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_alpha_diversity/AlphaDiversity_ASCC_ITS_bothdepths_August2024_ECM.csv",row.names = 1)

data_meta <- inner_join(data, metadata, by =c("samples"="samples"))

my_comparisons<-list(c('SF','TP'),c('SF','SJ'),c('TP','SJ'))

my_comparisons2<-list(c('0_5cm','5_15cm'),c('0_5cm','OM'),c('OM','5_15cm'))

*** LOOKING AT SITE ***

fam_colors <- c(
  "#3B1911",  # Dark Brown
  "#6D2F20",  # Red-Brown
  "#B75347",  # Light Red
  "#DF7666",  # Soft Red
  "#E09351",  # Orange-Brown
  "#EDC775",  # Yellow
  "#94B594",  # Soft Green
  "#6D928F",  # Blue-Green
  "#224B5E",  # Deep Blue
  "#11252E"   # Dark Blue
)
 
plot1 <- ggplot(data, aes(x=Site, y=Shannon)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Depth, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  scale_fill_manual(values = fam_colors[4:7]) +
  theme_bw()

ggsave(filename = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_alpha_diversity/shannons_ECMonly.svg",
       plot = plot1, width = 6, height = 6, dpi = 300, device = "svg")

plot2 <- ggplot(data, aes(x=Site, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Site), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Depth, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  scale_fill_manual(values = fam_colors[4:7]) +
  theme_bw()

ggsave(filename = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_alpha_diversity/speciesR_ECMonly.svg",
       plot = plot2, width = 6, height = 6, dpi = 300, device = "svg")

*** LOOKING AT DEPTH ***

plot3 <- ggplot(data, aes(x=Depth, y=Shannon)) + 
  geom_boxplot(aes(fill=Depth), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Site, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons2, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + scale_fill_manual(values = fam_colors[4:7]) +
  theme_bw()

ggsave(filename = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_alpha_diversity/shannons_ECMonly_bydepth.svg",
       plot = plot3, width = 6, height = 6, dpi = 300, device = "svg")

plot4 <- ggplot(data, aes(x=Depth, y=Species_Richness)) + 
  geom_boxplot(aes(fill=Depth), show.legend =TRUE) + 
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~Site, nrow = 1,scales = "free_x")+
  stat_compare_means(comparisons=my_comparisons2, method = "wilcox.test", label="p.signif", show.legend = FALSE, size = 6, hide.ns = FALSE, tip.length = 0.035, bracket.size = 0.5) + 
  scale_fill_manual(values = fam_colors[4:7]) +
  theme_bw()

ggsave(filename = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_alpha_diversity/speciesR_ECMonly_bydepth.svg",
       plot = plot4, width = 6, height = 6, dpi = 300, device = "svg")


print(plot1)
print(plot2)
print(plot3)
print(plot4)
