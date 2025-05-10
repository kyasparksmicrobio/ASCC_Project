
library(ggplot2)
library(tidyverse)
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
# library(httpgd)
# library(languageserver)

# hgd()
# hgd_view()
# hgd_browse()


data <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG_rem.txt")%>%
  pivot_longer(names_to = "samples", values_to = "Number_of ASVs",2:305)

write.csv(data,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG-TEST.csv")

taxa2 <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_FUNguild/ASCC_FUNGUILD_SIMPLE.txt") %>% 
  select(1:10)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_metadata/ASCC_ITS_metadata_rem.txt") 


data_meta <- inner_join(data, metadata, by =c("samples"="samples"))
data_meta_taxa <- inner_join(data_meta, taxa2, by =c("OTUID"="OTU"))


shallow <- subset(data_meta_taxa, Depth == "0_5cm")
deep <- subset(data_meta_taxa, Depth == "5_15cm")
OM <- subset(data_meta_taxa, Depth == "OM")

SF <- subset(data_meta_taxa, Site == "SF")
TP <- subset(data_meta_taxa, Site == "TP")
SJ <- subset(data_meta_taxa, Site == "SJ")

ggplot(data_meta_taxa, aes(x=Site, y=Number_of_ASVs, fill=guild)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(shallow, aes(x=Site, y=counts, fill=guild)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(deep, aes(x=Site, y=counts, fill=guild)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(OM, aes(x=Site, y=counts, fill=guild)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(SF, aes(x=Depth, y=counts, fill=guild)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(TP, aes(x=Depth, y=counts, fill=guild)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(SJ, aes(x=Depth, y=counts, fill=guild)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



