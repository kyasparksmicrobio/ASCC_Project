library(ggplot2)
library(tidyverse)
library(phyloseq)
library(dplyr)
library(readxl)

# data <- readx("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_feature_table_rem.txt")%>%
#   pivot_longer(names_to = "samples", values_to = "Number_of_ASVs",2:304)

data <- read_excel("/Users/kyasparks/Library/Mobile Documents/com~apple~CloudDocs/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/16S_Maaslin_TP.xlsx", sheet = "TP")%>%
  pivot_longer(names_to = "samples", values_to = "Relative Abundance",2:304)

write_csv(data,"/Users/kyasparks/Desktop/16S_Maaslin_TP.csv")

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_metadata_rem.txt") 
taxa <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_taxonomy/filtered/ASCC_16S_taxonomy_forprocessing_copy.txt")

data_meta <- inner_join(data, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, taxa, by =c("OTUID"="OTUID"))

#write.csv(data_meta_taxa,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_funguild_merged.csv")


shallow <- subset(data_meta_taxa, Depth == "0_5cm")
deep <- subset(data_meta_taxa, Depth == "5_15cm")
OM <- subset(data_meta_taxa, Depth == "OM")

SF <- subset(data_meta_taxa, Site == "SF")
TP <- subset(data_meta_taxa, Site == "TP")
SJ <- subset(data_meta_taxa, Site == "SJ")



### Filtering phyla column to only plot phyla with Number_of_ASVs ≥ 400
# filtered_data <- data_meta_taxa %>%
#   filter(Number_of_ASVs >= 1902) 

ggplot(otu_count_per_phylum, aes(x=Site, y=`Relative_Abundance`, fill=Phylum))+
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, size = ))

# ### Filtering Family column to only plot Family with Number_of_ASVs ≥ 1199
# 
# filtered_data2 <- data_meta_taxa %>%
#   filter(Number_of_ASVs >= 1199)
# 
# ggplot(filtered_data2, aes(x=Site, y=Number_of_ASVs, fill=Family)) +
#   geom_bar(stat="identity") +
#   #facet_wrap(~Site, nrow = 1, scales = "free_x")+
#   guides(fill=guide_legend(ncol = 1))+
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# ## 3 most common ASVs
# 
# ## multiple packages have functions with the same names, 
# ## which can lead to unexpected behavior if R is using the 
# ## wrong function
# 
# most_common_asvs <- data_meta_taxa %>%
#   dplyr::count(Site, OTUID) %>%
#   dplyr::arrange(Site, dplyr::desc(n)) %>%
#   dplyr::group_by(Site) %>%
#   dplyr::slice_max(order_by = n, n = 3)
# 
# most_common_asvs_taxa <- inner_join(most_common_asvs, taxa, by =c("OTUID"="OTUID"))
# most_common_asvs_taxa_meta <- inner_join(most_common_asvs_taxa, data_meta, by =c("OTUID"="OTUID"),relationship = "many-to-many")
# 
# # most_common_asvs_taxa_meta <- most_common_asvs_taxa_meta %>%
# #   mutate(
# #     OTUID = as.character(OTUID),
# #     Site = as.character(Site),
# #     Phylum = as.character(Phylum)
# #   )
# #colnames(most_common_asvs_taxa_meta) <- trimws(colnames(most_common_asvs_taxa_meta))
# #data_meta_taxa <- as_tibble(most_common_asvs_taxa_meta)
# 
# # cleaned_data <- most_common_asvs_taxa_meta %>%
# #   distinct(OTUID, Site.y, Phylum, .keep_all = TRUE)

ggplot(data_meta_taxa, aes(x=Site, y=`Relative Abundance`, fill=Phylum)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(shallow, aes(x=Site, y=`Relative Abundance`, fill=Phylum)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(deep, aes(x=Site, y=`Relative Abundance`, fill=Phylum)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(OM, aes(x=Site, y=`Relative Abundance`, fill=Phylum)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(SF, aes(x=Depth, y=`Relative Abundance`, fill=Phylum)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(TP, aes(x=Depth, y=`Relative Abundance`, fill=Phylum)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(SJ, aes(x=Depth, y=`Relative Abundance`, fill=Phylum)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data_meta_taxa, aes(x=Site, y=`Relative Abundance`, fill=Class)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# most_common_phyla <- most_common_asvs_taxa_meta %>%
#   count(Site, Phylum) %>%
#   arrange(Site, desc(n)) %>%
#   group_by(Site) %>%
#   slice_max(order_by = n, n = 3)
