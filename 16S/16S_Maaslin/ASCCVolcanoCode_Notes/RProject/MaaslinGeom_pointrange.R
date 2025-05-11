maaslin <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/ASCC_16S_Maaslin2_finalresultscombined.xlsx",sheet = "all")



# # Step 1: Identify Phyla Present in Every Site
phylum_site_count <- maaslin %>%
  group_by(Phylum) %>%
  summarise(Site_Count = n_distinct(Site)) %>%
  filter(Site_Count == length(unique(maaslin$Site))) %>%
  pull(Phylum)

# # Step 2: Subset the 'maaslin' Dataset
maaslin_subset <- maaslin %>%
  filter(Phylum %in% phylum_site_count)

# Count phylum occurrences grouped by site
phylum_counts <- maaslin_subset%>%
  group_by(Site, Phylum) %>%
  summarise(Phylum_Count = n()) %>%
  ungroup()

# Calculate total counts for each site
phylum_totals <- phylum_counts %>%
  group_by(Site) %>%
  summarise(Total_Count = sum(Phylum_Count))

# Add proportions by dividing phylum count by total count for each site
phylum_counts <- phylum_counts %>%
  left_join(phylum_totals, by = "Site") %>%
  mutate(Proportion = Phylum_Count / Total_Count)

# Join back to the original 'maaslin' dataframe
maaslin_subset <- maaslin_subset %>%
  left_join(phylum_counts, by = c("Site", "Phylum"))

maaslin_subset <- maaslin_subset %>%
  group_by(Site) %>%
  mutate(Site_Total = sum(Phylum_Count)) %>%
  ungroup() %>%
  mutate(Proportion = (Phylum_Count / Site_Total) * 100)

maaslin_subset <- maaslin_subset %>%
  group_by(Phylum, Site) %>%
  mutate(Proportion_min = Proportion - 0.1 * Proportion,  # Example: 5% lower than Proportion
         Proportion_max = Proportion + 0.1 * Proportion)

data <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_feature_table_rem.txt")%>%
  pivot_longer(names_to = "samples", values_to = "Number_of_ASVs",2:304)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_metadata_rem.txt") 
#taxa <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_taxonomy/filtered/ASCC_16S_taxonomy_forprocessing_copy.txt")

data_meta <- inner_join(data, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, maaslin_subset, by =c("OTUID"="OTUID", "Site"="Site"))


# phylum_stats <- maaslin_subset %>%
#   group_by(Phylum, Site) %>%
#   summarise(
#     Proportion_avg = mean(Proportion),
#     Proportion_min = min(Proportion),  
#     Proportion_max = max(Proportion))%>%
#   ungroup()
# 
# phylum_stats <- phylum_stats %>%
#   mutate(Proportion_min_adjusted = Proportion_min + -.1,
#          Proportion_max_adjusted = Proportion_max + .1)

# ggplot(phylum_counts, aes(x=Proportion, y=Phylum, fill=Phylum))+
#   geom_bar(stat="identity") +
#   #facet_wrap(~Site, nrow = 1, scales = "free_x")+
#   guides(fill=guide_legend(ncol = 1))+
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 75, hjust = 1, size = ))

ggplot(data_meta_taxa, aes(x = Proportion, y = Phylum, color = Site)) +
  geom_pointrange(aes(xmin = Proportion_min, xmax = Proportion_max), size = .3, position = position_dodge(width = 1)) +  # Position dodge for separation
  theme_minimal() +
  labs(title = "Phylum Proportions by Site", x = "Proportion", y = "Phylum") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("blue", "red", "green"))  # Customize colors for 3 sites

# ggplot(maaslin_subset, aes(x = Proportion, y = Phylum, color = Site)) +
#   geom_pointrange(aes(xmin = Proportion_min, xmax = Proportion_max), size = 0.3, position = position_dodge(width = 0.7)) +  # Position dodge for separation
#   theme_minimal() +
#   labs(title = "Phylum Proportions by Site",
#        x = "Average Proportion",
#        y = "Phylum") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_color_manual(values = c("blue", "red", "green")) 
#   #scale_x_continuous(labels = scales::percent) 





maaslin <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/ASCC_16S_Maaslin2_finalresultscombined.xlsx",sheet = "shallow")

maaslin <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/ASCC_16S_Maaslin2_finalresultscombined.xlsx",sheet = "deep")

maaslin <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/ASCC_16S_Maaslin2_finalresultscombined.xlsx",sheet = "OM")






maaslin.2 <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_maaslin/ASCC_ITS_Maaslin2_finalresultscombined.xlsx",sheet = "OM")

# # Step 1: Identify Phyla Present in Every Site
phylum_site_count <- maaslin.2 %>%
  group_by(Phyla) %>%
  summarise(Site_Count = n_distinct(Site)) %>%
  filter(Site_Count == length(unique(maaslin.2$Site))) %>%
  pull(Phyla)

# # Step 2: Subset the 'maaslin' Dataset
maaslin_subset <- maaslin.2 %>%
  filter(Phyla %in% phylum_site_count)

# Count phylum occurrences grouped by site
phylum_counts <- maaslin.2%>%
  group_by(Site, Phyla) %>%
  summarise(Phylum_Count = n()) %>%
  ungroup()

# Calculate total counts for each site
phylum_totals <- phylum_counts %>%
  group_by(Site) %>%
  summarise(Total_Count = sum(Phylum_Count))

# Add proportions by dividing phylum count by total count for each site
phylum_counts <- phylum_counts %>%
  left_join(phylum_totals, by = "Site") %>%
  mutate(Proportion = Phylum_Count / Total_Count)

# Join back to the original 'maaslin' dataframe
maaslin_subset <- maaslin.2 %>%
  left_join(phylum_counts, by = c("Site", "Phyla"))

maaslin_subset <- maaslin_subset %>%
  group_by(Site) %>%
  mutate(Site_Total = sum(Phylum_Count)) %>%
  ungroup() %>%
  mutate(Proportion = (Phylum_Count / Site_Total) * 100)

maaslin_subset <- maaslin_subset %>%
  group_by(Phyla, Site) %>%
  mutate(Proportion_min = Proportion - 0.1 * Proportion,  # Example: 5% lower than Proportion
         Proportion_max = Proportion + 0.1 * Proportion)

data <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG_rem.txt") %>%
  mutate(across(2:305, as.character)) %>%  # Convert columns 2 to 305 to character
  pivot_longer(cols = 2:305, names_to = "samples", values_to = "Number_of_ASVs")

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_metadata/ASCC_ITS_metadata_rem.txt") 
#taxa <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_taxonomy/filtered/ASCC_16S_taxonomy_forprocessing_copy.txt")

data_meta <- inner_join(data, metadata, by =c("samples"="samples"))
data_meta_taxa <- inner_join(data_meta, maaslin_subset, by =c("OTUID"="OTUID", "Site"="Site"))


# phylum_stats <- maaslin_subset %>%
#   group_by(Phylum, Site) %>%
#   summarise(
#     Proportion_avg = mean(Proportion),
#     Proportion_min = min(Proportion),  
#     Proportion_max = max(Proportion))%>%
#   ungroup()
# 
# phylum_stats <- phylum_stats %>%
#   mutate(Proportion_min_adjusted = Proportion_min + -.1,
#          Proportion_max_adjusted = Proportion_max + .1)

# ggplot(phylum_counts, aes(x=Proportion, y=Phylum, fill=Phylum))+
#   geom_bar(stat="identity") +
#   #facet_wrap(~Site, nrow = 1, scales = "free_x")+
#   guides(fill=guide_legend(ncol = 1))+
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 75, hjust = 1, size = ))

ggplot(data_meta_taxa, aes(x = Proportion, y = Phyla, color = Site)) +
  geom_pointrange(aes(xmin = Proportion_min, xmax = Proportion_max), size = .3, position = position_dodge(width = 1)) +  # Position dodge for separation
  theme_minimal() +
  labs(title = "Phylum Proportions by Site", x = "Proportion", y = "Phylum") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("blue", "red", "green"))  # Customize colors for 3 sites

# ggplot(maaslin_subset, aes(x = Proportion, y = Phylum, color = Site)) +
#   geom_pointrange(aes(xmin = Proportion_min, xmax = Proportion_max), size = 0.3, position = position_dodge(width = 0.7)) +  # Position dodge for separation
#   theme_minimal() +
#   labs(title = "Phylum Proportions by Site",
#        x = "Average Proportion",
#        y = "Phylum") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_color_manual(values = c("blue", "red", "green")) 
#   #scale_x_continuous(labels = scales::percent) 





maaslin <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/ASCC_16S_Maaslin2_finalresultscombined.xlsx",sheet = "shallow")

maaslin <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/ASCC_16S_Maaslin2_finalresultscombined.xlsx",sheet = "deep")

maaslin <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/ASCC_16S_Maaslin2_finalresultscombined.xlsx",sheet = "OM")


