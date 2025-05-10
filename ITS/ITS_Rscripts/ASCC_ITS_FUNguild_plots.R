######### FUNGuild Script ASCC 2024 #########
######## all depths, guild #################


library(ggplot2)
library(tidyverse)
library(vegan) # For general ecology functions
library(phyloseq)
library(dplyr)

data <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG_rem.txt")

data[2:305] <- lapply(data[2:305], as.character)

data <- data %>%
  pivot_longer(cols = 2:305, names_to = "samples", values_to = "counts")

#write.csv(data,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG-TEST.csv")

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

ggplot(data_meta_taxa, aes(x=Site, y=counts, fill=guild)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
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







######## all depths, LOOKING AT GENUS, & guild #################
######### LOOKING AT TOP 5 TAXA ############

# Step 3: Subset the data frame to keep only rows with the top 5 genera
data <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG_rem.txt")

data[2:305] <- lapply(data[2:305], as.character)

data <- data %>%
  pivot_longer(cols = 2:305, names_to = "samples", values_to = "counts")

#write.csv(data,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG-TEST.csv")

taxa2 <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_FUNguild/taxonomy_ITS_ASCC.guilds_new.txt") %>% 
  select(1:9)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_metadata/ASCC_ITS_metadata_rem.txt") 

data_meta <- inner_join(data, metadata, by =c("samples"="samples"))

data_meta_taxa <- inner_join(data_meta, taxa2, by =c("OTUID"="OTU"))
data_meta_taxa$counts <- as.numeric(as.character(data_meta_taxa$counts))

data_meta_taxa_clean <- data_meta_taxa %>%
  filter(Genus != "")

genus_summary <- data_meta_taxa_clean %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_5_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Genus)

taxa_genus <- data_meta_taxa_clean %>%
  filter(Genus %in% top_5_genus)

# Plot the top 5 genus
ggplot(taxa_genus, aes(x = Site, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Phylum", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 phyla are shown on x-axis

ggplot(taxa_genus, aes(x = Genus, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Phylum", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 phyla are shown on x-axis





### GENUS / GUILD LOOKING FOR TOP 10 TAXA
data <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG_rem.txt")

data[2:305] <- lapply(data[2:305], as.character)

data <- data %>%
  pivot_longer(cols = 2:305, names_to = "samples", values_to = "counts")

#write.csv(data,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG-TEST.csv")

taxa2 <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_FUNguild/taxonomy_ITS_ASCC.guilds_new.txt") %>% 
  select(1:9)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_metadata/ASCC_ITS_metadata_rem.txt") 

data_meta <- inner_join(data, metadata, by =c("samples"="samples"))

data_meta_taxa <- inner_join(data_meta, taxa2, by =c("OTUID"="OTU"))
data_meta_taxa$counts <- as.numeric(as.character(data_meta_taxa$counts))

data_meta_taxa_clean <- data_meta_taxa %>%
  filter(Genus != "")

genus_summary <- data_meta_taxa_clean %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_10_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

taxa_genus <- data_meta_taxa_clean %>%
  filter(Genus %in% top_10_genus)

# Plot the top 10 genus
ggplot(taxa_genus, aes(x = Site, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Phylum", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 phyla are shown on x-axis

ggplot(taxa_genus, aes(x = Genus, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Phylum", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 phyla are shown on x-axis











######## all depths, LOOKING AT PHYLUM, & guild #################
######### LOOKING AT TOP 5 TAXA ############

# Step 3: Subset the data frame to keep only rows with the top 5 genera
data <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG_rem.txt")

data[2:305] <- lapply(data[2:305], as.character)

data <- data %>%
  pivot_longer(cols = 2:305, names_to = "samples", values_to = "counts")

#write.csv(data,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG-TEST.csv")

taxa2 <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_FUNguild/taxonomy_ITS_ASCC.guilds_new.txt") %>% 
  select(1:9)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_metadata/ASCC_ITS_metadata_rem.txt") 

data_meta <- inner_join(data, metadata, by =c("samples"="samples"))

data_meta_taxa <- inner_join(data_meta, taxa2, by =c("OTUID"="OTU"))
data_meta_taxa$counts <- as.numeric(as.character(data_meta_taxa$counts))

sum(is.na(data_meta_taxa_clean$counts)) # Should be 0 if conversion was perfect

# Remove rows with NA counts
data_meta_taxa_clean <- data_meta_taxa %>%
  filter(!is.na(counts))

data_meta_taxa_clean <- data_meta_taxa %>%
  filter(Phylum != "")

phyla_summary <- data_meta_taxa_clean %>%
  filter(Phylum != "" & !is.na(Phylum)) %>%  # Exclude empty strings and NA
  group_by(Phylum) %>%
  summarize(total_counts = sum(counts), .groups = 'drop') %>%
  filter(total_counts > 0)  # Remove phyla with zero total counts

# Display the phyla summary to check results
print(phyla_summary)

top_5_phyla <- phyla_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Phylum)

# Print the top 5 phyla to verify
print(top_5_phyla)

taxa_phyla <- data_meta_taxa_clean %>%
  filter(Phylum %in% top_5_phyla)

print(summary(taxa_phyla))
print(unique(taxa_phyla$Phylum))

check_counts <- taxa_phyla %>%
  filter(Phylum %in% c("Mucoromycota", "Rozellomycota")) %>%
  summarize(total_counts = sum(counts, na.rm = TRUE))

head(taxa_phyla)

# Check unique Phylum values and their respective counts
unique(taxa_phyla$Phylum)
table(taxa_phyla$Phylum)

print(check_counts)

ggplot(taxa_phyla, aes(x = Phylum, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)

# Plot the top 5 Phylum
ggplot(taxa_phyla, aes(x = Site, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Phylum", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 phyla are shown on x-axis

ggplot(taxa_phyla, aes(x = Phylum, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Phylum", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 phyla are shown on x-axis

phyla_summary <- data_meta_taxa %>%
  filter(Phylum != "") %>%  # Remove empty strings
  group_by(Phylum) %>%
  summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)  # Remove phyla with total counts of zero

# Now, select the top 5 phyla by total counts
top_5_phyla <- phyla_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Phylum)

# Filter the original data for only these top 5 phyla
taxa_phyla <- data_meta_taxa %>%
  filter(Phylum %in% top_5_phyla)

# Check the top 5 phyla for confirmation
print(top_5_phyla)

# check_individual_counts <- taxa_phyla %>%
#   filter(Phylum == "Mucoromycota" | Phylum == "Rozellomycota") %>%
#   group_by(Phylum) %>%
#   summarize(total_counts = sum(counts, na.rm = TRUE))
# 
# str(taxa_phyla)
# head(taxa_phyla)
# table(taxa_phyla$Phylum)
# 
# print(check_individual_counts)
# unique(taxa_phyla$Phylum)

ggplot(taxa_phyla, aes(x = Phylum, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE) +  # Ensure all phyla are shown on the x-axis
  labs(x = "Phylum", y = "Counts", title = "Top 5 Phyla by Total Counts")
