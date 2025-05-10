library(ggplot2)
library(tidyverse)
library(phyloseq)
library(dplyr)

merged <- read_csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_feature_table_rem_RELABUND.csv") %>% 
  pivot_longer(names_to = "samples", values_to = "Relative_Abundance",2:304)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_metadata_rem.txt") 

taxa <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_taxonomy/filtered/ASCC_16S_taxonomy_forprocessing_copy.txt")

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, taxa, by =c("OTUID"="OTUID"))

write.csv(data_meta_taxa,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_funguild_merged.csv")

avg_rel_abund <- data_meta_taxa %>%
  filter(!is.na(Phylum) & !is.na(Relative_Abundance)) %>%  # Ensure no missing values
  group_by(Phylum) %>%  # Group by Phylum
  dplyr::summarize(
    avg_relative_abundance = mean(Relative_Abundance, na.rm = TRUE),  # Calculate average
    .groups = "drop"
  ) %>%
  arrange(desc(avg_relative_abundance))

top_10_phyla <- avg_rel_abund %>%
  slice_max(order_by = avg_relative_abundance, n = 10)

top_10_phyla_names <- top_10_phyla$Phylum

#Print the top 10 Phyla for verification
print(top_10_phyla)



shallow <- subset(data_meta_taxa, Depth == "0_5cm")
deep <- subset(data_meta_taxa, Depth == "5_15cm")
OM <- subset(data_meta_taxa, Depth == "OM")
# 
SF <- subset(data_meta_taxa, Site == "SF")
TP <- subset(data_meta_taxa, Site == "TP")
SJ <- subset(data_meta_taxa, Site == "SJ")


avg_rel_abund_per_otuid <- SF %>%
  filter(!is.na(OTUID) & !is.na(Relative_Abundance)) %>%  # Ensure no missing values
  group_by(OTUID, Phylum) %>%  # Group by OTUID and Phylum
  dplyr::summarize(
    total_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances across samples
    num_samples = n_distinct(samples),  # Count the number of unique samples
    avg_relative_abundance = total_relative_abundance / num_samples,  # Calculate average relative abundance
    .groups = "drop"
  )


avg_rel_abund_per_phylum <- avg_rel_abund_per_otuid %>%
  group_by(Phylum) %>%
  dplyr::summarize(
    avg_relative_abundance = sum(avg_relative_abundance, na.rm = TRUE),  # Sum average relabundances across OTUs
    .groups = "drop"
  ) %>%
  arrange(desc(avg_relative_abundance))  


top_10_phyla <- avg_rel_abund_per_phylum %>%
  slice_max(order_by = avg_relative_abundance, n = 10)

# Extract the names of the top 10 Phyla
top_10_phyla_names <- top_10_phyla$Phylum

print(top_10_phyla)


# # Calculate average relative abundance for each Phylum and select the top 10
# top_10_phyla <- data_meta_taxa %>%
#   group_by(Phylum) %>%
#   summarize(Avg_Relative_Abundance = mean(Relative_Abundance, na.rm = TRUE)) %>%
#   arrange(desc(Avg_Relative_Abundance)) %>%  # Sort by highest average relative abundance
#   slice_max(order_by = Avg_Relative_Abundance, n = 10)  

# # Extract the names of the top 10 Phyla
# top_10_phyla_names <- top_10_phyla$Phylum


# Keep only OTUs belonging to the top 10 most abundant Phyla
data_top_10_phyla <- SF %>%
  filter(Phylum %in% top_10_phyla_names)

data_top_10_phyla <-TP %>%
  filter(Phylum %in% top_10_phyla_names)

data_top_10_phyla <-SJ %>%
  filter(Phylum %in% top_10_phyla_names)

### Filtering phyla column to only plot phyla with Number_of_ASVs ≥ 400
filtered_data <- data_top_10_phyla %>%
  filter(Relative_Abundance >= .1)

# install.packages("paletteer")
library(paletteer)

# Use in a ggplot2 chart:

library(ggplot2)

# Define custom colors for the Phyla
phyla_colors <- c(
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

# Create bar plot with custom Phylum colors
ggplot(filtered_data, aes(x = samples, y = Average_relative_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Depth, nrow = 1, scales = "free_x") +
  
  # Apply custom colors to Phyla
  scale_fill_manual(values = phy_colors) +
  
  guides(fill = guide_legend(ncol = 1)) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 10)
  )

### Filtering Family column to only plot Family with Number_of_ASVs ≥ 1199

filtered_data2 <- data_meta_taxa %>%
  filter(Number_of_ASVs >= 1199)

ggplot(filtered_data2, aes(x=Site, y=Number_of_ASVs, fill=Family)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## 3 most common ASVs

## multiple packages have functions with the same names, 
## which can lead to unexpected behavior if R is using the 
## wrong function

most_common_asvs <- data_meta_taxa %>%
  dplyr::count(Site, OTUID) %>%
  dplyr::arrange(Site, dplyr::desc(n)) %>%
  dplyr::group_by(Site) %>%
  dplyr::slice_max(order_by = n, n = 3)

most_common_asvs_taxa <- inner_join(most_common_asvs, taxa, by =c("OTUID"="OTUID"))
most_common_asvs_taxa_meta <- inner_join(most_common_asvs_taxa, data_meta, by =c("OTUID"="OTUID"),relationship = "many-to-many")

# most_common_asvs_taxa_meta <- most_common_asvs_taxa_meta %>%
#   mutate(
#     OTUID = as.character(OTUID),
#     Site = as.character(Site),
#     Phylum = as.character(Phylum)
#   )
#colnames(most_common_asvs_taxa_meta) <- trimws(colnames(most_common_asvs_taxa_meta))
#data_meta_taxa <- as_tibble(most_common_asvs_taxa_meta)

# cleaned_data <- most_common_asvs_taxa_meta %>%
#   distinct(OTUID, Site.y, Phylum, .keep_all = TRUE)


ggplot(, aes(x=Site.y, y=Number_of_ASVs, fill=Phylum)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

most_common_phyla <- most_common_asvs_taxa_meta %>%
  count(Site, Phylum) %>%
  arrange(Site, desc(n)) %>%
  group_by(Site) %>%
  slice_max(order_by = n, n = 3)


avg_rel_abund_per_otuid <- data_meta_taxa %>%
  filter(!is.na(OTUID) & !is.na(Relative_Abundance)) %>%  # Ensure no missing values
  group_by(OTUID, Phylum) %>%  # Group by OTUID and Phylum
  dplyr::summarize(
    total_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances across samples
    num_samples = n_distinct(samples),  # Count the number of unique samples
    avg_relative_abundance = total_relative_abundance / num_samples,  # Calculate average relative abundance
    .groups = "drop"
  )


avg_rel_abund_per_phylum <- avg_rel_abund_per_otuid %>%
  group_by(Phylum) %>%
  dplyr::summarize(
    avg_relative_abundance = sum(avg_relative_abundance, na.rm = TRUE),  # Sum average relabundances across OTUs
    .groups = "drop"
  ) %>%
  arrange(desc(avg_relative_abundance))  


top_10_phyla <- avg_rel_abund_per_phylum %>%
  slice_max(order_by = avg_relative_abundance, n = 10)

# Extract the names of the top 10 Phyla
top_10_phyla_names <- top_10_phyla$Phylum

print(top_10_phyla)


# Filter the dataset to include only the top 10 phyla
filtered_data <- SF %>%
  filter(Phylum %in% top_10_phyla_names) %>%  # Keep only rows with top 10 phyla
  group_by(Phylum, samples) %>%  # Group by Phylum and samples
  summarize(
    avg_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances for each sample
    .groups = "drop"
  )

# Define custom colors for the Phyla (adjust as needed)
phyla_colors <- c(
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

# Create the bar plot
ggplot(filtered_data, aes(x = samples, y = avg_relative_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phyla_colors) +  # Apply custom colors
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 10),
    legend.position = "right"
  ) +
  labs(
    title = "Top 10 Phyla by Average Relative Abundance",
    x = "Samples",
    y = "Average Relative Abundance",
    fill = "Phylum"
  )




  # Filter the dataset to include only the top 10 phyla
filtered_data <- SF %>%
  filter(Phylum %in% top_10_phyla_names) %>%  # Keep only rows with top 10 phyla
  group_by(Phylum, samples, Depth) %>%  # Group by Phylum, samples, and Depth
  summarize(
    Average_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances for each sample
    .groups = "drop"
  )

# Define custom colors for the Phyla (adjust as needed)
phyla_colors <- c(
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

# Create the bar plot with facets by Depth
ggplot(filtered_data, aes(x = samples, y = Average_relative_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Depth, nrow = 1, scales = "free_x") +  # Facet by Depth
  
  # Apply custom colors to Phyla
  scale_fill_manual(values = phyla_colors) +
  
  guides(fill = guide_legend(ncol = 1)) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 10)
  ) +
  labs(
    title = "Top 10 Phyla by Average Relative Abundance Across Depths",
    x = "Samples",
    y = "Average Relative Abundance",
    fill = "Phylum"
  )









library(ggplot2)
library(tidyverse)
library(phyloseq)
library(dplyr)

# Load and merge data
merged <- read_csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_feature_table_rem_RELABUND.csv") %>% 
  pivot_longer(names_to = "samples", values_to = "Relative_Abundance", 2:304)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_metadata_rem.txt") 

taxa <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_taxonomy/filtered/ASCC_16S_taxonomy_forprocessing_copy.txt")

data_meta <- inner_join(merged, metadata, by = c("samples" = "Sample_ID"))
data_meta_taxa <- inner_join(data_meta, taxa, by = c("OTUID" = "OTUID"))

# Calculate average relative abundance per OTUID
avg_rel_abund_per_otuid <- data_meta_taxa %>%
  filter(!is.na(OTUID) & !is.na(Relative_Abundance)) %>%  # Ensure no missing values
  group_by(OTUID, Phylum) %>%  # Group by OTUID and Phylum
  summarize(
    total_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances across samples
    num_samples = n_distinct(samples),  # Count the number of unique samples
    avg_relative_abundance = total_relative_abundance / num_samples,  # Calculate average relative abundance
    .groups = "drop"
  )

# Calculate average relative abundance per Phylum
avg_rel_abund_per_phylum <- avg_rel_abund_per_otuid %>%
  group_by(Phylum) %>%
  summarize(
    avg_relative_abundance = sum(avg_relative_abundance, na.rm = TRUE),  # Sum average relabundances across OTUs
    .groups = "drop"
  ) %>%
  arrange(desc(avg_relative_abundance))

# Select the top 10 Phyla
top_10_phyla <- avg_rel_abund_per_phylum %>%
  slice_max(order_by = avg_relative_abundance, n = 10)

# Extract the names of the top 10 Phyla
top_10_phyla_names <- top_10_phyla$Phylum

# Filter the dataset to include only the top 10 Phyla
filtered_data <- data_meta_taxa %>%
  filter(Phylum %in% top_10_phyla_names) %>%  # Keep only rows with top 10 Phyla
  group_by(Phylum, samples, Depth) %>%  # Group by Phylum, samples, and Depth
  summarize(
    avg_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances for each sample
    .groups = "drop"
  )

# Define custom colors for the Phyla
phyla_colors <- c(
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

# Create the bar plot with facets by Depth
ggplot(filtered_data, aes(x = samples, y = avg_relative_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Depth, nrow = 1, scales = "free_x") +  # Facet by Depth
  
  # Apply custom colors to Phyla
  scale_fill_manual(values = phyla_colors) +
  
  guides(fill = guide_legend(ncol = 1)) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 10)
  ) +
  labs(
    title = "Top 10 Phyla by Average Relative Abundance Across Depths",
    x = "Samples",
    y = "Average Relative Abundance",
    fill = "Phylum"
  )








library(ggplot2)
library(tidyverse)
library(phyloseq)
library(dplyr)

merged <- read_csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_feature_table_rem_RELABUND.csv") %>% 
  pivot_longer(names_to = "samples", values_to = "Relative_Abundance",2:304)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_metadata_rem.txt") 

taxa <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_taxonomy/filtered/ASCC_16S_taxonomy_forprocessing_copy.txt")

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, taxa, by =c("OTUID"="OTUID"))

#write.csv(data_meta_taxa,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_funguild_merged.csv")

# avg_rel_abund <- data_meta_taxa %>%
#   filter(!is.na(Phylum) & !is.na(Relative_Abundance)) %>%  # Ensure no missing values
#   group_by(Phylum) %>%  # Group by Phylum
#   dplyr::summarize(
#     avg_relative_abundance = mean(Relative_Abundance, na.rm = TRUE),  # Calculate average
#     .groups = "drop"
#   ) %>%
#   arrange(desc(avg_relative_abundance))

# top_10_phyla <- avg_rel_abund %>%
#   slice_max(order_by = avg_relative_abundance, n = 10)

# top_10_phyla_names <- top_10_phyla$Phylum

# Print the top 10 Phyla for verification
# print(top_10_phyla)



shallow <- subset(data_meta_taxa, Depth == "0_5cm")
deep <- subset(data_meta_taxa, Depth == "5_15cm")
OM <- subset(data_meta_taxa, Depth == "OM")
# 
SF <- subset(data_meta_taxa, Site == "SF")
TP <- subset(data_meta_taxa, Site == "TP")
SJ <- subset(data_meta_taxa, Site == "SJ")


avg_rel_abund_per_otuid <- SF %>%
  filter(!is.na(OTUID) & !is.na(Relative_Abundance)) %>%  # Ensure no missing values
  group_by(OTUID, Phylum) %>%  # Group by OTUID and Phylum
  dplyr::summarize(
    total_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances across samples
    num_samples = n_distinct(samples),  # Count the number of unique samples
    avg_relative_abundance = total_relative_abundance / num_samples,  # Calculate average relative abundance
    .groups = "drop"
  )


avg_rel_abund_per_phylum <- avg_rel_abund_per_otuid %>%
  group_by(Phylum) %>%
  dplyr::summarize(
    avg_relative_abundance = sum(avg_relative_abundance, na.rm = TRUE),  # Sum average relabundances across OTUs
    .groups = "drop"
  ) %>%
  arrange(desc(avg_relative_abundance))  


top_10_phyla <- avg_rel_abund_per_phylum %>%
  slice_max(order_by = avg_relative_abundance, n = 10)

# Extract the names of the top 10 Phyla
top_10_phyla_names <- top_10_phyla$Phylum

print(top_10_phyla)


# # Calculate average relative abundance for each Phylum and select the top 10
# top_10_phyla <- data_meta_taxa %>%
#   group_by(Phylum) %>%
#   summarize(Avg_Relative_Abundance = mean(Relative_Abundance, na.rm = TRUE)) %>%
#   arrange(desc(Avg_Relative_Abundance)) %>%  # Sort by highest average relative abundance
#   slice_max(order_by = Avg_Relative_Abundance, n = 10)  

# # Extract the names of the top 10 Phyla
# top_10_phyla_names <- top_10_phyla$Phylum


# Keep only OTUs belonging to the top 10 most abundant Phyla
data_top_10_phyla <- SF %>%
  filter(Phylum %in% top_10_phyla_names)

data_top_10_phyla <-TP %>%
  filter(Phylum %in% top_10_phyla_names)

data_top_10_phyla <-SJ %>%
  filter(Phylum %in% top_10_phyla_names)

### Filtering phyla column to only plot phyla with Number_of_ASVs ≥ 400
filtered_data <- data_top_10_phyla %>%
  filter(Relative_Abundance >= .1)

# install.packages("paletteer")
library(paletteer)

# Use in a ggplot2 chart:

library(ggplot2)


filtered_data <- filtered_data %>%
  filter(Phylum %in% top_10_phyla_names) %>%  # Keep only rows with top 10 phyla
  group_by(Phylum, samples) %>%  # Group by Phylum and samples
  summarize(
    avg_relative_abundance = sum(avg_relative_abundance, na.rm = TRUE),  # Sum relative abundances for each sample
    .groups = "drop"
  )







library(ggplot2)
library(tidyverse)
library(phyloseq)
library(dplyr)

# Load and merge data
merged <- read_csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_feature_table_rem_RELABUND.csv") %>% 
  pivot_longer(names_to = "samples", values_to = "Relative_Abundance", 2:304)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_metadata_rem.txt") 

taxa <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_taxonomy/filtered/ASCC_16S_taxonomy_forprocessing_copy.txt")

data_meta <- inner_join(merged, metadata, by = c("samples" = "Sample_ID"))
data_meta_taxa <- inner_join(data_meta, taxa, by = c("OTUID" = "OTUID"))

# Step 1: Calculate average relative abundance per OTUID for the entire dataset
avg_rel_abund_per_otuid <- data_meta_taxa %>%
  filter(!is.na(OTUID) & !is.na(Relative_Abundance)) %>%  # Ensure no missing values
  group_by(OTUID, Phylum) %>%  # Group by OTUID and Phylum
  summarize(
    total_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances across samples
    num_samples = n_distinct(samples),  # Count the number of unique samples
    avg_relative_abundance = total_relative_abundance / num_samples,  # Calculate average relative abundance
    .groups = "drop"
  )

# Step 2: Calculate average relative abundance per Phylum for the entire dataset
avg_rel_abund_per_phylum <- avg_rel_abund_per_otuid %>%
  group_by(Phylum) %>%
  summarize(
    avg_relative_abundance = sum(avg_relative_abundance, na.rm = TRUE),  # Sum average relabundances across OTUs
    .groups = "drop"
  ) %>%
  arrange(desc(avg_relative_abundance))

# Step 3: Identify the top 10 phyla based on the entire dataset
top_10_phyla <- avg_rel_abund_per_phylum %>%
  slice_max(order_by = avg_relative_abundance, n = 10)

# Extract the names of the top 10 Phyla
top_10_phyla_names <- top_10_phyla$Phylum

# Step 4: Subset by site (e.g., SF, TP, SJ) and filter for the top 10 phyla
filtered_data <- data_meta_taxa %>%
  filter(Site %in% c("SF", "TP", "SJ") & Phylum %in% top_10_phyla_names) %>%  # Subset by site and top 10 phyla
  group_by(Site, Phylum, samples) %>%  # Group by Site, Phylum, and samples
  summarize(
    avg_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances for each sample
    .groups = "drop"
  )

# Define custom colors for the Phyla
phyla_colors <- c(
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

# Step 5: Create the bar plot with facets by Site
ggplot(filtered_data, aes(x = samples, y = avg_relative_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Site, nrow = 1, scales = "free_x") +  # Facet by Site
  
  # Apply custom colors to Phyla
  scale_fill_manual(values = phy_colors) +
  
  guides(fill = guide_legend(ncol = 1)) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 10)
  ) +
  labs(
    title = "Top 10 Phyla by Average Relative Abundance Across Sites",
    x = "Samples",
    y = "Average Relative Abundance",
    fill = "Phylum"
  )
















library(ggplot2)
library(tidyverse)
library(dplyr)

# Step 1: Subset data for each site and include Depth
sf_data <- data_meta_taxa %>%
  filter(Site == "SF" & Phylum %in% top_10_phyla_names) %>%  # Subset for SF and top 10 phyla
  mutate(Depth = factor(Depth, levels = c("OM", "0_5cm", "5_15cm"))) %>%  # Reorder Depth levels
  group_by(Phylum, samples, Depth) %>%  # Group by Phylum, samples, and Depth
  summarize(
    avg_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances for each sample
    .groups = "drop"
  )

tp_data <- data_meta_taxa %>%
  filter(Site == "TP" & Phyla %in% top_10_phyla_names) %>%  # Subset for TP and top 10 phyla
  mutate(Depth = factor(Depth, levels = c("OM", "0_5cm", "5_15cm"))) %>%  # Reorder Depth levels
  group_by(Phyla, samples, Depth) %>%  # Group by Phylum, samples, and Depth
  summarize(
    avg_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances for each sample
    .groups = "drop"
  )

sj_data <- data_meta_taxa %>%
  filter(Site == "SJ" & Phyla %in% top_10_phyla_names) %>%  # Subset for SJ and top 10 phyla
  mutate(Depth = factor(Depth, levels = c("OM", "0_5cm", "5_15cm"))) %>%  # Reorder Depth levels
  group_by(Phyla, samples, Depth) %>%  # Group by Phylum, samples, and Depth
  summarize(
    avg_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances for each sample
    .groups = "drop"
  )

# Step 2: Define custom colors for the Phyla
phyla_colors <- c(
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

# Step 3: Create separate plots for each site, faceted by Depth
sf_plot <- ggplot(sf_data, aes(x = samples, y = avg_relative_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Depth, nrow = 1, scales = "free_x") +  # Facet by Depth
  scale_fill_manual(values = phy_colors) +  # Apply custom colors
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 10),
    legend.position = "right"
  ) +
  labs(
    title = "Top 10 Phyla by Average Relative Abundance (SF)",
    x = "Samples",
    y = "Average Relative Abundance",
    fill = "Phylum"
  )

tp_plot <- ggplot(tp_data, aes(x = samples, y = avg_relative_abundance, fill = Phyla)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Depth, nrow = 1, scales = "free_x") +  # Facet by Depth
  scale_fill_manual(values = phyla_colors) +  # Apply custom colors
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 10),
    legend.position = "right"
  ) +
  labs(
    title = "Top 10 Phyla by Average Relative Abundance (TP)",
    x = "Samples",
    y = "Average Relative Abundance",
    fill = "Phyla"
  )

sj_plot <- ggplot(sj_data, aes(x = samples, y = avg_relative_abundance, fill = Phyla)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Depth, nrow = 1, scales = "free_x") +  # Facet by Depth
  scale_fill_manual(values = phyla_colors) +  # Apply custom colors
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 10),
    legend.position = "right"
  ) +
  labs(
    title = "Top 10 Phyla by Average Relative Abundance (SJ)",
    x = "Samples",
    y = "Average Relative Abundance",
    fill = "Phyla"
  )

# Step 4: Print the plots
print(sf_plot)
print(tp_plot)
print(sj_plot)














library(ggplot2)
library(tidyverse)
library(dplyr)



# Load and merge data
merged <- read_csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_rem_RELABUND.csv") %>% 
  pivot_longer(names_to = "samples", values_to = "Relative_Abundance", 2:305)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_metadata/ASCC_ITS_metadata_rem.txt") 

taxa <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_taxonomy/taxonomy_ITS_ASCC.txt")

data_meta <- inner_join(merged, metadata, by = c("samples" = "samples"))
data_meta_taxa <- inner_join(data_meta, taxa, by = c("OTUID" = "OTUID"))

# Step 1: Subset data for each site and include Depth
sf_data <- data_meta_taxa %>%
  filter(Site == "SF" & Phylum %in% top_10_phyla_names) %>%  # Subset for SF and top 10 phyla
  mutate(Depth = factor(Depth, levels = c("OM", "0_5cm", "5_15cm"))) %>%  # Reorder Depth levels
  group_by(Phylum, samples, Depth) %>%  # Group by Phylum, samples, and Depth
  summarize(
    avg_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances for each sample
    .groups = "drop"
  )

tp_data <- data_meta_taxa %>%
  filter(Site == "TP" & Phyla %in% top_10_phyla_names) %>%  # Subset for TP and top 10 phyla
  mutate(Depth = factor(Depth, levels = c("OM", "0_5cm", "5_15cm"))) %>%  # Reorder Depth levels
  group_by(Phyla, samples, Depth) %>%  # Group by Phylum, samples, and Depth
  summarize(
    avg_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances for each sample
    .groups = "drop"
  )

sj_data <- data_meta_taxa %>%
  filter(Site == "SJ" & Phyla %in% top_10_phyla_names) %>%  # Subset for SJ and top 10 phyla
  mutate(Depth = factor(Depth, levels = c("OM", "0_5cm", "5_15cm"))) %>%  # Reorder Depth levels
  group_by(Phyla, samples, Depth) %>%  # Group by Phylum, samples, and Depth
  summarize(
    avg_relative_abundance = sum(Relative_Abundance, na.rm = TRUE),  # Sum relative abundances for each sample
    .groups = "drop"
  )

# Step 2: Define custom colors for the Phyla
phyla_colors <- c(
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

# Step 3: Create separate plots for each site, faceted by Depth
sf_plot <- ggplot(sf_data, aes(x = samples, y = avg_relative_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Depth, nrow = 1, scales = "free_x") +  # Facet by Depth
  scale_fill_manual(values = phy_colors) +  # Apply custom colors
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 10),
    legend.position = "right"
  ) +
  labs(
    title = "Top 10 Phyla by Average Relative Abundance (SF)",
    x = "Samples",
    y = "Average Relative Abundance",
    fill = "Phylum"
  )

tp_plot <- ggplot(tp_data, aes(x = samples, y = avg_relative_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Depth, nrow = 1, scales = "free_x") +  # Facet by Depth
  scale_fill_manual(values = phyla_colors) +  # Apply custom colors
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 10),
    legend.position = "right"
  ) +
  labs(
    title = "Top 10 Phyla by Average Relative Abundance (TP)",
    x = "Samples",
    y = "Average Relative Abundance",
    fill = "Phylum"
  )

sj_plot <- ggplot(sj_data, aes(x = samples, y = avg_relative_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Depth, nrow = 1, scales = "free_x") +  # Facet by Depth
  scale_fill_manual(values = phyla_colors) +  # Apply custom colors
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 10),
    legend.position = "right"
  ) +
  labs(
    title = "Top 10 Phyla by Average Relative Abundance (SJ)",
    x = "Samples",
    y = "Average Relative Abundance",
    fill = "Phylum"
  )

# Step 4: Print the plots
print(sf_plot)
print(tp_plot)
print(sj_plot)