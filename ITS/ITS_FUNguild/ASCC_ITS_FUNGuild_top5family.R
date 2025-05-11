######### LOOKING AT TOP 5 TAXA ############
data <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_rem_RELABUND.txt")

data[2:305] <- lapply(data[2:305], as.character)

data <- data %>%
  pivot_longer(cols = 2:305, names_to = "samples", values_to = "relative_abundance")

#write.csv(data,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG-TEST.csv")

taxa2 <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_FUNguild/taxonomy_ITS_ASCC.guilds_new.txt") %>% 
  select(1:9)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_metadata/ASCC_ITS_metadata_rem.txt") 

data_meta <- inner_join(data, metadata, by =c("samples"="samples"))
data_meta_taxa <- inner_join(data_meta, taxa2, by =c("OTUID"="OTU"))
data_meta_taxa$relative_abundance <- as.numeric(as.character(data_meta_taxa$relative_abundance)) 


# Remove rows with NA counts
data_meta_taxa_clean <- data_meta_taxa %>%
  filter(!is.na(relative_abundance))

data_meta_taxa_clean <- data_meta_taxa %>%
  filter(Family != "")

family_summary <- data_meta_taxa_clean %>%
  filter(Family != "" & !is.na(Family)) %>%  # Exclude empty strings and NA
  group_by(Family) %>%
  summarize(relative_abundance = sum(relative_abundance), .groups = 'drop') %>%
  filter(relative_abundance > 0)  # Remove phyla with zero total counts

top_5_family <- family_summary %>%
  arrange(desc(relative_abundance)) %>%
  slice_head(n = 5) %>%
  pull(Family)

taxa_family <- data_meta_taxa %>%
  filter(Family %in% top_5_family)

ggplot(taxa_family, aes(x = Site, y = relative_abundance, fill = guild)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)


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

# Remove rows with NA counts
data_meta_taxa_clean <- data_meta_taxa %>%
  filter(!is.na(counts))

data_meta_taxa_clean <- data_meta_taxa %>%
  filter(Family != "")

family_summary <- data_meta_taxa_clean %>%
  filter(Family != "" & !is.na(Family)) %>%  # Exclude empty strings and NA
  group_by(Family) %>%
  summarize(total_counts = sum(counts), .groups = 'drop') %>%
  filter(total_counts > 0)  # Remove phyla with zero total counts

top_10_family <- family_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Family)

taxa_family <- data_meta_taxa %>%
  filter(Family %in% top_10_family)

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)

ggplot(taxa_family, aes(x = Site, y = counts, fill = guild)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)




#######################################################
#######################################################
#######################################################
#######################################################
######## LOOKING AT INDIVIDUAL DEPTHS #################
#######################################################
#######################################################
#######################################################
#######################################################





######## all depths, LOOKING AT family, & guild #################
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
  filter(Family != "")

shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
# deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")

family_summary <- shallow %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_5_family <- family_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Family)

taxa_family <- shallow %>%
  filter(Family %in% top_5_family)

# Plot the top 5 family
ggplot(taxa_family, aes(x = Site, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis





### family / GUILD LOOKING FOR TOP 10 TAXA
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
  filter(Family != "")

shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")



family_summary <- shallow %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_10_family <- family_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Family)

taxa_family <- shallow %>%
  filter(Family %in% top_10_family)

# Plot the top 10 family
ggplot(taxa_family, aes(x = Site, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis



######## all depths, LOOKING AT family, & guild #################
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
  filter(Family != "")

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")

family_summary <- deep %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_5_family <- family_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Family)

taxa_family <- deep %>%
  filter(Family %in% top_5_family)

# Plot the top 5 family
ggplot(taxa_family, aes(x = Site, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis


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
  filter(Family != "")

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")



family_summary <- deep %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_10_family <- family_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Family)

taxa_family <- deep %>%
  filter(Family %in% top_10_family)

# Plot the top 10 family
ggplot(taxa_family, aes(x = Site, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis








######## all depths, LOOKING AT FAMILY, & guild #################
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
  filter(Family != "")

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")

family_summary <- OM %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_5_family <- family_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Family)

taxa_family <- OM %>%
  filter(Family %in% top_5_family)

# Plot the top 5 family
ggplot(taxa_family, aes(x = Site, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis


### FAMILY / GUILD LOOKING FOR TOP 10 TAXA
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
  filter(Family != "")

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
OM <- subset(data_meta_taxa_clean, Depth == "OM")

# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")

family_summary <- OM %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_10_family <- family_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Family)

taxa_family <- OM %>%
  filter(Family %in% top_10_family)

# Plot the top 10 family
ggplot(taxa_family, aes(x = Site, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

#######################################################
#######################################################
#######################################################
#######################################################
######## LOOKING AT INDIVIDUAL SITES #################
#######################################################
#######################################################
#######################################################
#######################################################



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
  filter(Family != "")

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
# deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")

family_summary <- SF %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_5_family <- family_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Family)

taxa_family <- SF %>%
  filter(Family %in% top_5_family)

# Plot the top 5 family
ggplot(taxa_family, aes(x = Depth, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis





### FAMILY / GUILD LOOKING FOR TOP 10 TAXA
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
  filter(Family != "")

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")


family_summary <- SF %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_10_family <- family_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Family)

taxa_family <- SF %>%
  filter(Family %in% top_10_family)

# Plot the top 10 family
ggplot(taxa_family, aes(x = Depth, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis



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
  filter(Family != "")

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")

family_summary <- TP %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_5_family <- family_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Family)

taxa_family <- TP %>%
  filter(Family %in% top_5_family)

# Plot the top 5 family
ggplot(taxa_family, aes(x = Depth, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis


### FAMILY / GUILD LOOKING FOR TOP 10 TAXA
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
  filter(Family != "")

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")



family_summary <- TP %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_10_family <- family_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Family)

taxa_family <- TP %>%
  filter(Family %in% top_10_family)

# Plot the top 10 family
ggplot(taxa_family, aes(x = Depth, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis








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
  filter(Family != "")

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
#OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
SJ <- subset(data_meta_taxa_clean, Site == "SJ")

family_summary <- SJ %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_5_family <- family_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Family)

taxa_family <- SJ %>%
  filter(Family %in% top_5_family)

# Plot the top 5 family
ggplot(taxa_family, aes(x = Depth, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis


### FAMILY / GUILD LOOKING FOR TOP 10 TAXA
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
  filter(Family != "")

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
#OM <- subset(data_meta_taxa_clean, Depth == "OM")

# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
SJ <- subset(data_meta_taxa_clean, Site == "SJ")

family_summary <- SJ %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_10_family <- family_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Family)

taxa_family <- SJ %>%
  filter(Family %in% top_10_family)

# Plot the top 10 family
ggplot(taxa_family, aes(x = Depth, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

ggplot(taxa_family, aes(x = Family, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  #labs(x = "Family", y = "Counts", title = "Top 5 Phyla by Total Counts") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)  # Ensure all top 5 family are shown on x-axis

