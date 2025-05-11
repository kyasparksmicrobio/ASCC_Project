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




#######################################################
#######################################################
#######################################################
#######################################################
######## LOOKING AT INDIVIDUAL DEPTHS #################
#######################################################
#######################################################
#######################################################
#######################################################





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

shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
# deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")

genus_summary <- shallow %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_5_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Genus)

taxa_genus <- shallow %>%
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

shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")



genus_summary <- shallow %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_10_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

taxa_genus <- shallow %>%
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

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")

genus_summary <- deep %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_5_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Genus)

taxa_genus <- deep %>%
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

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")



genus_summary <- deep %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_10_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

taxa_genus <- deep %>%
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

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")

genus_summary <- OM %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_5_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Genus)

taxa_genus <- OM %>%
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

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
#OM <- subset(data_meta_taxa_clean, Depth == "OM")

# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")

genus_summary <- OM %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_10_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

taxa_genus <- OM %>%
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
  filter(Genus != "")

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
# deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")

genus_summary <- SF %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_5_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Genus)

taxa_genus <- SF %>%
  filter(Genus %in% top_5_genus)

# Plot the top 5 genus
ggplot(taxa_genus, aes(x = Depth, y = counts, fill = guild)) +
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

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")


genus_summary <- SF %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_10_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

taxa_genus <- SF %>%
  filter(Genus %in% top_10_genus)

# Plot the top 10 genus
ggplot(taxa_genus, aes(x = Depth, y = counts, fill = guild)) +
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

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")

genus_summary <- TP %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_5_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Genus)

taxa_genus <- TP %>%
  filter(Genus %in% top_5_genus)

# Plot the top 5 genus
ggplot(taxa_genus, aes(x = Depth, y = counts, fill = guild)) +
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

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
# OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
TP <- subset(data_meta_taxa_clean, Site == "TP")
# SJ <- subset(data_meta_taxa_clean, Site == "SJ")



genus_summary <- TP %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_10_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

taxa_genus <- TP %>%
  filter(Genus %in% top_10_genus)

# Plot the top 10 genus
ggplot(taxa_genus, aes(x = Depth, y = counts, fill = guild)) +
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

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
#OM <- subset(data_meta_taxa_clean, Depth == "OM")
# 
# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
SJ <- subset(data_meta_taxa_clean, Site == "SJ")

genus_summary <- SJ %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_5_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 5) %>%
  pull(Genus)

taxa_genus <- SJ %>%
  filter(Genus %in% top_5_genus)

# Plot the top 5 genus
ggplot(taxa_genus, aes(x = Depth, y = counts, fill = guild)) +
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

#shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
#deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
#OM <- subset(data_meta_taxa_clean, Depth == "OM")

# SF <- subset(data_meta_taxa_clean, Site == "SF")
# TP <- subset(data_meta_taxa_clean, Site == "TP")
SJ <- subset(data_meta_taxa_clean, Site == "SJ")

genus_summary <- SJ %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(total_counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_counts > 0)

top_10_genus <- genus_summary %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

taxa_genus <- SJ %>%
  filter(Genus %in% top_10_genus)

# Plot the top 10 genus
ggplot(taxa_genus, aes(x = Depth, y = counts, fill = guild)) +
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
