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

shallow <- subset(data_meta_taxa_clean, Depth == "0_5cm")
deep <- subset(data_meta_taxa_clean, Depth == "5_15cm")
OM <- subset(data_meta_taxa_clean, Depth == "OM")

SF <- subset(data_meta_taxa_clean, Site == "SF")
TP <- subset(data_meta_taxa_clean, Site == "TP")
SJ <- subset(data_meta_taxa_clean, Site == "SJ")



ggplot(data_meta_taxa_clean, aes(x = Site, y = relative_abundance, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)

ggplot(shallow, aes(x = Site, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)

ggplot(deep, aes(x = Site, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)

ggplot(OM, aes(x = Site, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)

ggplot(SF, aes(x = Depth, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)

ggplot(TP, aes(x = Depth, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)

ggplot(SJ, aes(x = Depth, y = counts, fill = guild)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(drop = FALSE)
