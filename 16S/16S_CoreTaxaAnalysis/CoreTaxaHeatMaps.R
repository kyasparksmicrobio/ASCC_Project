library(tidyverse)
library(ggplot2)
library(grid)
library(dplyr)
library(vegan)
library(Hmisc)
library(readxl)
install.packages("readxl")
install.packages("readxl")
library(readxl)

core_taxa_60_data_16S <- read_excel("/Users/kyasparks/Desktop/Book1.xlsx",, sheet = d70, row.names=1) #this is a list of the ASV ids and the corresponding taxa strings
core_taxa_60_data_16S <- read_excel("/Users/kyasparks/Desktop/Book1.xlsx", sheet = "d70")


core <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_CoreTaxaAnalysis/16S_CoreTaxaAnalysis.xlsx", sheet = d70)


merged <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_denoise_files/ASCC_16S_feature_table_rem.txt") %>% 
  pivot_longer(names_to = "samples", values_to = "Number_of_ASVs",2:304)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_metadata_rem.txt") 

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, core_taxa_60_data_16S, by =c("OTUID"="OTUID"))

ggplot(data_meta_taxa, aes(x = Site, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_colour_paletteer_d("lisa::BridgetRiley")+
  scale_fill_paletteer_d("lisa::BridgetRiley")+
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Site, y = Order, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


core_taxa_data_70_16S <- read_csv("/Users/kyasparks/Desktop/Desktop - MacBook Pro (2)/Kyas_PhD/CSU_2023_2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_CoreTaxaAnalysis/CoreTaxa70_16S_ASCC.csv")

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, core_taxa_data_70_16S, by =c("OTUID"="OTUID"))

ggplot(data_meta_taxa, aes(x = Site, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Site, y = Order, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


core_taxa_data_70_5_16S <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_CoreTaxaAnalysis/CoreTaxa_70_16S_ASCC.xlsx",sheet = "d70_shallow")

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, core_taxa_data_70_5_16S, by =c("OTUID"="OTUID"))

ggplot(data_meta_taxa, aes(x = Site, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Site, y = Class, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




core_taxa_data_70_SF_16S <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_CoreTaxaAnalysis/CoreTaxa_70_16S_ASCC.xlsx",sheet = "d70_SF")

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, core_taxa_data_70_16S, by =c("OTUID"="OTUID"))

ggplot(data_meta_taxa, aes(x = Depth, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(data_meta_taxa, aes(x = Depth, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_colour_paletteer_d("lisa::BridgetRiley")+
  scale_fill_paletteer_d("lisa::BridgetRiley")+
  scale_x_continous()+
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Depth, y = Order, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


core_taxa_data_70_TP_16S <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_CoreTaxaAnalysis/CoreTaxa_70_16S_ASCC.xlsx",sheet = "d70_TP")

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, core_taxa_data_70_TP_16S, by =c("OTUID"="OTUID"))

ggplot(data_meta_taxa, aes(x = Depth, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Depth, y = Order, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


core_taxa_data_70_SJ_16S <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_CoreTaxaAnalysis/CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70_SJ")

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, core_taxa_data_70_SJ_16S, by =c("OTUID"="OTUID"))

ggplot(data_meta_taxa, aes(x = Depth, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Depth, y = Order, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

core_taxa_data_70_15_16S <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_CoreTaxaAnalysis/CoreTaxa_70_16S_ASCC.xlsx",sheet = "d70_deep")

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, core_taxa_data_70_15_16S, by =c("OTUID"="OTUID"))

ggplot(data_meta_taxa, aes(x = Site, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Site, y = Order, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




core_taxa_data_70_15_16S <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_CoreTaxaAnalysis/CoreTaxa_70_16S_ASCC.xlsx",sheet = "d70_deep")

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, core_taxa_data_70_15_16S, by =c("OTUID"="OTUID"))

ggplot(data_meta_taxa, aes(x = Site, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Site, y = Order, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))








core_taxa_data_70_OM_16S <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_CoreTaxaAnalysis/CoreTaxa_70_16S_ASCC.xlsx",sheet = "d70_OM")

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, core_taxa_data_70_OM_16S, by =c("OTUID"="OTUID"))

ggplot(data_meta_taxa, aes(x = Site, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Site, y = Order, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






core_taxa_data_70_OM_16S <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_CoreTaxaAnalysis/CoreTaxa_70_16S_ASCC.xlsx",sheet = "d70_OM")

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, core_taxa_data_70_OM_16S, by =c("OTUID"="OTUID"))

ggplot(data_meta_taxa, aes(x = Site, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Site, y = Order, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))











core_taxa_data_80_16S <- read_csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_CoreTaxaAnalysis/CoreTaxa80_16S_ASCC.csv")

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, core_taxa_data_80_16S, by =c("OTUID"="OTUID"))

ggplot(data_meta_taxa, aes(x = Site, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Site, y = Class, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##############################
########## ITS ###############
##############################

core_taxa_60_data_ITS <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_CoreTaxaAnalysis/CoreTaxa_60_ITS_ASCC.xlsx",sheet = "d60")


merged <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG_rem.txt") %>% 
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(names_to = "samples", values_to = "Number_of_ASVs",2:305)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_metadata/ASCC_ITS_metadata_rem.txt") 

data_meta <- inner_join(merged, metadata, by =c("samples"="samples"))
data_meta_taxa <- inner_join(data_meta, core_taxa_60_data_ITS, by =c("OTUID"="OTUID"))

data_meta_taxa$Number_of_ASVs <- as.numeric(data_meta_taxa$Number_of_ASVs)

ggplot(data_meta_taxa, aes(x = Depth, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 100),  # Set the limits of the scale
    breaks = seq(0, 100, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Site, y = Class, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


core_taxa_data_60_SJ_ITS <- read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_CoreTaxaAnalysis/CoreTaxa_60_ITS_ASCC.xlsx", sheet = "d60_SJ")

data_meta <- inner_join(merged, metadata, by =c("samples"="samples"))
data_meta_taxa <- inner_join(data_meta, core_taxa_data_60_SJ_ITS, by =c("OTUID"="OTUID"))
data_meta_taxa$Number_of_ASVs <- as.numeric(data_meta_taxa$Number_of_ASVs)

ggplot(data_meta_taxa, aes(x = Site, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 100),  # Set the limits of the scale
    breaks = seq(0, 100, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Depth, y = Class, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


core_taxa_data_80_ITS <- read_csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_CoreTaxaAnalysis/CoreTaxa80_ITS_ASCC.csv")

data_meta <- inner_join(merged, metadata, by =c("samples"="Sample_ID"))
data_meta_taxa <- inner_join(data_meta, core_taxa_data_80_ITS, by =c("OTUID"="OTUID"))

ggplot(data_meta_taxa, aes(x = Site, y = Phyla, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_meta_taxa, aes(x = Site, y = Class, fill = Number_of_ASVs)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),  # Customize colors
    limits = c(0, 1000),  # Set the limits of the scale
    breaks = seq(0, 1000, by = 200),  # Set the breaks for the color scale
    name = "Number of ASVs"  # Label for the color scale
  ) +
  theme_minimal() +
  labs(title = "Core Taxa Heatmap", x = "Site", y = "Phyla") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



