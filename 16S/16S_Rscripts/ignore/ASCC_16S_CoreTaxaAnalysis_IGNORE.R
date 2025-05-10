## Code template from Dr. Emily Bechtold / Julie Fowler
## Doing these analyses by both taxa (collapse to genus or species, rows combined so each taxa only represented once) and by unique ASV - one file for each of these analyses (genus, species, ASVs)
## This is based on presence/absence, rather than any abundance calculations
## Went forward using 50% ASV cutoff based on “Defining and quantifying the core microbiome: Challenges and prospects” by Neu et al. for the 50%, ASV as that works for our questions
## I have an R markdown file with more notes in one of my folders

### ITS ###

library(tidyverse)
library(ggplot2)
library(grid)
library(dplyr)
library(vegan)
library(Hmisc)

merged <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG_rem.txt")

merged <- merged %>%
  mutate(across(everything(), as.character))

# Function to insert "0_" after the second underscore
update_column_names <- function(col_name) {
  if (str_detect(col_name, "SF")) {
    parts <- str_split(col_name, "_", simplify = TRUE)
    
    if (length(parts) >= 3) {
      # Ensure no extra underscores by combining parts with the right placement of "0_"
      new_name <- str_c(parts[1], parts[2], "0", paste(parts[3:length(parts)], collapse = "_"), sep = "_")
      return(new_name)
    } else {
      return(col_name)
    }
  } else {
    return(col_name)
  }
}

# Update column names
colnames(merged) <- sapply(colnames(merged), update_column_names)

##metadata <- read_csv(“metadata_genus_all_merged.csv”)

##metadata$Depth <- as.character(metadata$Depth)
# d1 <- merged %>%
#   pivot_longer(2:305, names_to = "Samples", values_to = "Counts") %>%
#   separate(Samples, into = c("SampleID", "Plot", "Unit", "Date", "Site", "NumOfCores", "Depth"), sep = "_", fill = "right") %>%
#   mutate(Value = ifelse(Counts > 0, 1, 0))

# d1 <- merged %>%
#   pivot_longer(cols = 2:305, names_to = "Samples", values_to = "Counts") %>%
#   separate(Samples, into = c("SampleID", "Plot", "Unit", "Date", "Site", "NumOfCores", "Depth"), sep = "_") %>%
#   mutate(Counts = as.numeric(Counts)) %>%  # Convert Counts to numeric
#   mutate(Value = ifelse(is.na(Counts), 0, 1))

d1 <- merged %>%
  pivot_longer(2:305, names_to = "Samples", values_to = "Counts") %>%
  separate(Samples, into = c("SampleID", "Plot", "Unit", "Date", "Site", "NumOfCores", "Depth"), sep = "_") %>%
  mutate(Counts = as.numeric(Counts)) %>% 
  mutate(Value=Counts/Counts)

# d1 <- merged %>%
#   pivot_longer(2:305, names_to = "Samples", values_to = "Counts") %>%
#   separate(Samples, into = c("SampleID", "Plot", "Unit", "Date", "Site", "NumOfCores", "Depth"), sep="_") %>%
#   mutate(Value = ifelse(Counts > 0, 1, 0))
#   
#d1 <- d1[,c("SampleID", "Plot", "Unit", "Date", "Site", "NumOfCores", "Depth")]
d1 <- d1 %>%
  mutate_all(~ ifelse(is.na(.), 0, .))

d1$NewColumn <- 1

## Different Cutoffs to be Considered Core Microbiome
## May need to use different cutoffs for different depths, etc. due to sample size differences?
## Ash: sample size of 3 common in each site, so 66% or 100% is logical
## See what other related papers do (notes in PowerPoint) - using 100% and 50% as high/low demonstrations and as the two that seem common in the literature

d30 <- d1 %>%
  group_by(OTUID,Site, Depth) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= 0.3)

d30b <- d1 %>%
  group_by(OTUID,Site) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= .3)

d50 <- d1 %>%
  group_by(OTUID, Site, Depth) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= 0.5)

d50b <- d1 %>%
  group_by(OTUID, Site) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= 0.5)

d60 <- d1 %>%
  group_by(OTUID, Site, Depth) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= 0.6)

d60b <- d1 %>%
  group_by(OTUID, Site) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= 0.6)

## 66% in a given combo of site/depth
d66 <- d1 %>%
  group_by(OTUID, Site, Depth) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= 0.66)

d66b <- d1 %>%
  group_by(OTUID, Site) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= 0.66)

d70 <- d1 %>%
  group_by(OTUID, Site, Depth) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= 0.7)

d70b <- d1 %>%
  group_by(OTUID, Site) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= 0.7)

## 80% in a given combo of site/depth
d80 <- d1 %>%
  group_by(OTUID, Site, Depth) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= 0.8)

d80b <- d1 %>%
  group_by(OTUID, Site) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= 0.8)

d90 <- d1 %>%
  group_by(OTUID, Site, Depth) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= 0.9)

# d90b <- d1 %>%
#   group_by(OTUID, Site) %>%
#   summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
#   arrange(desc(presence)) %>%
#   filter(presence >= 0.9)


## 100% in a given combo of site/depth
d100 <- d1 %>%
  group_by(OTUID,Site, Depth) %>%
  summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
  arrange(desc(presence)) %>%
  filter(presence >= 1)

# d100b <- d1 %>%
#   group_by(OTUID, Site) %>%
#   summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
#   arrange(desc(presence)) %>%
#   filter(presence >= 1)

# d01 <- d1 %>%
#   group_by(OTUID, Site, Depth) %>%
#   summarise(presence = sum(Value == 1) / n(), .groups = "drop") %>%
#   arrange(desc(presence)) %>%
#   filter(presence >= .01)



## 100% Cutoff - Looking for taxa that are present in at least 100% of a given site/depth combo in 100% of the sites
SF_5 <- d100 %>%
  filter(Site == "SF", Depth == "5")

TP_5 <- d100 %>%
  filter(Site == "TP", Depth == "5")

SJ_5 <- d100 %>%
  filter(Site == "SJ", Depth == "5")

SFTPSJ_5 <- rbind(SF_5,TP_5,SJ_5)


write.csv(SFTPSJ_5,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_5_100_coretaxa.csv")


SF_15 <- d100 %>%
  filter(Site == "SF", Depth == "15")

TP_15 <- d100 %>%
  filter(Site == "TP", Depth == "15")

SJ_15 <- d100 %>%
  filter(Site == "SJ", Depth == "15")

SFTPSJ_15 <- rbind(SF_15,TP_15,SJ_15)

write.csv(SFTPSJ_15,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_15_100_coretaxa.csv")


SF_OM <- d100 %>%
  filter(Site == "SF", Depth == "OM")

TP_OM <- d100 %>%
  filter(Site == "TP", Depth == "OM")

SJ_OM <- d100 %>%
  filter(Site == "SJ", Depth == "OM")

SFTPSJ_OM <- rbind(SF_OM,TP_OM,SJ_OM)

write.csv(SFTPSJ_OM,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_OM_100_coretaxa.csv")











## 80% Cutoff - Looking for taxa that are present in at least 80% of a given site/depth combo in 100% of the sites
SF_5 <- d80 %>%
  filter(Site == "SF", Depth == "5")

TP_5 <- d80 %>%
  filter(Site == "TP", Depth == "5")

SJ_5 <- d80 %>%
  filter(Site == "SJ", Depth == "5")

SFTPSJ_5 <- rbind(SF_5,TP_5,SJ_5)

write.csv(SFTPSJ_5,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_5_coretaxa_80.csv")


SF_15 <- d80 %>%
  filter(Site == "SF", Depth == "15")

TP_15 <- d80 %>%
  filter(Site == "TP", Depth == "15")

SJ_15 <- d80 %>%
  filter(Site == "SJ", Depth == "15")

SFTPSJ_15 <- rbind(SF_15,TP_15,SJ_15)

write.csv(SFTPSJ_15,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_15_80_coretaxa.csv")


SF_OM <- d80 %>%
  filter(Site == "SF", Depth == "OM")

TP_OM <- d80 %>%
  filter(Site == "TP", Depth == "OM")

SJ_OM <- d80 %>%
  filter(Site == "SJ", Depth == "OM")

SFTPSJ_5 <- rbind(SF_OM,TP_OM,SJ_OM)

write.csv(SFTPSJ_OM,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_OM_80_coretaxa.csv")




## 70% Cutoff - Looking for taxa that are present in at least 70% of a given site/depth combo in 100% of the sites
SF_5 <- d70 %>%
  filter(Site == "SF", Depth == "5")

TP_5 <- d70 %>%
  filter(Site == "TP", Depth == "5")

SJ_5 <- d70 %>%
  filter(Site == "SJ", Depth == "5")

SFTPSJ_5 <- rbind(SF_5,TP_5,SJ_5)

write.csv(SFTPSJ_5,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_5_coretaxa_70.csv")


SF_15 <- d70 %>%
  filter(Site == "SF", Depth == "15")

TP_15 <- d70 %>%
  filter(Site == "TP", Depth == "15")

SJ_15 <- d70 %>%
  filter(Site == "SJ", Depth == "15")

SFTPSJ_5 <- rbind(SF_15,TP_15,SJ_15)

write.csv(SFTPSJ_15,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_15_70_coretaxa.csv")


SF_OM <- d70 %>%
  filter(Site == "SF", Depth == "OM")

TP_OM <- d70 %>%
  filter(Site == "TP", Depth == "OM")

SJ_OM <- d70 %>%
  filter(Site == "SJ", Depth == "OM")

SFTPSJ_OM <- rbind(SF_OM,TP_OM,SJ_OM)

write.csv(SFTPSJ_OM,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_OM_70_coretaxa.csv")




## 60% Cutoff - Looking for taxa that are present in at least 60% of a given site/depth combo in 100% of the sites
SF_5 <- d60 %>%
  filter(Site == "SF", Depth == "5")

TP_5 <- d60 %>%
  filter(Site == "TP", Depth == "5")

SJ_5 <- d60 %>%
  filter(Site == "SJ", Depth == "5")

SFTPSJ_5 <- rbind(SF_5,TP_5,SJ_5)

write.csv(SFTPSJ_5,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_5_coretaxa_60.csv")


SF_15 <- d60 %>%
  filter(Site == "SF", Depth == "15")

TP_15 <- d60 %>%
  filter(Site == "TP", Depth == "15")

SJ_15 <- d60 %>%
  filter(Site == "SJ", Depth == "15")

SFTPSJ_15 <- rbind(SF_15,TP_15,SJ_15)

write.csv(SFTPSJ_15,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_15_60_coretaxa.csv")


SF_OM <- d60 %>%
  filter(Site == "SF", Depth == "OM")

TP_OM <- d60 %>%
  filter(Site == "TP", Depth == "OM")

SJ_OM <- d60 %>%
  filter(Site == "SJ", Depth == "OM")

SFTPSJ_OM <- rbind(SF_OM,TP_OM,SJ_OM)

write.csv(SFTPSJ_OM,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_OM_60_coretaxa.csv")





## 50% Cutoff - Looking for taxa that are present in at least 50% of a given site/depth combo in 100% of the sites
SF_5 <- d50 %>%
  filter(Site == "SF", Depth == "5")

TP_5 <- d50 %>%
  filter(Site == "TP", Depth == "5")

SJ_5 <- d50 %>%
  filter(Site == "SJ", Depth == "5")

SFTPSJ_5 <- rbind(SF_5,TP_5,SJ_5)

write.csv(SFTPSJ_5,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_5_coretaxa_50.csv")


SF_15 <- d50 %>%
  filter(Site == "SF", Depth == "15")

TP_15 <- d50 %>%
  filter(Site == "TP", Depth == "15")

SJ_15 <- d50 %>%
  filter(Site == "SJ", Depth == "15")

SFTPSJ_15 <- rbind(SF_15,TP_15,SJ_15)

write.csv(SFTPSJ_15,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_15_50_coretaxa.csv")


SF_OM <- d50 %>%
  filter(Site == "SF", Depth == "OM")

TP_OM <- d50 %>%
  filter(Site == "TP", Depth == "OM")

SJ_OM <- d50 %>%
  filter(Site == "SJ", Depth == "OM")

SFTPSJ_OM <- rbind(SF_OM,TP_OM,SJ_OM)

write.csv(SFTPSJ_OM,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_OM_50_coretaxa.csv")






## 30% Cutoff - Looking for taxa that are present in at least 30% of a given site/depth combo in 100% of the sites
SF_5 <- d30 %>%
  filter(Site == "SF", Depth == "5")

TP_5 <- d30 %>%
  filter(Site == "TP", Depth == "5")

SJ_5 <- d30 %>%
  filter(Site == "SJ", Depth == "5")

SFTPSJ_5 <- rbind(SF_5,TP_5,SJ_5)

write.csv(SFTPSJ_5,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_5_coretaxa_30.csv")


SF_15 <- d30 %>%
  filter(Site == "SF", Depth == "15")

TP_15 <- d30 %>%
  filter(Site == "TP", Depth == "15")

SJ_15 <- d30 %>%
  filter(Site == "SJ", Depth == "15")

SFTPSJ_15 <- rbind(SF_15,TP_15,SJ_15)

write.csv(SFTPSJ_15,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_15_30_coretaxa.csv")


SF_OM <- d30 %>%
  filter(Site == "SF", Depth == "OM")

TP_OM <- d30 %>%
  filter(Site == "TP", Depth == "OM")

SJ_OM <- d30 %>%
  filter(Site == "SJ", Depth == "OM")

SFTPSJ_OM <- rbind(SF_OM,TP_OM,SJ_OM)

write.csv(SFTPSJ_OM,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/SFTPSJ_OM_30_coretaxa.csv")
