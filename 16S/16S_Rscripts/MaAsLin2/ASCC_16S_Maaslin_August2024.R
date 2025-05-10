# load libraries
##############
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Maaslin2")

library(Maaslin2)
library(devtools)
library(tidyverse) # tidyverse has it all! (ggplot2, dplyr, tidyr, etc...)
library(readxl)

otus <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_denoise_files/ASCC_16S_feature_table_rem_1.txt", row.names=1,header=TRUE)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_metadata_rem.txt",row.names=1) 

taxa <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_taxonomy/filtered/ASCC_16S_taxonomy_forprocessing.txt",row.names=1)


transposed_table <- t(otus)
transposed_table1 <- as.data.frame(transposed_table)

#write.csv(transposed_table1,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_feature_table_rem_transposed.csv")


# Extract the 'taxonomy' column from 'taxa'
new_colnames <- taxa[, "taxonomy"]

# Assign these values as the new column names for 'transposed_table1'
colnames(transposed_table1) <- as.character(new_colnames)

meta_shallow_subset <- subset(metadata, Depth == "0_5cm")
meta_deep_subset <- subset(metadata, Depth == "5_15cm")
meta_OM_subset <- subset(metadata, Depth == "OM")

#######

metadata_SF_TP <- subset(metadata, Site %in% c("SF","TP"))
metadata_TP_SJ <- subset(metadata, Site %in% c("TP","SJ"))
metadata_SJ_SF <- subset(metadata, Site %in% c("SJ","SF"))

metadata_SF_TP_shallow <- subset(meta_shallow_subset, Site %in% c("SF","TP"))
metadata_TP_SJ_shallow <- subset(meta_shallow_subset, Site %in% c("TP","SJ"))
metadata_SJ_SF_shallow <- subset(meta_shallow_subset, Site %in% c("SJ","SF"))

metadata_SF_TP_deep <- subset(meta_deep_subset, Site %in% c("SF","TP"))
metadata_TP_SJ_deep <- subset(meta_deep_subset, Site %in% c("TP","SJ"))
metadata_SJ_SF_deep <- subset(meta_deep_subset, Site %in% c("SJ","SF"))

metadata_SF_TP_OM <- subset(meta_OM_subset, Site %in% c("SF","TP"))
metadata_TP_SJ_OM <- subset(meta_OM_subset, Site %in% c("TP","SJ"))
metadata_SJ_SF_OM <- subset(meta_OM_subset, Site %in% c("SJ","SF"))

# metadata_SF_shallow_deep <- subset(meta_SF_sub, Depth %in% c("0_5cm","5_15cm"))
# metadata_TP_deep_OM <- subset(meta_TP_sub, Depth %in% c("5_15cm","OM"))
# metadata_SJ_OM_shallow <- subset(meta_SJ_sub, Depth %in% c("OM","0_5cm"))

# metadata_SF_deep <- subset(meta_deep_subset, Depth %in% c("0_5cm","5_15cm"))
# metadata_TP_deep<- subset(meta_deep_subset, Depth %in% c("5_15cm","OM"))
# metadata_SJ_deep <- subset(meta_deep_subset, Depth %in% c("OM","0_5cm"))

# metadata_SF_OM <- subset(meta_OM_subsetF, Depth %in% c("0_5cm","5_15cm"))
# metadata_TP_OM<- subset(meta_OM_subset, Depth %in% c("5_15cm","OM"))
# metadata_SJ_OM <- subset(meta_OM_subset, Depth %in% c("OM","0_5cm"))

# metadata_shallow_deep <- subset(metadata, Site %in% c("0_5cm","5_15cm"))
# metadata_deep_OM <- subset(metadata, Site %in% c("5_15cm","OM"))
# metadata_OM_shallow <- subset(metadata, Site %in% c("OM","0_5cm"))




#extracts the sample names from the metadataShallow dataframe
sample_names1 <- rownames(metadata_SF_TP)
sample_names2 <- rownames(metadata_TP_SJ)
sample_names3 <- rownames(metadata_SJ_SF)
sample_names4 <- rownames(metadata_SF_TP_shallow)
sample_names5 <- rownames(metadata_TP_SJ_shallow)
sample_names6 <- rownames(metadata_SJ_SF_shallow)
sample_names7 <- rownames(metadata_SF_TP_deep)
sample_names8 <- rownames(metadata_TP_SJ_deep)
sample_names9 <- rownames(metadata_SJ_SF_deep)
sample_names10 <- rownames(metadata_SF_TP_OM)
sample_names11 <- rownames(metadata_TP_SJ_OM)
sample_names12 <- rownames(metadata_SJ_SF_OM)


otus_SF_TP <- transposed_table1[sample_names1, ]
otus_TP_SJ <- transposed_table1[sample_names2, ]
otus_SJ_SF <- transposed_table1[sample_names3, ]
otus_SF_TP_shallow <- transposed_table1[sample_names4, ]
otus_TP_SJ_shallow <- transposed_table1[sample_names5, ]
otus_SJ_SF_shallow <- transposed_table1[sample_names6, ]
otus_SF_TP_deep <- transposed_table1[sample_names7, ]
otus_TP_SJ_deep <- transposed_table1[sample_names8, ]
otus_SJ_SF_deep <- transposed_table1[sample_names9, ]
otus_SF_TP_OM <- transposed_table1[sample_names10, ]
otus_TP_SJ_OM <- transposed_table1[sample_names11, ]
otus_SJ_SF_OM <- transposed_table1[sample_names12, ]



# write.csv(dataShallow, "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/Church Park/16S/Maaslin2 test/\\datashallow.csv", row.names=FALSE)
# write.csv(metadataShallow, "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/Church Park/16S/Maaslin2 test/\\metadatashallow.csv", row.names=FALSE)

#Now, dataCT should contain only the rows corresponding to the samples present in metadataCT, forming a new subsetted data matrix

#Run maaslin again but change 'data' and 'metadata' to whatever new matrices is named!
Maaslin2(input_data = otus_SF_TP,
          input_metadata = metadata_SF_TP,
          output = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SF_TP",
          fixed_effects = "Site",
          max_significance = 0.05)


Maaslin2(input_data = otus_TP_SJ,
          input_metadata = metadata_TP_SJ,
          output = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_TP_SJ",
          fixed_effects = "Site", 
          max_significance = 0.05)
 
Maaslin2(input_data = otus_SJ_SF,
          input_metadata = metadata_SJ_SF,
          output = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SJ_SF",
          fixed_effects = "Site", 
          max_significance = 0.05)

Maaslin2(input_data = otus_SF_TP_shallow,
          input_metadata = metadata_SF_TP_shallow,
          output = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SF_TP_shallow",
          fixed_effects = "Site", 
          max_significance = 0.05)
 
Maaslin2(input_data = otus_TP_SJ_shallow,
          input_metadata = metadata_TP_SJ_shallow,
          output = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_TP_SJ_shallow",
          fixed_effects = "Site", 
          max_significance = 0.05)

Maaslin2(input_data = otus_SJ_SF_shallow,
          input_metadata = metadata_SJ_SF_shallow,
          output = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SJ_SF_shallow",
          fixed_effects = "Site", 
          max_significance = 0.05)

Maaslin2(input_data = otus_SF_TP_deep,
          input_metadata = metadata_SF_TP_deep,
          output = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SF_TP_deep",
          fixed_effects = "Site", 
          max_significance = 0.05)
 
Maaslin2(input_data = otus_TP_SJ_deep,
          input_metadata = metadata_TP_SJ_deep,
          output = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_TP_SJ_deep",
          fixed_effects = "Site", 
          max_significance = 0.05)

Maaslin2(input_data = otus_SJ_SF_deep,
          input_metadata = metadata_SJ_SF_deep,
          output = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SJ_SF_deep",
          fixed_effects = "Site", 
          max_significance = 0.05)

Maaslin2(input_data = otus_SF_TP_OM,
          input_metadata = metadata_SF_TP_OM,
          output = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SF_TP_OM",
          fixed_effects = "Site", 
          max_significance = 0.05)

Maaslin2(input_data = otus_TP_SJ_OM,
          input_metadata = metadata_TP_SJ_OM,
          output = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_TP_SJ_OM",
          fixed_effects = "Site", 
          max_significance = 0.05)

Maaslin2(input_data = otus_SJ_SF_OM,
          input_metadata = metadata_SJ_SF_OM,
          output = "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SJ_SF_OM",
          fixed_effects = "Site", 
          max_significance = 0.05)

#########
# Graph results
#########
# read significant table in

data.sftp <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SF_TP/significant_results.tsv")
data.tpsj <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_TP_SJ/significant_results.tsv")
data.sjsf <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SJ_SF/significant_results.tsv")
data.sftp.shal <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SF_TP_shallow/significant_results.tsv")
data.tpsj.shal <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_TP_SJ_shallow/significant_results.tsv")
data.sjsf.shal <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SJ_SF_shallow/significant_results.tsv")
data.sftp.deep <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SF_TP_deep/significant_results.tsv")
data.tpsj.deep <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_TP_SJ_deep/significant_results.tsv")
data.sjsf.deep <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SJ_SF_deep/significant_results.tsv")
data.sftp.OM <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SF_TP_OM/significant_results.tsv")
data.tpsj.OM <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_TP_SJ_OM/significant_results.tsv")
data.sjsf.OM <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SJ_SF_OM/significant_results.tsv")


# basic plot
plot <- ggplot(data.sftp, aes(x = reorder(feature, -coef), y = coef)) + geom_bar(stat = "identity")
plot

data.sftp <- readxl::read_excel("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_Maaslin/output_16S_SF_TPoutput_shallow/significant_taxaonly.xlsx", 
                             sheet = "Sheet10")

data.1$color <- ifelse(data.1$coef_lda >= 0, "Positive", "Negative")

plot2 <- ggplot(data.1, aes(x = coef_lda, y = feature, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Positive" = "green", "Negative" = "red")) +
  theme_minimal()
plot2

ggsave("church_park_lda.pdf",plot = plot2, width=15,height = 20,dpi = 300)


