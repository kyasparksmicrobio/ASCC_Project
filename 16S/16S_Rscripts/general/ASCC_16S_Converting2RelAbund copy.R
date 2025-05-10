# Read in your file
#Here, the skip flag tells R to read in the file while skipping over the first row
otus <- read.csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_alpha_diversity/AlphaDiversity_ASCC_16S_bothdepths_July2024_rarefied.csv",row.names=1,header=TRUE)

#Remove the 'taxaonomy' column if present
#data <- subset(data, select = -c(taxonomy))

# Calculate sum of counts for each sample (column)
sample_sums <- colSums(otus)

# Calculate relative abundance - which is the proportion of counts for a AVS out of the total for each sample
relative_abundance <- sweep(otus, 2, sample_sums, FUN="/") * 100

#Save file as .csv
write.csv(relative_abundance, "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/16S/16S_metadata/ASCC_16S_feature_table_rarified_rem_RELABUND.csv")

