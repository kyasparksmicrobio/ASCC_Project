# Read in your file
#Here, the skip flag tells R to read in the file while skipping over the first row
alpha2 <- read.csv("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_alpha_diversity/AlphaDiversity_ASCC_ITS_bothdepths_August2024_rarefied.csv", row.names=1,header=TRUE)

#Remove the 'taxaonomy' column if present
#data <- subset(data, select = -c(taxonomy))

# Calculate sum of counts for each sample (column)
sample_sums <- colSums(alpha2)

# Calculate relative abundance - which is the proportion of counts for a AVS out of the total for each sample
relative_abundance <- sweep(alpha2, 2, sample_sums, FUN="/") * 100

#Save file as .csv
write.csv(relative_abundance, "/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_rem_RELABUND_rarefied_RELABUND.csv")

