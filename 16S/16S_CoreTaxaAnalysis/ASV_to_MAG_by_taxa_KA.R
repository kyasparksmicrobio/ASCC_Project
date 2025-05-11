# Load packages
library(tidyverse)
library(ggplot2)

setwd("/Users/kkamundson/Library/Mobile Documents/com~apple~CloudDocs/Documents/csu/projects/helping_others/kya/16S_matching_code")

#### READ IN FILES ####
# read in MAG table - this will always be the same (unless updated in the future by me)
mags <- read.delim("coniferous_MAGdb_1172MAGs.txt", header = TRUE)

# read in user feature table 
taxa <- read.delim("taxa.txt", header = TRUE)


################ FORMAT MAG TABLE ##############################################
# This is truncating the MAG tax string at different levels for a 'lookup' later on 
# Technically only needs to be done once per environment load 
# GENUS - truncate taxonomic string for genus level 
mags <- mags %>%
  mutate(MAG_genus = sub("(;s__.*)$", "", MAG_FULL_tax))

# FAMILY - truncate the taxonomy string 
mags <- mags %>%
  mutate(MAG_family = sub("(;g__.*)$", "", MAG_FULL_tax))

# ORDER - truncate the taxonomy string
mags <- mags %>%
  mutate(MAG_order = sub("(;f__.*)$", "", MAG_FULL_tax))

# CLASS - truncate the taxonomy string
mags <- mags %>%
  mutate(MAG_class = sub("(;o__.*)$", "", MAG_FULL_tax))




################### FORMAT ASV TABLE ###########################################
taxa <- as.data.frame(taxa)
colnames(taxa)[2] <- "ASV_FULL_tax"

taxa <- taxa %>%
  separate(ASV_FULL_tax, into = c("d","p","c","o","f","g","s"), sep = ";", remove = FALSE)

# GENUS - truncate the taxonomy string to everything after 'genus' by removing all species
taxa <- taxa %>%
  mutate(ASV_genus = sub("(;s__.*)$", "", ASV_FULL_tax))

# FAMILY - truncate the taxonomy string to family level
taxa <- taxa %>%
  mutate(ASV_family = sub("(;g__.*)$", "", ASV_FULL_tax))

# ORDER - truncate the taxonomy string to order level
taxa <- taxa %>%
  mutate(ASV_order = sub("(;f__.*)$", "", ASV_FULL_tax))

# CLASS - truncate the taxonomy string to class level
taxa <- taxa %>%
  mutate(ASV_class = sub("(;o__.*)$", "", ASV_FULL_tax))



################################################################################
######################## START OF COMPARING ####################################
################################################################################


########## COUNTING ASVs THAT ARE INSUFFICIENTLY CLASSIFIED ####################
# Counting the number of ASVs that are Unassigned at the domain (d) level
insuf_ASV_UNASSIGNED <- sum(taxa$d == "Unassigned", na.rm = TRUE)

# Counting the number of ASVs that do not classify past the domain level 
insuf_ASV_justDOMAIN <- sum(is.na(taxa$p))

# Counting the number of ASVs that do not classify past the phylum level 
insuf_ASV_justPHYLUM <- sum(is.na(taxa$c))

# Counting the number of ASVs that do not classify past the class level 
insuf_ASV_justCLASS <- sum(is.na(taxa$o))

# Make a new column 'match_level' and mark all of the above as having insufficient data to continue forward 
merged_data <- taxa %>%
  mutate(match_level = ifelse(is.na(o), "insufficient ASV classification", NA))

# Count the number of ASVs with insufficient classification to move forward
count_ASV_insuf <- sum(merged_data$match_level == "insufficient ASV classification", na.rm = TRUE)



#################### FULL TAX MATCH ############################################
# Assuming 'taxa' is the dataframe with the 'ASV_FULL_tax' column
# Assuming 'mags' is the dataframe with the 'MAG_FULL_tax', 'MAG', and 'compl' columns

# Select rows with the highest completion value for each 'MAG_FULL_tax'
# Need to do this as there might be more than one MAG with identical taxanomic classification - so in those cases it should choose the highest completion MAG
mags_filtered_FULL <- mags %>%
  group_by(MAG_FULL_tax) %>%
  slice_max(order_by = Completeness, with_ties = FALSE) %>%
  ungroup()

# Count the number of unique MAGs at full classification in dataset
MAGs_unique_FULL <- nrow(mags_filtered_FULL)

# Filter rows with NA in merged_data 
na_rows_full <- merged_data %>% filter(is.na(match_level))

# left join to look up MAG for matching full classificaiton between ASV and MAG
merged_data_full <- left_join(na_rows_full, mags_filtered_FULL %>% select(MAG_FULL_tax, MAG), by = c("ASV_FULL_tax" = "MAG_FULL_tax"))

# Add a new column based on condition
merged_data_full <- merged_data_full %>%
  mutate(match_level = ifelse(!is.na(MAG), "full tax", NA))

# Combine with original dataset
merged_data <- bind_rows(merged_data %>% filter(!is.na(match_level)), merged_data_full)

# CHECK SOME NUMBERS !
## First - merged_data should be the same as 'feat' 
nrow(merged_data)

# Second - how many ASVs match at the full taxanomic string match level?
count_ASV_FULL <- sum(merged_data$match_level == "full tax", na.rm = TRUE)



#################### GENUS LEVEL MATCH #########################################
# Assuming 'taxa' is the dataframe with the 'ASV_genus' column
# Assuming 'mags' is the dataframe with the 'MAG_genus' column
# Assuming 'merged_data' is the dataframe from the previous step

# Filter mags to just unique with highest completion 
mags_filtered_GENUS <- mags %>%
  group_by(MAG_genus) %>%
  slice_max(order_by = Completeness, with_ties = FALSE) %>%
  ungroup()

# Count the number of unique MAGs at genus level classification in dataset
MAGs_unique_GENUS <- nrow(mags_filtered_GENUS)

# Filter rows with NA in taxa match (MAG/level match columns)
na_rows_genus <- merged_data %>% filter(is.na(match_level))

# Left join on filtered rows
merged_data_genus <- left_join(na_rows_genus, mags_filtered_GENUS %>% select(MAG_genus, MAG), by = c("ASV_genus" = "MAG_genus"))

# Need to rename the column due to the default naming
merged_data_genus <- merged_data_genus %>%
  rename(MAG = MAG.y)

# Fill in new_column with 'genus' for newly added rows
merged_data_genus <- merged_data_genus %>%
  mutate(match_level = ifelse(!is.na(MAG), "genus", NA))

# Combine original merged_data with new_merged_data
merged_data <- bind_rows(merged_data %>% filter(!is.na(match_level)), merged_data_genus)

# Need to remove random MAG.x column it makes
merged_data <- merged_data %>%
  select(-MAG.x)

# Again - CHECK some numbers!
# check merged data
nrow(merged_data)

# Count how many ASVs matched at the genus level
count_ASV_GENUS <- sum(merged_data$match_level == "genus", na.rm = TRUE)




########################### FAMILY LEVEL #######################################
# Assuming 'taxa' is your dataframe with the 'ASV_genus' column
# Assuming 'mags' is your dataframe with the 'MAG_genus' column
# Assuming 'merged_data' is your dataframe from the previous step
# Filter mags to just unique with highest completion 

# Filter mags to just unique with highest completion at family level
mags_filtered_FAMILY <- mags %>%
  group_by(MAG_family) %>%
  slice_max(order_by = Completeness, with_ties = FALSE) %>%
  ungroup()

# Count the number of unique MAGs at family level classification in dataset
MAGs_unique_FAMILY <- nrow(mags_filtered_FAMILY)

# Filter rows with NA in taxa match (MAG/level match columns)
na_rows_family <- merged_data %>% filter(is.na(match_level))

# Left join on filtered rows
merged_data_family <- left_join(na_rows_family, mags_filtered_FAMILY %>% select(MAG_family, MAG), by = c("ASV_family" = "MAG_family"))

# Need to rename the column due to the default naming
merged_data_family <- merged_data_family %>%
  rename(MAG = MAG.y)

# Fill in new_column with 'family' for newly added rows
merged_data_family <- merged_data_family %>%
  mutate(match_level = ifelse(!is.na(MAG), "family", NA))

# Combine original merged_data with new_merged_data
merged_data <- bind_rows(merged_data %>% filter(!is.na(match_level)), merged_data_family)

# Need to remove random MAG.x column it makes
merged_data <- merged_data %>%
  select(-MAG.x)

# Again - CHECK some numbers!
# check merged data
nrow(merged_data)

# Count how many ASVs matched at the genus level
count_ASV_FAMILY <- sum(merged_data$match_level == "family", na.rm = TRUE)



######################### ORDER LEVEL ##########################################
# Assuming 'taxa' is the dataframe with the 'ASV_order' column
# Assuming 'mags' is the dataframe with the 'MAG_order' column
# Assuming 'merged_data' is the dataframe from the previous step
# Filter mags to just unique with highest completion 

# Filter mags to just unique with highest completion at order level
mags_filtered_ORDER <- mags %>%
  group_by(MAG_order) %>%
  slice_max(order_by = Completeness, with_ties = FALSE) %>%
  ungroup()

# Count the number of unique MAGs at order level classification in dataset
MAGs_unique_ORDER <- nrow(mags_filtered_ORDER)

# Filter rows with NA in taxa match (MAG/level match columns)
na_rows_order <- merged_data %>% filter(is.na(match_level))

# Left join on filtered rows
merged_data_order <- left_join(na_rows_order, mags_filtered_ORDER %>% select(MAG_order, MAG), by = c("ASV_order" = "MAG_order"))

# Need to rename the column due to the default naming
merged_data_order <- merged_data_order %>%
  rename(MAG = MAG.y)

# Fill in new_column with 'family' for newly added rows
merged_data_order <- merged_data_order %>%
  mutate(match_level = ifelse(!is.na(MAG), "order", NA))

# Combine original merged_data with new_merged_data
merged_data <- bind_rows(merged_data %>% filter(!is.na(match_level)), merged_data_order)

# Need to remove random MAG.x column it makes
merged_data <- merged_data %>%
  select(-MAG.x)

# Again - CHECK some numbers!
# check merged data
nrow(merged_data)

# Count how many ASVs matched at the genus level
count_ASV_ORDER <- sum(merged_data$match_level == "order", na.rm = TRUE)


################################################################################
################################################################################
################################################################################


###############################################################################
###################### FINAL SAVED OUTPUT TABLES ###############################

# First need to collect the MAG data and add it to the ASV (merged_data) matrix
merged_data_FINAL <- merged_data %>%
  left_join(
    mags %>% select(MAG, fire, year_post_burn, severity, Completeness, Contamination, MAG_FULL_tax, MAG_size_bp, numb_contigs),
    by = "MAG")

# Need to establish an order to the data so that it presents in a logical order downstream
match_order <- c("full tax", "genus", "family", "order", "insufficient ASV classification", "NA")

# Clean up final merged data matrix for user output table
merged_data_FINAL <- merged_data_FINAL %>%
  arrange(match_level)

remove_columns <- c("d", "p", "c", "o", "f", "g", "s", "ASV_genus", "ASV_family", "ASV_order", "ASV_class")

merged_data_OUTPUT <- merged_data_FINAL %>%
  select(-all_of(remove_columns))

##### SAVE FINAL DATA FRAME FOR USER ####
write_csv(merged_data_OUTPUT, "ASV_to_MAGs_output.csv")


# Summarize the data to get counts of each match_level
match_level_counts <- merged_data_FINAL %>%
  group_by(match_level) %>%
  summarize(count = n()) %>%
  ungroup()

# Also want to add percentages to this - setting up for ggplot pie chart etc below
match_level_counts <- match_level_counts %>%
  mutate(percentage = count / sum(count) * 100) %>%
  mutate(percentage = round(count / sum(count) * 100, 2)) %>%
  mutate(label = paste0(match_level, "\n(n=", count, ", ", percentage, "%)"))

##### SAVE COUNTS AS OUTPUT ####
write.csv(match_level_counts, "ASV_matching_percents.csv")












