library(ggplot2)

# reading in files # 
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/csu/projects/hydrogen_extraction/16S/Round1_Feb2023/stacked_barcharts")
data<-read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_FUNguild/funguild_ASCC_ITS3.txt",header=T, check.names=FALSE)

# plotting
ggplot(data, aes(x=guild, y=relabund, fill=guiid)) +
  geom_bar(stat="identity") +
  theme_bw() + 
  #theme(legend.position = "none") +
  #facet_wrap(~basin, scales="free_x") +
  scale_fill_viridis_d(option = "plasma") +
  guides(fill = guide_legend(ncol = 1)) +
  xlab("Sample") +
  ylab("% Relative Abundance")
#theme(axis.text.x=element_text(angle = 45, hjust = 1))

ggsave("All_barchart.pdf")


ggplot(data, aes(x=Site, y=counts, fill=guild)) +
  geom_bar(stat="identity") +
  #facet_wrap(~Treatment, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
 
ggsave("facet_by_trt.pdf")
 




library(tidyverse)

library(readxl)
data <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG_rem.txt") %>%
  pivot_longer(names_to = "samples", values_to = "counts", 3:50)

write.csv(data,"/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_denoise_files/ASCC_ITS_feature_table_OG-TEST.csv")



taxa2 <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_FUNguild/ASCC_FUNGUILD_SIMPLE.txt") %>% 
  select(1:10)

metadata <- read.delim("/Users/kyasparks/Desktop/CSU 2023-2024/Projects/ASCC/ASCC_July2024_USETHIS_ks/ITS/ITS_metadata/ASCC_ITS_metadata_rem.txt") 


data_meta <- inner_join(data, metadata, by =c("samples"="samples")) %>% 
  pivot_wider(names_from = samples, values_from = counts)


ggplot(data_meta, aes(x=Site, y=counts, fill=Family)) +
  geom_bar(stat="identity") +
  facet_wrap(~Site, nrow = 1, scales = "free_x")+
  guides(fill=guide_legend(ncol = 1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
