#### Ternary Plots w/ Core Microbiome Info ####
library(tidyverse)
library(ggtern)
library(ggplot2)
## will need to restart R and get rid of ggthemr to run this correctly
## Colors for Core: “#DB735C”, “#E58E65", “#EFA86E”, “#F3C57B”, “#D5AE71”, “#B79667", “#997F5D”, “#7A6752"
## Colors for none: “#CDCBCB”
## Used the website Coolors to make the color palette wider
## Issue: Unburned? Right now I have unburned included in the percentages along the side
## Make sure to note this when telling Mike, also write down my methods
Data <- read.delim(“Input_Data_GGTern.txt”,header = T,check.names=FALSE)

# Data_Actinobacteria <- read.delim(“Input_Data_GGTern.txt”,header = T,check.names=FALSE) %>%
#   filter(Class==“Actinobacteria”)
# Data_Gammaproteobacteria <- read.delim(“Input_Data_GGTern.txt”,header = T,check.names=FALSE) %>%
#   filter(Class==“Gammaproteobacteria”)
# Data_Alphaproteobacteria <- read.delim(“Input_Data_GGTern.txt”,header = T,check.names=FALSE) %>%
#   filter(Class==“Alphaproteobacteria”)
# Data_Bacteroidia <- read.delim(“Input_Data_GGTern.txt”,header = T,check.names=FALSE) %>%
#   filter(Class==“Bacteroidia”)
# Data_Blastocatellia <- read.delim(“Input_Data_GGTern.txt”,header = T,check.names=FALSE) %>%
#   filter(Class==“Blastocatellia”)
# Data_Bacilli <- read.delim(“Input_Data_GGTern.txt”,header = T,check.names=FALSE) %>%
#   filter(Class==“Bacilli”)
# Data_Verrucomicrobiae <- read.delim(“Input_Data_GGTern.txt”,header = T,check.names=FALSE) %>%
#   filter(Class==“Verrucomicrobiae”)
# Data_Acidimicrobiia <- read.delim(“Input_Data_GGTern.txt”,header = T,check.names=FALSE) %>%
#   filter(Class==“Acidimicrobiia”)
# Data_Planctomycetes <- read.delim(“Input_Data_GGTern.txt”,header = T,check.names=FALSE) %>%
#   filter(Class==“Planctomycetes”)
unique(Data$Core)
Data[“Core”][Data[“Core”] == “Unburned Mineral “] <- “Unburned Mineral”
Data[“Core”][Data[“Core”] == “Organic & Mineral “] <- “Organic & Mineral”
Data[“Core”][Data[“Core”] == “Organic  “] <- “Organic”
Data[“Core”][Data[“Core”] == “Organic “] <- “Organic”
Data[“Core”][Data[“Core”] == “Ash & Organic “] <- “Ash & Organic”
unique(Data$Core)
Data$Core <- factor(Data$Core, levels = c(“Ash”, “Organic”, “Mineral”, “Ash & Organic”,
                                          “Ash, Organic, & Mineral”, “Ash & Organic; Unburned Mineral”,
                                          “Organic & Mineral”, “Unburned Mineral”, “None”))
## 8 x 13 in
Plot_All <- ggtern(data = Data, aes(x = Ash, y = Charred_Org_C, z = Mineral, color=Core)) +
  geom_point(aes(size=(Relative_Abundance))) +
  scale_size(name=“Mean Relative Abundance”)+ ##?
  theme_classic(base_size = 20) +
  facet_wrap(.~Class) +
  scale_color_manual(values = c(“#DB735C”, “#E58E65”, “#EFA86E”, “#F3C57B”,
                                “#D5AE71", “#B79667”, “#997F5D”, “#7A6752”, “#CDCBCB”)) +
  theme(legend.text=element_text(size=16),
        legend.title=element_text(size=16, face=“bold”),
        strip.text.x = element_text(size = 13, face = “bold”)) +
  Tlab(“Organic”) + Llab(“Ash”) + Rlab(“Min.“) +
  theme(axis.text = element_text(size= 10)) +
  theme(axis.title = element_text(size = 12)) +
  theme_nomask() +
  theme(legend.position=“right”) +
  scale_size_continuous(breaks=c(0.001,0.005,0.01,0.02,0.03,0.04,0.05,0.06))
Plot_All
#### Accompanying Venn Diagram ####
install.packages(“ggvenn”)
library(“ggvenn”)
Data_VennDiagram <- read.delim(“Input_Data_GGVenn.txt”,header = T,check.names=FALSE)
ggvenn(Data_VennDiagram, c(“Ash”, “Organic”,“Mineral”,“Unburned Mineral”),
       fill_color = c(“#DB735C”, “#EFA86E”, “#B79667”, “#7A6752"),
       text_size = 4,
       set_name_size = 5.5) +
  scale_x_continuous(expand = expansion(mult = 0.5))