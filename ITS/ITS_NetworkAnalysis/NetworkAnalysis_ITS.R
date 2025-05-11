# install.packages("igraph") # uncomment this line in order to install this package
library(igraph)  
# install.packages("Hmisc") # uncomment this line in order to install this package
library(Hmisc)  
# install.packages("Matrix") # uncomment this line in order to install this package
library(Matrix) 
library(tidyverse)
#devtools::install_github("RMHogervorst/gephi")
library(gephi)
To download gephi

asv.table<-read.csv("networks/genus_table_filtered_p8.csv", header=T, row.names = 1)
tax<-read.csv("networks/genus_id_gtdb.csv",header=T, row.names = 1)

##Check the dimensions of the data frame

dim(asv.table)

##Filter out low abundance asvs (less than 10 reads)

asv.table.filter <- asv.table[ ,colSums(asv.table) >= 10]
print(c(ncol(asv.table),"versus",ncol(asv.table.filter))) #compare initial and filtered counts
We use rcorr to calculate the Spearman correlation coefficient between ASVs. This creates a list with 3 values:

# "r" (rho): Correlation coefficient

# "n": number of observations

# "P": p-value

asv.cor <- rcorr(as.matrix(asv.table.filter), type="spearman")
asv.cor$P will pull the p-value info from the list above and forceSymmetric() will assign NA to self-correlations

asv.pval <- forceSymmetric(asv.cor$P) 

##Filter taxa to select only the ASVs retained after filtering for low read counts

sel.tax <- tax[rownames(asv.cor$P),,drop=FALSE]


##Make sure your filtered tables match

all.equal(rownames(sel.tax), rownames(asv.pval))

##Filter to retain only significant associations

p.yes <- asv.pval<0.05
r.val = asv.cor$r # select all the correlation values
p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion 


##Select ASVs based on Spearman Correlation. Here we are keeping only values with correlation coefficients higher than 0.75, but this value can be adjusted to retain stronger or weaker correlations.

p.yes.r <- abs(p.yes.r)>0.75 # output is logical vector
p.yes.rr <- p.yes.r*r.val # use logical vector for subscripting.


##Create an adjacency matrix

adjm <- as.matrix(p.yes.rr)


##The next step creates an object from the adjacency matrix. Weight represents the level of correlation

net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE)

##You can use this to pull out edge weights

edgew<-E(net.grph)$weight

##Creating a vector to remove the isolated nodes (nodes with no interactions) and then remove those nodes from the object

bad.vs<-V(net.grph)[degree(net.grph) == 0] 

net.grph <-delete.vertices(net.grph, bad.vs)

##Write out the edge files

gephi_write_edges(net.grph, "networks/p8_filtered_0.75_edge_all_edge.csv")

##Create a nodes file

taxa<-read_csv("networks/genus_id_gtdb.csv")
edges <- read_csv("networks/p8_filtered_0.75_edge_all_edge.csv")
edges_s <- edges %>% 
  select(1)
edges_t <- edges %>% 
  select(2)
edges_all <- edges_s %>% 
  full_join(edges_t,by = c("Source" = "Target")) %>% 
  distinct() %>% 
  rename("GenID" = "Source") %>% 
  left_join(taxa)