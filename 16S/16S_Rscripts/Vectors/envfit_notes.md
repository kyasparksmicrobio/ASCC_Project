#### R Code for Envfit - Adding Vectors to NMDS

## Generate NMDS first, then run the following code to add vectors

```r
# First, some data checks
env_data = as.data.frame(sample_tab) # Assigning your data (sample_tab file with metadata) a new variable name 
all.equal(row.names(sample_tab), row.names(env_data)) # Just checking all is aligned

# Run envfit
fit = envfit(mgd_relabund.bray.nmds, env_data, perm = 10000) # Function envfit finds vectors or factor averages of environmental variables/metadata
fit.scrs = as.data.frame(scores(fit, display = "vectors")) # Add those scores as vectors to the plot
fit.scrs = cbind(fit.scrs, Species = rownames(fit.scrs)) # Adding this info to the matrix 'fit.scrs'

# This is your ggplot code from your normal NMDS - adjust as needed so yours looks similar
# Note here we are assigning the plot as 'p' - so it will not show the plot until you run line 20 below which is just 'p'
p = ggplot(sample_tab, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(size = day, color = basin)) +
  theme_classic() +
  labs(size = "Day", col = "Basin")   
p 

# ADD VECTORS TO THE PLOT
# Now we are telling ggplot to add more 'stuff' to the plot - the envfit information 
# So we say our 'p' is the old plot + new stuff (everything after geom_segment)  
# Again, we run 'p' to see what this has added, hopefully vectors! 
p = p +
  geom_segment(data = fit.scrs, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), arrow.fill = "black")
p

# LABEL VECTORS
# Again, we are building off the previous plots (as saved by 'p')
# This time the ordination should be updated to have vectors that are labeled 
p = p + geom_text(data = fit.scrs, aes(x = NMDS1, y = NMDS2, label = Species), size = 3.5) 
p 

# Save the ordination + vectors 
ggsave("gross_NMDS-with-vectors.pdf")

# The following commands will list which variables have a significant p value (first line), 
# show the p-values (second line), and then show all data (third line)
which(fit$vectors$pvals <= 0.05)
fit$vectors$pvals
fit$vectors

# Write the p-values for the variables to a .csv file and save it in your working directory
write.table(fit$vectors$pvals, file = 'pvals_vectors_sulfur.txt', quote = FALSE, sep = '/t')
```
