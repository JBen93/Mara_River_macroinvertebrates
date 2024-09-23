# Comapr

# clear everything in memory (of R)
remove(list=ls())

renv::restore()

library(vegan) 
library(psych) 
library(tidyverse)

#database source
##browseURL("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/edit?usp=sharing")
# read the vegetation data
combinedmacros<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1573424697&single=true&output=csv") |> 
tibble::column_to_rownames(var="Location_ID") # #names as variable to use in nmds otherwise it will not work, only works with numbers. 

#The NMDS plot will be used to determine the temporal changes in the macroinvertebrates community structure in the Mara river.

# Calculate the Jaccard dissimilarity matrix since it's a presence/absence data
jmacros <- vegdist(combinedmacros, method = "jaccard")

# Perform NMDS using the Jaccard dissimilarity matrix
nmds_jaccard <- metaMDS(jmacros, k = 2, trace =T)
stressplot(nmds_jaccard)
#plot out our results and see if there is a difference between the historical and current
#Identify the time period as groups:
treat=as.matrix(c(rep("Historical",4),rep("Current",4)))

#Plot out the points (Location_ID):
 
ordiplot(nmds_jaccard,type="n",xlim=c(-.5,.5),ylim=c(-.5,.5))

## species scores not available
orditorp(nmds_jaccard,display="sites",col=c(rep("green",4),rep("blue",4)),air=0.01,cex=1.25)
legend(-.55,.5, c("Historical","Current"), cex=0.8,
       col=c("green","blue"), pch=15:15)
#Add a convex hull around each group:

ordihull(nmds_jaccard, treat, display="si",lty=1, col="green", show.groups="Historical")

ordihull(nmds_jaccard, treat, display="si", lty=1, col="blue", show.groups="Current")





#For plotting the NMDS, create the groups to assign different colors to each time period. You will also use these groups for the ANOSIM:
# Adjust margins (bottom, left, top, right)



