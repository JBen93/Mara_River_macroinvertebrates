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
set.seed(123)
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
orditorp(nmds_jaccard,display="sites",col=c(rep("black",4),rep("brown",4)),air=0.01,cex=1.25)
legend(-1,.5, c("Historical","Current"), cex=0.8,
       col=c("black","brown"), pch=15:15)
#Add  a dotted eclipse around the points to show the groups:



ordiellipse(nmds_jaccard, treat, display="si",lty=3,kind = "sd",conf = 0.75, col="black", show.groups="Historical")

ordiellipse(nmds_jaccard, treat, display="si", lty=3,kind = "sd",conf = 0.75, col="brown", show.groups="Current")

#add the stress value to the plot
text(-0.5,-0.5, paste("Stress = ", round(nmds_jaccard$stress,2)))


#####################################################################################################
#using PCoA to determine the temporal changes in the macroinvertebrates community structure in the Mara river.
#####################################################################################################



