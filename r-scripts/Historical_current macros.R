# Comapr

# clear everything in memory (of R)
remove(list=ls())

renv::restore()

library(vegan) 
library(psych) 
library(tidyverse)

# read the vegetation data
combinedmacros<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1573424697&single=true&output=csv") |>
#names as variable to use
tibble::column_to_rownames(var="Location_ID") # convert distance_m to the row names of the tibble

#The NMDS plot will be used to determine the temporal changes in the macroinvertebrates community structure in the Mara river.

#NMDS plot of the macros_PCoA data
nmdsmacros<-metaMDS(combinedmacros,k=2, trace=T)

#For plotting the NMDS, create the groups to assign different colors to each time period. You will also use these groups for the ANOSIM:
group=as.matrix(c(rep("Historical",4),rep("Current",4)))
ordiplot(nmdsmacros,type="n",xlim=c(-0.5,1),ylim=c(-.5,.5))
## species scores not available for this method
orditorp(nmdsmacros,display="sites",col=c(rep("green",4),rep("blue",4)),air=0.01,cex=1.25)
# Add the legend on top right
legend("topleft",legend=c("Historical","Current"),fill=c("green","blue"),bty="n",cex=0.8)

# Add the title
title("NMDS plot")

