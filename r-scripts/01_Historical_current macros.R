# Compare the historical and current macroinvertebrate community structure in the Mara River using NMDS, DCA, and PCoA to determine the temporal changes in the community structure.

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
set.seed(123)
jmacros <- vegdist(combinedmacros, method = "jaccard")

#You are going to use the cmdscale function in the stats package to run the PCoA:
cmd<-cmdscale(jmacros, k=5, eig=TRUE) 
cmd

#The “points” are the coordinates of each sites. They are the eigenvectors scaled by the square root of their eigenvalues (i.e. the standard deviation):
cmd$points

#make a PCoA table to look at the eigenvalues, and the proportional and cumulative variance:
eigenvalues<-cmd$eig[1:5]
propVar<-eigenvalues/sum(eigenvalues)
cumVar<-cumsum(propVar)
PCoA_Table<-cbind(eigenvalues,propVar,cumVar)
PCoA_Table

#How many axes should you keep? You can use a scree plot to help you decide.
#Scree plot:
plot(eigenvalues)
lines(lowess(eigenvalues))

#The scree plot shows that the first two axes are the most important.

#Now, plot the first two PCoA axes:
x<-cmd$points[,1]
y<-cmd$points[,2]
plot(x,y,xlab= "PCoA1(72%)", ylab="PCoA2 (13.5%)", xlim=range(x)*1.2,ylim=range(y)*1.2, type="n")
text(x,y,labels=rownames(cmd$points), cex=.9)
# Add the title:
title(main="Macroinvertebrate taxa composition", cex.main=1)

#add dotted lines to the inside of the plot to show the axes:
#abline(v=0, lty=3)
#abline(h=0, lty=3)

#use different colors for the historical and current data points using ecllipes
ordiellipse(cmd, treat, display="sites",lty=3,kind = "sd",conf = 0.75, col="black", show.groups="Historical")
ordiellipse(cmd, treat, display="sites", lty=3,kind = "sd",conf = 0.60, col="brown", show.groups="Current")



