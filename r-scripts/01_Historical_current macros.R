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
print(jmacros)
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

#leged 
legend(-0.3,0.3, c("Historical","Current"), cex=0.8,
       col=c("black","brown"), pch=15:15)

#####################################################################################################
#PerMANOVA test to determine if there is a significant difference between the historical and current macroinvertebrate community structure in the Mara river.
#####################################################################################################
#next step is to determine if the historical and current sites are significantly different from each other using a permutation test.
#You will use the adonis function in the vegan package to run the permutation test:
set.seed(11) #set seed for reproducibility
permmacros<-adonis2(jmacros ~ treat, permutations=999)
permmacros

#The p-value is 0.024*, which means that the historical and current sites are significantly different from each other.
#Number of permutations: 999
#F-statistic: 15.369
#R-squared: 0.719


#####################################################################################################
#betadisper test 
#use the betadisper function in the vegan package :to test the differences betweeen the current sites
#betadisper test is used to determine the homogeneity of variance between the groups.
#####################################################################################################
set.seed(11) #set seed for reproducibility
#use the combined data-historical and current (presence, absence data)
# Euclidean distances between sites
dis<- vegdist(combinedmacros, "euclidean")
# groups are the 2 different years (historical and current)
groups <- as.factor(c(rep("Historical",4),rep("Current",4)))
#multivariate dispersions
MVdisp <- betadisper(dis, groups)
MVdisp

#Perform parametric test to determine if the multivariate dispersions are significantly different between the historical and current sites.
disp_aov<-anova(MVdisp)
disp_aov

# Tukey's Honest Significant Differences (HSD) test to determine which groups are significantly different from each other.
MVdisp.HSD <- TukeyHSD(MVdisp)
MVdisp.HSD

## non-parametric test: Permutation test for F-statistic to determine if the multivariate dispersions are significantly different between the historical and current sites.
perm_MVdisp <- permutest(MVdisp, permutations = 99, pairwise = TRUE)
perm_MVdisp


#####################################################################################################
#Simpler analysis to identify the taxa that contribute the most to the dissimilarity between the historical and current sites.
#####################################################################################################
#The next step is to determine which taxa are driving the differences between the historical and current sites.
#You will use the similarity percentage (SIMPER) analysis to identify the taxa that contribute the most to the dissimilarity between the historical and current sites.

# Perform SIMPER analysis to identify the taxa that contribute the most to the dissimilarity between the historical and current sites

names(combinedmacros)
group<- as.factor (c("Leptoceridae", "Hydropsychidae", "Ecnomidae", "Perlidae", "Libellulidae", "  Gomphidae", "Coenagrionidae", "Pyralidae", "Veliidae", "Notonectidae"," Nepidae","Corixidae","Belostomatidae"," Leptophlebiidae","Heptageniidae"," Bivalvia"," Naucoridae","Baetidae","Simuliidae","Muscidae","Ceratopogonidae"," Scirtidae"," Hydrophilidae","Elmidae","Dytiscidae","Ephemereliidae","Caenidae","Planaria"," Philopotamidae","Tabanidae"," Tipulidae","Tricorythidae"," Chironomidae","Corydalidae","Dicercomyzidae"," Gyrinidae","Lepidostomatidae","Naididae"))
group
# Run the SIMPER analysis
simper_result <- simper(combinedmacros, group, permutations = 999)
# View the results
print(simper_result)

