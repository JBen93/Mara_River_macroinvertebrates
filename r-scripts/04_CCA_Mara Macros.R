# clear everything in memory (of R)
remove(list=ls())
setwd("/Users/Joshua/Library/CloudStorage/OneDrive-UniversityofFlorida/R input")
#call in biostats data from your working folder
source("biostats.r")
renv::restore()

library(vegan) 
library(psych) 
library(tidyverse)
library(ggpubr)
library(terra)
library(raster)
library(ade4)
library(ca)
library(VennDiagram)

#database source 
#browseURL("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/edit?usp=sharing")
#call in data cca_ macros data from the Google link for the macros abundance data
maramacros<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=510356961&single=true&output=csv")|> 
tibble::column_to_rownames(var="Site")#use the first column as the row names

#call in the cca_environmental data from the Google link
maraenv<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1414987851&single=true&output=csv")|> 
tibble::column_to_rownames(var="Site")#use the first column as the row names

#check the structure of this data set
str(maramacros)
str(maraenv)
summary(maramacros)
summary(maraenv)

##Data selection, transformation and standardization
#Species within community data sets vary greatly in there occurrence, abundance, and habitat specificity. 
#Species that are common, widespread and extremely abundant can obscure patterns in the ordination. 
#Species that are rare and have few occurrences in a data set may not be accurately placed in ecological space. 
#You must decide which species are “rare” and which are super abundant.
##Selecting Species##
#To explore patterns of rarity and commonness, you will use the foa function from the Biostats package. 
#This function will give you a whole series of plots that allow you to explore the occurrence and abundance patterns of the species in your data. 
#The second plot, Empirical Distribution of Species Relative Occurrence, will be the one we use to remove common and/or rare species. 

occur<-foa.plots(maramacros)
occur

#determine the rare species are those that occur in less than 5 samples
rare <- which(occur[,2]<5)
rare

#determine the common species that occur in more than 95 samples
common<- which(occur[,2]>95)
common

#determine the species needed for the cca analysis
reduced<-maramacros[,-c(rare,common)]
str(reduced)

##Species transformations and standardizations
#First, check if species abundances are normally distributed across sites for the reduced species data set.
#First, check if species abundances are normally distributed across sites:
mapply(hist,as.data.frame(reduced),main=colnames(reduced),xlab="abundance")

##As you can see, most of the macros taxa distributions are right skewed. Use the log transformation (logx+1) to transform the species distributions for both the full and reduced datasets:
#log transformation of the total data set
log.full<-log1p(maramacros)
log.full
#log transformation of the reduced data set
log.red<-log1p(reduced)
log.red

####Next, check the row and column sum variability using the coefficient of variation (cv) for both data sets:

#Full data set:
rsum<-rowSums(log.full)
csum<-colSums(log.full)
cv(csum)

#Reduced data set:
rsumRed<-rowSums(log.red)
csumRed<-colSums(log.red)
cv(rsumRed)

#If either the row or column sums have cv >50, standardize by the total:

#If either the row or column sums have cv >50, standardize by the total:
cSpec<-sweep(log.full,2,csum,"/")
cSpec

#Reduced data set:
cSpecRed<-sweep(log.red,2,csumRed,"/")
cSpecRed

##Now that the data are reduced, transformed and standardized, you need to determine if species abundances show a linear (RDA) or a unimodal (CCA) relationship with the underlying gradient.
##First, use Detrended Correspondence Analysis (DCA) to determine the length of the canonical axes. 
##You will use the decorana function in the vegan Library. While DCA is a separate analysis with its own assumptions and multifaceted output, you will focus on axis length. 
##An axis length > 3 is evidence of a unimodal relationship(CCA). An axis length of <3 is evidence of a linear relationship(RDA). 
#total species 
decorana(cSpec)

#Reduced species
decorana(cSpecRed)

#Explanatory Variables
#Constrained ordination affords you the ability to include explanatory variables in the ordination. 
#You want to avoid mullitcolinearity among explanatory variables and check if they are measured on the same scale. 
#First look at all of the pairwise correlations between these variables: 
Vars<-maraenv[,c(1,3,4,5)]
Vars
round(as.dist(cor(Vars)),2)

#Do the environmental variables look like they are measured on different scale?
#Check the cv to see if you need to z-standardize them:
cv(colSums(Vars))

#scale the data
env<-as.data.frame(scale(Vars))
env

##Constrained Ordination of the reduced species data set
##Run CCA:
sp.CCA<-cca(cSpecRed~.,data=env)
summary(sp.CCA)

##Function for plotting species abundances vs. CCA Axis 1:
f2 <- function(x) {
  plot(x~sp.CCA$CC$wa[,1],xlab="CCA AXIS 1", ylab= "Abundance ")
}
#You will run the constrained ordination using the cca in the vegan library.
#Unconstrained Ordination (CA)
#Before running the constrained model, run an unconstrained ordination (i.e. a regular Correspondence Analysis. 
#CA will give you a measure of the amount of variation in the site by species matrix that you will try to explain with the explanatory variables (i.e. constraints). 
#Full Data
caRed<-cca(cSpecRed)
summary(caRed)

# Reduce margins (bottom, left, top, right)
par(mar = c(3, 3, 2, 2))  # Smaller margins
plot(caRed)

#Running the CCA using reduced data set
#The first thing you should focus on in the summary is the proportion of “inertia” (i.e. variance) explained by the Constrained Ordination. 
#Notice that the total amount of inertia is the same as the Unconstrained Ordination you just ran.
#Now look at the eigenvalue and proportion and cumulative amount of variation.
sp.CCA<-cca(cSpecRed~.,data=env)
summary(sp.CCA)

#Monte Carlo testing of the significance of the constrained axis.
#The permutation allows you to test if the constrained axes explain more variation than would be expected randomly. 
#You will use the anova.cca function in vegan to conduct the permutation. 
#It is “anova-like” but not an anova.  Global Test (i.e. all variables together):
anova(sp.CCA)

##Axes Tests (i.e. each axis individually):
anova(sp.CCA,by='axis')

##Variable Tests (i.e. each variable individually):
anova(sp.CCA,by='terms')

#Observed (F matrix) and Predicted (Z Matrix) Site Scores
#Now look back at you cca summary again:
summary(sp.CCA)

#Variable Selection: forward selection using ordiR2step in vegan
#Here we want to reduce the number of variables while maintaining a model that describes as much variance as possible. Here selection is based on increasing and R2 and p-value of the permutation test. See settings in ordiR2step for more options. 
step.forward <- ordiR2step(cca(cSpecRed ~ 1,data=env), 
                           scope=formula(sp.CCA), R2scope = F, direction="forward", pstep=1000)
#Partial CCA
partial.elev <- cca( cSpecRed~ elevation + Condition(velocity + temperature + DO) , data=env)
summary(partial.elev)

#The matrix labeled “Site scores (weighted averages of species scores)” is the F matrix and the matrix labeled “Site constraints (linear combinations of constraining variables)”is the Z matrix.

#Look at these two sets of site scores projected in ordination space:

plot(sp.CCA$CC$wa[,1],sp.CCA$CC$wa[,2],xlab="CCA AXIS 1", ylab= "CCA AXIS 2")
plot(sp.CCA$CC$u[,1],sp.CCA$CC$u[,2],xlab="CCA AXIS 1", ylab= "CCA AXIS 2")

#These correlations can lend insight as to how well the predicted site locations match the observed ones. 
#However, they are not to be trusted as the only line of evidence.
spenvcor(sp.CCA)

#Intra-set correlations and biplot scores for the constraining variables.
#Correlations between the Z matrix (predicted site scores) and the environmental variables provide information on which variables have the largest influence on the constrained ordination. These also denote the placement of the environmental variables as vectors on the CCA tri-plot.
sp.CCA$CCA$biplot
###############################################################################
######The Tri-Plot (using the site scores from the F matrix):
#adjust the magrin for the CCA1 axis to -5,5
par(mar=c(4,4,2,2))

plot(sp.CCA,choices=c(1,2),display=c('wa','sp','bp'),scaling=2)
title("CCA Tri-Plot")



#and using the site scores from the Z matrix:
par(mar=c(5,5,2,2))
plot(sp.CCA,choices=c(1,2),display=c('lc','sp','bp'),scaling=2)
#add title to the plot
title("CCA Tri-Plot")

#Steps to Create a Partitioned Variance Diagram:
#Draw the Diagram
library(VennDiagram)

# Variance Partitioning between constrained and unconstrained variance 
venn.plot <- draw.pairwise.venn(
  area1 = 52.79,  # Elevation and Water quality 
  area2 = 47.21,  # Unconstrained variance
  cross.area = 0, # Shared variance (set to 0 if no overlap)
  category = c("Elevation and Water Quality", "Unconstrained"),
  fill = c("grey", "black"),         # Keep distinct colors
  lty = "blank",                     # No border lines
  cat.cex = 1,                     # Increase the font size of titles
  cat.pos = c(-4, 4),              # Adjust title positions
  cat.dist = c(0.07, 0.07)           # Increase distance from circles for visibility
)
grid.draw(venn.plot)


#Forward selection between the constrained variables
venn.plot <- draw.pairwise.venn(
  area1 = 15.7,  # Elevation and Temperature 
  area2 = 12.68,  # DO and Velocity
  cross.area = 0, # Shared variance (set to 0 if no overlap)
  category = c("Elevation and Temperature", "DO and Velocity"),
  fill = c("grey", "lightgrey"),         # Keep distinct colors
  lty = "blank",                     # No border lines
  cat.cex = 1,                     # Increase the font size of titles
  cat.pos = c(-4, 4),              # Adjust title positions
  cat.dist = c(0.07, 0.07)           # Increase distance from circles for visibility
)
grid.draw(venn.plot)



