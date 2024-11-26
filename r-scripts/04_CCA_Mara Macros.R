# clear everything in memory (of R)
remove(list=ls())

renv::restore()

library(vegan) 
library(psych) 
library(tidyverse)
library(ggpubr)
library(terra)
library(raster)
source("biostats.r") # this is a custom script that contains the functions used in this script

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
