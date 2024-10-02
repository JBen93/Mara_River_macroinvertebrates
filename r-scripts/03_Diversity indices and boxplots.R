# clear everything in memory (of R)
remove(list=ls())

renv::restore()

library(vegan) 
library(psych) 
library(tidyverse)
library(ggpubr)

#calculate the diversity indices for the historical macros dataset 2008-2009
#database source 
#browseURL("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/edit?usp=sharing")
#load data filter 2008,2009 and also group by Location_ID, month, year,River_reach and Family
macros<-readr::read_csv("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/pub?gid=1151562191&single=true&output=csv") |> 
  dplyr::filter(year %in% c(2008, 2009))
macros2<-macros %>% select(-c(Sample_ID, month, year, River_reach)) #remove the columns that are not needed for the analysis
head(macros2)
# Convert data from wide to long format
macros_long<- tidyr::pivot_longer(macros2, -Location_ID, names_to="Taxa", values_to="Count")
head(macros_long)

#simpson_diversity_index function
simpson_diversity_index <- function(species_counts) {
  total_count <- sum(species_counts)
  pi_squared_sum <- sum((species_counts / total_count) ^ 2)
D <- 1 - pi_squared_sum
  return(D)}

#shannon_weiner_index function
shannon_weiner_index <- function(species_counts) {
  total_count <- sum(species_counts)  # Total number of individuals
  proportions <- species_counts / total_count  # Proportion of each species
  
  # Calculate the Shannon index
  shannon_index <- -sum(proportions * log(proportions), na.rm = TRUE)
  
  return(shannon_index)
}
# Calculate Simpson Diversity Index for each site
D_index <- tapply(data_long$Count, data_long$Site, simpson_diversity_index)
print(D_index)
# For a comprehensive view, including evenness, you can create a summary table
# Calculate alpha diversity metrics
# Create a summary table
alpha_diversity_summary <- data.frame(Site=unique(data_long$Site))
# Calculate Simpson Diversity Index for each site
#Calculate Simpson Diversity Index for each site

alpha_diversity_summary$Simpson_Diversity_Index <- tapply(data_long$Count, data_long$Site, simpson_diversity_index)
# Calculate Shannon Diversity Index for each site
alpha_diversity_summary$Shannon_Diversity_Index <- tapply(data_long$Count, data_long$Site, function(x) diversity(x, index="shannon"))
# Calculate species evenness
alpha_diversity_summary$Evenness <- alpha_diversity_summary$Simpson_Diversity_Index / length(unique(data_long$Taxa))
print(alpha_diversity_summary)
