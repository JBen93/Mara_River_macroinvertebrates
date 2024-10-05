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
macros2<-macros %>% select(-c(Location_ID, month, year, River_reach)) #remove the columns that are not needed for the analysis
head(macros2)
# Convert data from wide to long format
macros_long<- tidyr::pivot_longer(macros2, -Sample_ID, names_to="Taxa", values_to="Count")
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
D_index <- tapply(macros_long$Count, macros_long$Sample_ID, simpson_diversity_index)
print(D_index)
# For a comprehensive view, including evenness, you can create a summary table
# Calculate alpha diversity metrics
# Create a summary table
alpha_diversity_summary <- data.frame(Site=unique(macros_long$Sample_ID))
# Calculate Simpson Diversity Index for each site
#Calculate Simpson Diversity Index for each site

alpha_diversity_summary$Simpson_Diversity_Index <- tapply(macros_long$Count, macros_long$Sample_ID, simpson_diversity_index)
# Calculate Shannon Diversity Index for each site
alpha_diversity_summary$Shannon_Diversity_Index <- tapply(macros_long$Count, macros_long$Sample_ID, function(x) diversity(x, index="shannon"))
# Calculate species evenness
alpha_diversity_summary$Evenness <- alpha_diversity_summary$Simpson_Diversity_Index / length(unique(macros_long$Taxa))
print(alpha_diversity_summary)

#Using box plot to plot the Shannon diversity index across the river reaches

# make a data frame for each historical site (M2, M3, M5, M9) with the Shannon diversity index
  # For reproducibility
  set.seed(123)
  # Load necessary libraries

  # Shannon diversity data for each site
  M2_shannon <- c(1.189, 0.44, 1.202, 0.496, 1.712, 1.586, 1.680, 1.161)
  M3_shannon <- c(0.767, 0.405, 1.066, 1.007, 1.555, 0.794, 0.492, 0.857)
  M5_shannon <- c(0.613, 0.708, 0.768, 0.858, 0.775, 0.776, 1.10, 0.614)
  M9_shannon <- c(1.169, 0.574, 1.272, 0.757, 1.079, 0.529, 1.123, 0.485)
  
  # Combine the data into a data frame
  historical_shannon <- data.frame(
    Index = c(M2_shannon, M3_shannon, M5_shannon, M9_shannon),
    Reach = factor(rep(c("M2", "M3", "M5", "M9"), each = length(M2_shannon)))
  )
head(historical_shannon)


#Create the boxplot with mean points and error bars
ggplot(historical_shannon, aes(x = Reach, y = Index)) +
  geom_boxplot(outlier.shape = NA, fill = "grey") +  # Boxplot with grey fill
  
  # Add labels and adjust the axis labels
  labs(
    title = "2008-2009 Macroinvertebrate taxa",
    x = "Site",
    y = "Shannon-Wiener Index"
  ) +
  
  # Set the base theme and customize the font sizes
  theme_minimal() +  # Use theme_minimal for a clean look
  
# Customize the theme: remove panel background, adjust text size, and center the title
  theme(
    panel.background = element_blank(),   # Remove background
    axis.text = element_text(size = 14),  # Axis text size (e.g., tick labels)
    axis.title = element_text(size = 14),  # Axis title size
    plot.title = element_text(size = 14, hjust = 0.5))  # Title size and center the title

# Perform one-way ANOVA to test for differences in Shannon-Wiener Index across the sites
histoanova_result <- aov(Index ~ Reach, data = historical_shannon)

# Summary of the ANOVA test
summary(histoanova_result)

# Test for normality using the Shapiro-Wilk test to see if ANOVA assumptions are met
histoshapiro_test <- shapiro.test(residuals(histoanova_result))

# Print the result
print(histoshapiro_test)

# Perform Tukey's HSD test if ANOVA is significant
histotukey_result <- TukeyHSD(histoanova_result)

# Print the result
print(histotukey_result)
########################################################################################################################################
#calculate the species accumulation curve for the historical data
#Generate the species accumulation curve
histospecies_data<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1460274792&single=true&output=csv")

# Remove species with all-zero occurrences (no occurrences across any samples)
histospecies_data_clean <- histospecies_data[, colSums(histospecies_data) > 0]


# Plot the species accumulation curve without points or crosses
# Load necessary library
library(vegan)

# Run the species accumulation curve analysis
species_accum <- specaccum(histospecies_data_clean)

# Plot the species accumulation curve with only a black line
plot(species_accum, 
     main = "2008-2009 Species Accumulation Curve",  # Title of the plot
     xlab = "Number of Samples",           # Label for x-axis
     ylab = "Cumulative Number of Species",# Label for y-axis
     col = "black",                        # Set the line color to black
     lwd = 2)                              # Line width for the curve

# Overlay a smoothed line using the spline function
lines(spline(1:length(species_accum$richness), species_accum$richness), col = "black", lwd = 2)

#####################################################################################################
#Current data
#####################################################################################################
#calculate the diversity indices for the current macros dataset 2021,2022,2023
#database source 
#browseURL("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/edit?usp=sharing")
currentmacros<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=798918621&single=true&output=csv") |> 
  dplyr::filter(year %in% c(2021, 2022, 2023))
currentmacros2<-currentmacros %>% select(-c(observation_ID, month, year, River_reach)) #remove the columns that are not needed for the analysis
head(currentmacros2)
# Convert data from wide to long format
currentmacros_long<- tidyr::pivot_longer(currentmacros2, -sample_ID, names_to="Taxa", values_to="Count")
head(currentmacros_long)

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
D_index <- tapply(currentmacros_long$Count, currentmacros_long$sample_ID, simpson_diversity_index)
print(D_index)
# For a comprehensive view, including evenness, you can create a summary table
# Calculate alpha diversity metrics
# Create a summary table
alpha_diversity_summary <- data.frame(Site=unique(currentmacros_long$sample_ID))
# Calculate Simpson Diversity Index for each site
#Calculate Simpson Diversity Index for each site

alpha_diversity_summary$Simpson_Diversity_Index <- tapply(currentmacros_long$Count, currentmacros_long$sample_ID, simpson_diversity_index)
# Calculate Shannon Diversity Index for each site
alpha_diversity_summary$Shannon_Diversity_Index <- tapply(currentmacros_long$Count, currentmacros_long$sample_ID, function(x) diversity(x, index="shannon"))
# Calculate species evenness
alpha_diversity_summary$Evenness <- alpha_diversity_summary$Simpson_Diversity_Index / length(unique(currentmacros_long$Taxa))
print(alpha_diversity_summary)

#Using box plot to plot the Shannon diversity index across the river reaches

# make a data frame for each historical site (M2, M3, M5, M9) with the Shannon diversity index
# For reproducibility
set.seed(123)
# Shannon diversity data for each site
currentM2_shannon <- c(1.147,1.336,1.293,0.966,1.024)
currentM3_shannon <- c(1.641,1.728,1.737,1.507,0.515)
currentM5_shannon <- c(1.175,2.149,1.396,1.369,1.387)
currentM9_shannon <- c(1.581,1.444,0.96,1.447,0.816)

# Combine the data into a data frame
current_shannon <- data.frame(
  Index = c(currentM2_shannon, currentM3_shannon, currentM5_shannon, currentM9_shannon),
  Reach = factor(rep(c("M2", "M3", "M5", "M9"), each = length(currentM2_shannon)))
)
head(current_shannon)


#Create the boxplot with mean points and error bars
ggplot(current_shannon, aes(x =Reach, y =Index)) +
  geom_boxplot(outlier.shape = NA, fill = "grey") +  # Boxplot with grey fill
  
  # Add labels and adjust the axis labels
  labs(
    title = "2021-2023 Macroinvertebrate taxa",
    x = "Site",
    y = "Shannon-Wiener Index"
  ) +
  
  # Set the base theme and customize the font sizes
  theme_minimal() +  # Use theme_minimal for a clean look
  
  # Customize the theme: remove panel background, adjust text size, and center the title
  theme(
    panel.background = element_blank(),   # Remove background
    axis.text = element_text(size = 14),  # Axis text size (e.g., tick labels)
    axis.title = element_text(size = 14),  # Axis title size
    plot.title = element_text(size = 14, hjust = 0.5))  # Title size and center the title

# Perform one-way ANOVA to test for differences in Shannon-Wiener Index across the sites
anova_result <- aov(Index ~ Reach, data = current_shannon)

# Summary of the ANOVA test
summary(anova_result)

# Test for normality using the Shapiro-Wilk test to see if ANOVA assumptions are met
shapiro_test <- shapiro.test(residuals(anova_result))

# Print the result
print(shapiro_test)

# Perform Tukey's HSD test if ANOVA is significant
tukey_result <- TukeyHSD(anova_result)

# Print the result
print(tukey_result)

#######################################################################################################################################
# Load necessary library
library(vegan)
#calculate the species accumulation curve for the current data
#Generate the species accumulation curve
currentspecies_data<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1215461950&single=true&output=csv")

# Remove species with all-zero occurrences (no occurrences across any samples)
currentspecies_data_clean <- currentspecies_data[, colSums(currentspecies_data) > 0]


# Run the species accumulation curve analysis
species_accum <- specaccum(currentspecies_data_clean)

# Plot the species accumulation curve with only a black line
plot(species_accum, 
     main = "2021-2023 Species Accumulation Curve",  # Title of the plot
     xlab = "Number of Samples",           # Label for x-axis
     ylab = "Cumulative Number of Species",# Label for y-axis
     col = "black",                        # Set the line color to black
     lwd = 2)                              # Line width for the curve

#################################################################################################
#determining the species diversity difference between the historical and current data
# Assuming both historical_shannon and current_shannon have the same structure
historical_shannon$Period <- "Historical"
current_shannon$Period <- "Current"

# Combine the datasets
combined_shannon <- rbind(historical_shannon, current_shannon)
combined_anova_result <- aov(Index ~ Reach * Period, data = combined_shannon)

# Summary of the ANOVA test
summary(combined_anova_result)

# Post-hoc Tukey test for Reach/site differences
tukey_reach <- TukeyHSD(combined_anova_result, "Reach")
print(tukey_reach)
