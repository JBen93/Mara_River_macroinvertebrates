# Load necessary libraries
remove(list=ls())
setwd("/Users/Joshua/Library/CloudStorage/OneDrive-UniversityofFlorida/R input")
library(dplyr)
library(readxl)

# Load necessary libraries
library(dplyr)
library(readxl)

# Load the data from the Excel file
df <- read_excel("2008-2009 macros.xlsx")

# Clean column names to remove any leading/trailing spaces
colnames(df) <- trimws(colnames(df))

# Define the EPT orders
ept_orders <- c("Ephemeroptera", "Plecoptera", "Trichoptera")

# Filter the data to include only EPT families
ept_data <- df %>% filter(Order %in% ept_orders)

# Aggregate the total abundance across both years for each site
site_data <- df %>%
  group_by(Site) %>%
  summarise(Total_Abundance = sum(Abundance, na.rm = TRUE))

# Aggregate the EPT abundance across both years for each site
site_ept_data <- ept_data %>%
  group_by(Site) %>%
  summarise(EPT_Abundance = sum(Abundance, na.rm = TRUE))

# Merge the total and EPT abundance data
site_comparison <- merge(site_data, site_ept_data, by = "Site", all.x = TRUE)

# Calculate %EPT for each site
site_comparison <- site_comparison %>%
  mutate(Percent_EPT = (EPT_Abundance / Total_Abundance) * 100)

# Display the results
print(site_comparison)


#######################################################

# Load the data from the Excel file
# Load necessary libraries
library(dplyr)
library(readxl)

# Load the data from the Excel file
df <- read_excel("2021_2022 macros.xlsx", sheet = "FactBenthos")

# Clean column names
colnames(df) <- trimws(colnames(df))

# Define EPT orders
ept_orders <- c("Ephemeroptera", "Plecoptera", "Trichoptera")

# Filter the data for EPT families
ept_data <- df %>% filter(Order %in% ept_orders)

# Aggregate the data by Site, summing the abundance (Count) across both years
site_data <- df %>%
  filter(Site %in% c("M2", "M3", "M5", "M9")) %>%
  group_by(Site) %>%
  summarise(Total_Abundance = sum(Count, na.rm = TRUE))

# Aggregate the EPT data by Site, summing the abundance across both years
site_ept_data <- ept_data %>%
  filter(Site %in% c("M2", "M3", "M5", "M9")) %>%
  group_by(Site) %>%
  summarise(EPT_Abundance = sum(Count, na.rm = TRUE))

# Merge the total and EPT abundance data
site_comparison <- merge(site_data, site_ept_data, by = "Site", all.x = TRUE)

# Calculate %EPT for each site
site_comparison <- site_comparison %>%
  mutate(Percent_EPT = (EPT_Abundance / Total_Abundance) * 100)

# Display the results
print(site_comparison)
