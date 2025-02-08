# Load necessary libraries
# Clear the environment
# Clear the environment
remove(list = ls())

# Set working directory
setwd("/Users/Joshua/Library/CloudStorage/OneDrive-UniversityofFlorida/R input")
getwd()

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Load the CSV file
file_path <- "~/OneDrive - University of Florida/R input/Landcover.csv"
df <- read.csv(file_path)

# Remove NA columns
df <- df %>% select(where(~ !all(is.na(.)))) 

# Ensure Year column exists
if (!"Year" %in% colnames(df)) {
  stop("Error: 'Year' column is missing in the dataset.")
}

# Convert data from wide to long format
df_long <- df %>%
  pivot_longer(cols = Agriculture:Grassland, 
               names_to = "Land_Cover_Type", 
               values_to = "Area") %>%
  drop_na(Area)  # Remove rows where Area is NA

# Check if historical and current datasets exist
historical <- df_long %>% filter(Year == 2009)
current <- df_long %>% filter(Year >= 2021)

if (nrow(historical) == 0) stop("Error: No historical data (Year 2009) found.")
if (nrow(current) == 0) stop("Error: No current data (Year >= 2021) found.")

# Merge historical and current data by Location and Land Cover Type
df_merged <- merge(historical, current, 
                   by = c("Location_ID", "Land_Cover_Type"),
                   suffixes = c("_2009", "_Current"))

# Ensure df_merged contains data
if (nrow(df_merged) == 0) {
  stop("Error: No matching data found between historical and current records.")
}

# Compute differences
df_merged$Diff <- df_merged$Area_Current - df_merged$Area_2009

# Check normality of differences
shapiro_test <- shapiro.test(df_merged$Diff)

# Choose appropriate test based on normality
if (shapiro_test$p.value > 0.05) {
  test_result <- t.test(df_merged$Area_2009, df_merged$Area_Current, paired = TRUE)
  test_type <- "Paired t-test"
} else {
  test_result <- wilcox.test(df_merged$Area_2009, df_merged$Area_Current, paired = TRUE)
  test_type <- "Wilcoxon signed-rank test"
}

# Print results
print(paste("Using:", test_type))
print(test_result)

# Visualize land cover changes
ggplot(df_merged, aes(x = Land_Cover_Type, y = Diff, fill = Land_Cover_Type)) +
  geom_boxplot() +
  labs(title = "Land Cover Change Between 2009 and 2021–2023",
       y = "Change in Area",
       x = "Land Cover Type") +
  theme_minimal()
