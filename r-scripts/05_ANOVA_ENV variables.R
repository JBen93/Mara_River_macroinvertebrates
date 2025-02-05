#Environmental variables/ water quality data ~ ANOVA and Boxplots for different sites and years

# clear everything in memory (of R)
remove(list=ls())

renv::restore()

library(vegan) 
library(psych) 
library(tidyverse)
library(ggpubr)
library(car)

#database source
#browseURL("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/edit?usp=sharing")
#Historical environmental variables 

# Read and process the dataset
histmaraenv <- readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=238930453&single=true&output=csv") |>
  dplyr::filter(year %in% c(2008, 2009)) |>                  # Filter for years 2008 and 2009
  dplyr::mutate(year = factor(year)) |>                     # Convert year to a factor
  dplyr::filter(Site %in% c("M2", "M3", "M5", "M9")) |>     # Filter for specific sites
  dplyr::distinct(Sample_ID, .keep_all = TRUE)              # Remove duplicate rows by Sample_ID

# Convert Sample_ID to row names and back to a tibble
histmaraenv <- histmaraenv %>%
  tibble::column_to_rownames(var = "Sample_ID") %>%          # Make Sample_ID row names
  tibble::rownames_to_column(var = "Sample_ID")             # Bring Sample_ID back as a column

# Select numeric columns for ANOVA
numeric_cols <- histmaraenv %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

# View numeric column names
print(numeric_cols)

# Initialize a results list to store ANOVA results for each parameter
histanova_results <- list()

# Loop through each numeric column and perform ANOVA by 'Site' and 'year'
for (param in numeric_cols) {
  cat("Performing ANOVA for parameter:", param, "\n")
  
  # One-way ANOVA for each parameter grouped by 'Site'
  model_site <- aov(as.formula(paste(param, "~ Site")), data = histmaraenv)
  summary_site <- summary(model_site)
  
  # One-way ANOVA for each parameter grouped by 'year'
  model_year <- aov(as.formula(paste(param, "~ year")), data = histmaraenv)
  summary_year <- summary(model_year)
  
  # Store results in the list
  histanova_results[[param]] <- list(
    "Site" = summary_site,
    "Year" = summary_year
  )
}


# Print the ANOVA results for each parameter
for (param in names(histanova_results)) {
  cat("\n\nANOVA Results for:", param, "\n")
  cat("Grouped by Site:\n")
  print(histanova_results[[param]]$Site)
  cat("\nGrouped by Year:\n")
  print(histanova_results[[param]]$Year)
}
# Initialize a list to store Post Hoc results
histpost_hoc_results <- list()

# Loop through numeric columns for ANOVA and Post Hoc
for (param in numeric_cols) {
  cat("Performing ANOVA and Post Hoc for parameter:", param, "\n")
  
  # One-way ANOVA grouped by 'Site'
  histmodel_site <- aov(as.formula(paste(param, "~ Site")), data = histmaraenv)
  histsummary_site <- summary(histmodel_site)
  
  # One-way ANOVA grouped by 'year'
  histmodel_year <- aov(as.formula(paste(param, "~ year")), data = histmaraenv)
  histsummary_year <- summary(histmodel_year)
  
  # Perform Tukey's HSD post hoc test if ANOVA is significant
  if (histsummary_site[[1]][["Pr(>F)"]][1] < 0.05) {
    tukey_site <- TukeyHSD(histmodel_site)
  } else {
    tukey_site <- "No significant difference by Site"
  }
  
  if (histsummary_year[[1]][["Pr(>F)"]][1] < 0.05) {
    tukey_year <- TukeyHSD(histmodel_year)
  } else {
    tukey_year <- "No significant difference by Year"
  }
  
  # Store results in the Post Hoc results list
  histpost_hoc_results[[param]] <- list(
    "ANOVA_Site" = histsummary_site,
    "PostHoc_Site" = tukey_site,
    "ANOVA_Year" = histsummary_year,
    "PostHoc_Year" = tukey_year
  )
}

# Print ANOVA and Post Hoc results
for (param in names(histpost_hoc_results)) {
  cat("\n\nResults for Parameter:", param, "\n")
  cat("ANOVA by Site:\n")
  print(histpost_hoc_results[[param]]$ANOVA_Site)
  cat("Post Hoc by Site:\n")
  print(histpost_hoc_results[[param]]$PostHoc_Site)
  cat("\nANOVA by Year:\n")
  print(histpost_hoc_results[[param]]$ANOVA_Year)
  cat("Post Hoc by Year:\n")
  print(histpost_hoc_results[[param]]$PostHoc_Year)
}
# Calculate the mean for each parameter grouped by Site
histmaraenv_mean <- histmaraenv %>%
  group_by(Site) %>%
  summarize(
    across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)), # Calculate mean with NA handling
    .groups = "drop"                                        # Drop grouping after summarization
  ) %>%
  pivot_longer(
    cols = all_of(numeric_cols),    # Use all numeric columns for pivoting
    names_to = "Parameter",         # Name of the parameter column
    values_to = "MeanValue"         # Column for mean values
  )

# View the calculated means
print(histmaraenv_mean)

# Convert the data to long format for visualization
histmaraenv_long <- histmaraenv %>%
  pivot_longer(
    cols = all_of(numeric_cols),  # Only numeric columns
    names_to = "Parameter",       # Column for parameter names
    values_to = "Value"           # Column for parameter values
  )

# View the long-format data
print(histmaraenv_long)



##########call in data for the current macros.

#database source
#browseURL("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/edit?usp=sharing")
#Historical environmental variables 

# Read and process the dataset

#Call the water quality data from the google sheets link
# Read and process the dataset
maraenv <- readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1736847188&single=true&output=csv") |>
  dplyr::filter(year %in% c(2021, 2022)) |>
  dplyr::mutate(year = factor(year)) |>
  dplyr::filter(Site %in% c("M2", "M3", "M5", "M9"))

# Convert Sample_ID to row names and back to a tibble
maraenv <- maraenv %>%
  tibble::column_to_rownames(var = "Sample_ID") %>% 
  tibble::rownames_to_column(var = "Sample_ID")  # Bring back Sample_ID as a column to keep it as a tibble

# Select numeric columns for ANOVA
numeric_cols <- maraenv %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

# View numeric column names
print(numeric_cols)

# Initialize a results list to store ANOVA results for each parameter
anova_results <- list()

# Loop through each numeric column and perform ANOVA by 'Site' and 'year'
for (param in numeric_cols) {
 cat("Performing ANOVA for parameter:", param, "\n")
  # One-way ANOVA for each parameter grouped by 'Site'
  model_site <- aov(as.formula(paste(param, "~ Site")), data = maraenv)
  summary_site <- summary(model_site)
 
   # One-way ANOVA for each parameter grouped by 'year'
  model_year <- aov(as.formula(paste(param, "~ year")), data = maraenv)
  summary_year <- summary(model_year)
 
   # Store results in the list
  anova_results[[param]] <- list(
    "Site" = summary_site,
    "Year" = summary_year
  )
}

# Print the ANOVA results for each parameter
for (param in names(anova_results)) {
  cat("\n\nANOVA Results for:", param, "\n")
  cat("Grouped by Site:\n")
  print(anova_results[[param]]$Site)
  cat("\nGrouped by Year:\n")
  print(anova_results[[param]]$Year)
}

# Loop through numeric columns for ANOVA and Post Hoc
post_hoc_results <- list()

for (param in numeric_cols) {
  cat("Performing ANOVA and Post Hoc for parameter:", param, "\n")
  
  # One-way ANOVA grouped by 'Site'
  model_site <- aov(as.formula(paste(param, "~ Site")), data = maraenv)
  summary_site <- summary(model_site)
  
  # One-way ANOVA grouped by 'year'
  model_year <- aov(as.formula(paste(param, "~ year")), data = maraenv)
  summary_year <- summary(model_year)
  
  # Perform Tukey's HSD post hoc test if ANOVA is significant
  if (summary_site[[1]][["Pr(>F)"]][1] < 0.05) {
    tukey_site <- TukeyHSD(model_site)
  } else {
    tukey_site <- "No significant difference by Site"
  }
  
  if (summary_year[[1]][["Pr(>F)"]][1] < 0.05) {
    tukey_year <- TukeyHSD(model_year)
  } else {
    tukey_year <- "No significant difference by Year"
  }
  
  # Store results
  post_hoc_results[[param]] <- list(
    "ANOVA_Site" = summary_site,
    "PostHoc_Site" = tukey_site,
    "ANOVA_Year" = summary_year,
    "PostHoc_Year" = tukey_year
  )
}

# Print ANOVA and post hoc results
for (param in names(post_hoc_results)) {
  cat("\n\nResults for Parameter:", param, "\n")
  cat("ANOVA by Site:\n")
  print(post_hoc_results[[param]]$ANOVA_Site)
  cat("Post Hoc by Site:\n")
  print(post_hoc_results[[param]]$PostHoc_Site)
  cat("\nANOVA by Year:\n")
  print(post_hoc_results[[param]]$ANOVA_Year)
  cat("Post Hoc by Year:\n")
  print(post_hoc_results[[param]]$PostHoc_Year)
}


# Calculate the mean for each parameter grouped by Site
maraenv_mean <- maraenv %>%
  group_by(Site) %>%
  summarize(across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(
    cols = all_of(numeric_cols), # Use all_of() for external vectors
    names_to = "Parameter",
    values_to = "MeanValue"
  )

# Convert the data to long format for visualization
# Convert the data to long format for visualization
maraenv_long <- maraenv %>%
  tidyr::pivot_longer(
    cols = all_of(numeric_cols),  # Only numeric columns
    names_to = "Parameter",       # Column for parameter names
    values_to = "Value"           # Column for parameter values
  )


# Create a boxplot with whiskers (min/max values) for each parameter, merging years
merged_boxplot <- ggplot(maraenv_long, aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~ Parameter, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Site",
    y = "Mean Parameter Value"
  ) +
  theme(
    strip.text = element_text(size = 10),  # Customize facet labels
    axis.text.x = element_text(angle = 45, hjust = 1)  # Adjust x-axis text
  )


# Print the plot
print(merged_boxplot)

########################### Boxplot for the two periods ############################

histmaraenv_long <- histmaraenv_long %>%
  mutate(Period = "2008-2009")

maraenv_long <- maraenv_long %>%
  mutate(Period = "2021-2022")

# Combine the two datasets into one
combined_long <- bind_rows(histmaraenv_long, maraenv_long)

# Create a combined boxplot using facet_wrap to compare both periods for each parameter
combined_boxplot <- ggplot(combined_long, aes(x = Site, y = Value, fill = Period)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~ Parameter, scales = "free_y") +
  scale_fill_manual(values = c("2008-2009" = "skyblue", "2021-2022" = "orange")) +  # Distinct colors for each period
  theme_minimal() +
  labs(
    x = "Site",
    y = "Mean Parameter Value"
  ) +
  theme(
    strip.text = element_text(size = 10),  # Customize facet labels
    axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust x-axis text
    legend.title = element_blank()  # Remove legend title for clarity
  )

# Print the combined plot
print(combined_boxplot)


##########################################################################################
  #Test for the differences in water quality parameters histmaraenv and maraenv
  #using ANOVA and Tukey's HSD post hoc test

# Function to perform ANOVA and Post Hoc Testing
test_anova_posthoc <- function(data, dataset_name) {
  numeric_cols <- data %>% select(where(is.numeric)) %>% colnames()
  
  results <- list()
  
  for (param in numeric_cols) {
    cat("\nPerforming ANOVA for parameter:", param, "in", dataset_name, "\n")
    
    # Two-way ANOVA with interaction
    model <- aov(as.formula(paste(param, "~ Site * year")), data = data)
    summary_model <- summary(model)
    
    # Perform Tukey's HSD post hoc test if ANOVA is significant
    tukey_results <- list()
    if (summary_model[[1]][["Pr(>F)"]][1] < 0.05) {
      tukey_results[["Site"]] <- TukeyHSD(model, "Site")
    }
    if (summary_model[[1]][["Pr(>F)"]][2] < 0.05) {
      tukey_results[["Year"]] <- TukeyHSD(model, "year")
    }
    if (summary_model[[1]][["Pr(>F)"]][3] < 0.05) {
      tukey_results[["Interaction"]] <- TukeyHSD(model, "Site:year")
    }
    
    results[[param]] <- list(
      "ANOVA" = summary_model,
      "PostHoc" = tukey_results
    )
  }
  
  return(results)
}

# Read and process historical dataset
histmaraenv <- readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=238930453&single=true&output=csv") %>%
  filter(year %in% c(2008, 2009), Site %in% c("M2", "M3", "M5", "M9")) %>%
  mutate(year = factor(year)) %>%
  distinct(Sample_ID, .keep_all = TRUE) %>%
  column_to_rownames(var = "Sample_ID") %>%
  rownames_to_column(var = "Sample_ID")

# Perform ANOVA and Post Hoc tests on historical data
hist_results <- test_anova_posthoc(histmaraenv, "Historical Data")

# Read and process recent dataset
maraenv <- readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1736847188&single=true&output=csv") %>%
  filter(year %in% c(2021, 2022), Site %in% c("M2", "M3", "M5", "M9")) %>%
  mutate(year = factor(year)) %>%
  column_to_rownames(var = "Sample_ID") %>%
  rownames_to_column(var = "Sample_ID")

# Perform ANOVA and Post Hoc tests on recent data
mara_results <- test_anova_posthoc(maraenv, "Recent Data")

# Function to print results
print_results <- function(results, dataset_name) {
  for (param in names(results)) {
    cat("\n\nResults for", dataset_name, "Parameter:", param, "\n")
    print(results[[param]]$ANOVA)
    cat("Post Hoc Tests:\n")
    print(results[[param]]$PostHoc)
  }
}

# Print results
print_results(hist_results, "Historical Data")
print_results(mara_results, "Recent Data")

##########################################################################################

#compare between the two datasets( histmaraenv and maraenv)
# Combine both datasets for comparison
combined_data <- bind_rows(histmaraenv %>% mutate(Period = "Historical"), 
                           maraenv %>% mutate(Period = "Recent"))
#remove variable pH from the combined data
combined_data <- combined_data %>% select(-pH)
# Perform ANOVA comparing historical vs recent data
test_combined_anova_posthoc <- function(data) {
  numeric_cols <- data %>% select(where(is.numeric)) %>% colnames()
  
  combined_results <- list()
  
  for (param in numeric_cols) {
    cat("\nPerforming ANOVA for parameter:", param, "across Periods\n")
    
    model <- aov(as.formula(paste(param, "~ Period * Site")), data = data)
    summary_model <- summary(model)
    
    tukey_results <- list()
    if (summary_model[[1]][["Pr(>F)"]][1] < 0.05) {
      tukey_results[["Period"]] <- TukeyHSD(model, "Period")
    }
    if (summary_model[[1]][["Pr(>F)"]][2] < 0.05) {
      tukey_results[["Site"]] <- TukeyHSD(model, "Site")
    }
    if (summary_model[[1]][["Pr(>F)"]][3] < 0.05) {
      tukey_results[["Interaction"]] <- TukeyHSD(model, "Period:Site")
    }
    
    combined_results[[param]] <- list(
      "ANOVA" = summary_model,
      "PostHoc" = tukey_results
    )
  }
  
  return(combined_results)
}

# Perform and print ANOVA for combined dataset
combined_results <- test_combined_anova_posthoc(combined_data)
print_results(combined_results, "Combined Historical vs Recent Data")

