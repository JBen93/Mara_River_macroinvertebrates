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

#Call the water quality data from the google sheets link
maraenv<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1736847188&single=true&output=csv")|> 
dplyr::filter(year %in% c(2021, 2022))|>
dplyr::mutate(year=factor(year))|>  # make year a factor
dplyr::filter(Site %in% c("M2", "M3", "M5", "M9"))|> # filter for specific sites
tibble::column_to_rownames(var="Sample_ID") # make Sample_ID the rownames

# Select numeric columns for ANOVA (excluding non-numeric columns like Sample_ID, Site, year)
numeric_cols <- maraenv %>% select(where(is.numeric)) %>% colnames()

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
maraenv_long <- maraenv %>%
  rownames_to_column(var = "Sample_ID") %>%
  pivot_longer(
    cols = all_of(numeric_cols), # Only numeric columns
    names_to = "Parameter",
    values_to = "Value"
  )

# Create a boxplot with whiskers (min/max values) for each parameter, merging years
merged_boxplot <- ggplot(maraenv_long, aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~ Parameter, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Box Plots of Environmental variables by Site (2021-2022)",
    x = "Site",
    y = "Parameter Value"
  ) +
  theme(
    strip.text = element_text(size = 10),  # Customize facet labels
    axis.text.x = element_text(angle = 45, hjust = 1)  # Adjust x-axis text
  )

# Print the plot
print(merged_boxplot)

####################################################################################################
#Historical environmental variables 

#Call the water quality data from the google sheets link
histmaraenv<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=238930453&single=true&output=csv")|> 
  dplyr::filter(year %in% c(2008, 2009))|>
  dplyr::mutate(year=factor(year))|>  # make year a factor
  dplyr::filter(Site %in% c("M2", "M3", "M5", "M9"))|> # filter for specific sites
  tibble::column_to_rownames(var="Sample_ID") # make Sample_ID the rownames

# Select numeric columns for ANOVA (excluding non-numeric columns like Sample_ID, Site, year)
numeric_cols <- histmaraenv %>% select(where(is.numeric)) %>% colnames()

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
  summarize(across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(
    cols = all_of(numeric_cols), # Use all_of() for external vectors
    names_to = "Parameter",
    values_to = "MeanValue"
  )

# Convert the data to long format for visualization
histmaraenv_long <- histmaraenv %>%
  rownames_to_column(var = "Sample_ID") %>%
  pivot_longer(
    cols = all_of(numeric_cols), # Only numeric columns
    names_to = "Parameter",
    values_to = "Value"
  )

# Create a boxplot with whiskers (min/max values) for each parameter, merging years
histmerged_boxplot <- ggplot(histmaraenv_long, aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~ Parameter, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Box Plots of Environmental variables by Site (2008-2009)",
    x = "Site",
    y = "Parameter Value"
  ) +
  theme(
    strip.text = element_text(size = 10),  # Customize facet labels
    axis.text.x = element_text(angle = 45, hjust = 1)  # Adjust x-axis text
  )

# Print the plot
print(histmerged_boxplot)





