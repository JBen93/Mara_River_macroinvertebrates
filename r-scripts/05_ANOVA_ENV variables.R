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

# Loop through numeric columns and generate boxplots
for (param in numeric_cols) {
  # Generate the plot
  plot <- ggplot(maraenv, aes(x = Site, y = !!sym(param), fill = Site)) +
    geom_boxplot(alpha = 0.7) +
    facet_wrap(~year) +
    labs(
      title = paste("Boxplot of", param, "by Site and Year"),
      x = "Site",
      y = param
    ) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3")
  
  # Save the boxplot with a unique filename
  ggsave(
    filename = paste0("plots/Boxplot_of_", param, ".png"),
    plot = plot,
    width = 8, 
    height = 6, 
    dpi = 300
  )
}
