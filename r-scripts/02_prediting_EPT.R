# Multiple regression in a glm context:Quadratic effect of elevation on EPT distribution
# Q: How does the EPTA distribution along the Mara River gradient change between years and elevation?

# clear everything in memory (of R and restore libraries)
remove(list=ls())

#restore libraries to the environment if shared with other collaborators 

renv::restore()

# load the tidyverse libraries since it was in the R environment

library(tidyverse)

#database used  (remove hashtag)
#browseURL("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/edit?usp=sharing")

# read the macros data abundance data, call the dataset macrodat
# filter to use only years (2021,2022,2023), since 2024 sampling explored new sites which are not comparable to the previous years. 
# group by year and Location_Code
# calculate the sum of the number of EPT found per year and Location_ID
# do all the above in one pipeline

macrosdat <- readr::read_csv ("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1254679428&single=true&output=csv")|>
  dplyr::filter(year %in% c(2021, 2022, 2023)) |>
  dplyr::group_by(year, Location_ID,Order) |>
  dplyr::summarise(CountSum = sum(Count, na.rm = TRUE))
print(macrosdat)

# read the Fact Elevation data, filter and select the elevation in 3 years (2021,2022,2023) only
elevdat<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=11740542&single=true&output=csv") |># read the data
  dplyr::filter(year %in% c(2021, 2022, 2023)) # filter to retain only 3 years (2021,2022,2023)
print(elevdat)

# join  the elevation with  the macrosdata data by year and Location_Code, and filter to retain only 3 years (2021,2022,2023) and remove rows with missing values
# make year a factor
# call the new dataset macroselev

macroselev<-macrosdat|>
  dplyr::left_join(elevdat, by=c("year","Location_ID"))|># join the two datasets
  dplyr::filter(!is.na(elevation))|># remove rows with missing values
  dplyr::mutate(year=factor(year)) # make year a factor
names(macroselev)
# explore how Ephemeroptera abundance changes along the Mara River gradient in a scatter plot
macroselev|>  dplyr::filter(Order %in% c("Ephemeroptera","Plecoptera","Trichoptera")) |>
  ggplot2::ggplot(mapping=aes(x=elevation, y=CountSum)) + # x is the location, y is the count, group is the Year
  geom_point(aes(shape=year),size=3) +
  geom_smooth(method = "glm", , se = F, formula = y ~ x+I(x^2),
              method.args = list(family = "poisson"), col="black") +
  ylab("Count") + xlab("Elevation (m.a.s.l)") +
  facet_wrap(~Order,ncol=1,scales="free") # what is in rows ~what is in columns

#make the plot background white with no grid lines:
macroselev |>
  dplyr::filter(Order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera")) |>
  ggplot2::ggplot(mapping = aes(x = elevation, y = CountSum)) + 
  geom_point(aes(shape = year), size = 3) +
  geom_smooth(
    method = "glm",
    formula = y ~ x + I(x^2),
    se = FALSE,
    method.args = list(family = "poisson"),
    color = "black"
  ) +
  ylab("Count") + 
  xlab("Elevation (m.a.s.l)") +
  facet_wrap(~ Order, ncol = 1, scales = "free") +
  theme_minimal() +  # Start with a minimal theme
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # White background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    strip.background = element_rect(fill = "white"),  # White background for facet labels
    strip.text = element_text(size = 12, face = "bold")  # Customize facet label text
  )

# calculate what is the best model for Ephemeroptera
model_Ephemeroptera_lin<-glm(CountSum~elevation, 
                         data=macroselev |>dplyr::filter(Order %in% c("Ephemeroptera")),
        family = poisson(link=log)) # glm model
anova(model_Ephemeroptera_lin,test="Chisq")
model_Ephemeroptera_qua<-glm(CountSum~elevation+I(elevation^2), 
                             data=macroselev |>dplyr::filter(Order %in% c("Ephemeroptera")),
                             family = poisson(link=log)) # Quadratic model
anova(model_Ephemeroptera_qua,model_Ephemeroptera_lin,test="Chisq")

#calculate what is the best model for Plecoptera
model_Plecoptera_lin<-glm(CountSum~elevation, 
                             data=macroselev |>dplyr::filter(Order %in% c("Plecoptera")),
                             family = poisson(link=log)) # glm model
anova(model_Plecoptera_lin,test="Chisq")
model_Plecoptera_qua<-glm(CountSum~elevation+I(elevation^2), 
                             data=macroselev |>dplyr::filter(Order %in% c("Plecoptera")),
                             family = poisson(link=log)) # Quadratic model
anova(model_Plecoptera_qua,model_Plecoptera_lin,test="Chisq")

#calculate what is the best model for Trichoptera
model_Trichoptera_lin<-glm(CountSum~elevation, 
                          data=macroselev |>dplyr::filter(Order %in% c("Trichoptera")),
                          family = poisson(link=log)) # glm model
anova(model_Trichoptera_lin,test="Chisq")

model_Trichoptera_qua<-glm(CountSum~elevation+I(elevation^2), 
                          data=macroselev |>dplyr::filter(Order %in% c("Trichoptera")),
                          family = poisson(link=log)) # Quadratic model
anova(model_Trichoptera_qua,model_Trichoptera_lin,test="Chisq")

# plot the best model for EPT
ggsave("plots/Fig_best EPT model_Mara River.png", width = 6, height = 4, dpi=300, units = "in")


#%EPT (Ephemeroptera, Plecoptera, Trichoptera) against elevation

#EPT Taxa agaisnt elevation GLM)
# Load necessary libraries
library(readr)
library(dplyr)
library(ggplot2)

# Load and process macroinvertebrate data
macrosdat <- read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1254679428&single=true&output=csv") |>
  filter(year %in% c(2021, 2022, 2023)) |>
  group_by(year, Location_ID, Order) |>
  summarise(CountSum = sum(Count, na.rm = TRUE), .groups = "drop")
print(macrosdat)

# Load and process elevation data
elevdat <- read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=11740542&single=true&output=csv") |>
  filter(year %in% c(2021, 2022, 2023))
print(elevdat)

# Join datasets and calculate percentage EPT
macroselev <- macrosdat |>
  left_join(elevdat, by = c("year", "Location_ID")) |>
  group_by(year, Location_ID) |>
  mutate(
    TotalCount = sum(CountSum, na.rm = TRUE),
    PercentEPT = if_else(Order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera"),
                         (CountSum / TotalCount) * 100, 0)
  ) |>
  filter(!is.na(elevation)) |>
  mutate(year = factor(year)) |>
  ungroup()
print(macroselev)

# Aggregate to calculate total %EPT per site per year
macroselev_summary <- macroselev |>
  filter(Order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera")) |>
  group_by(year, Location_ID, elevation) |>
  summarise(PercentEPT = sum(PercentEPT, na.rm = TRUE), .groups = "drop")



# Plot %EPT against elevation with continuous black margin lines
ggplot(macroselev_summary, aes(x = elevation, y = PercentEPT)) +
  geom_point(aes(shape = year), size = 3) +
  geom_smooth(
    method = "glm",
    formula = y ~ x + I(x^2),
    method.args = list(family = "poisson"),
    se = FALSE,
    color = "black"
  ) +
  ylab("% EPT") +
  xlab("Elevation (m.a.s.l)") +
  ggtitle("Percentage of EPT Taxa vs Elevation") +
  theme_minimal()


# Save the final plot
ggsave("plots/Fig_Percent_EPT_vs_Elevation.png", width = 6, height = 4, dpi = 300, units = "in")
########################################################################################

#database used  (remove hashtag)
#browseURL("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/edit?usp=sharing")

# read the macros data abundance data, call the dataset macrodat
# filter to use only years (2008,2009).
# group by year and Location_Code
# calculate the sum of the number of EPT found per year and Location_ID
# do all the above in one pipeline

histmacrosdat <- readr::read_csv ("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1970181164&single=true&output=csv")|>
  dplyr::filter(year %in% c(2008,2009)) |>
  dplyr::group_by(year, Location_ID,Order) |>
  dplyr::summarise(CountSum = sum(Count, na.rm = TRUE))
print(histmacrosdat)

# read the Fact Elevation data, filter and select the elevation in 3 years (2021,2022,2023) only
elevdat<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=11740542&single=true&output=csv") |># read the data
  dplyr::filter(year %in% c(2021, 2022, 2023)) # filter to retain only 3 years (2021,2022,2023)
print(elevdat)

# join  the elevation with  the macrosdata data by year and Location_Code, and filter to retain only 3 years (2021,2022,2023) and remove rows with missing values
# make year a factor
# call the new dataset macroselev

macroselev<-macrosdat|>
  dplyr::left_join(elevdat, by=c("year","Location_ID"))|># join the two datasets
  dplyr::filter(!is.na(elevation))|># remove rows with missing values
  dplyr::mutate(year=factor(year)) # make year a factor
names(macroselev)
# explore how Ephemeroptera abundance changes along the Mara River gradient in a scatter plot
macroselev|>  dplyr::filter(Order %in% c("Ephemeroptera","Plecoptera","Trichoptera")) |>
  ggplot2::ggplot(mapping=aes(x=elevation, y=CountSum)) + # x is the location, y is the count, group is the Year
  geom_point(aes(shape=year),size=3) +
  geom_smooth(method = "glm", , se = F, formula = y ~ x+I(x^2),
              method.args = list(family = "poisson"), col="black") +
  ylab("Count") + xlab("Elevation (m.a.s.l)") +
  facet_wrap(~Order,ncol=1,scales="free") # what is in rows ~what is in columns

#make the plot background white with no grid lines:
macroselev |>
  dplyr::filter(Order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera")) |>
  ggplot2::ggplot(mapping = aes(x = elevation, y = CountSum)) + 
  geom_point(aes(shape = year), size = 3) +
  geom_smooth(
    method = "glm",
    formula = y ~ x + I(x^2),
    se = FALSE,
    method.args = list(family = "poisson"),
    color = "black"
  ) +
  ylab("Count") + 
  xlab("Elevation (m.a.s.l)") +
  facet_wrap(~ Order, ncol = 1, scales = "free") +
  theme_minimal() +  # Start with a minimal theme
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # White background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    strip.background = element_rect(fill = "white"),  # White background for facet labels
    strip.text = element_text(size = 12, face = "bold")  # Customize facet label text
  )

# calculate what is the best model for Ephemeroptera
model_Ephemeroptera_lin<-glm(CountSum~elevation, 
                             data=macroselev |>dplyr::filter(Order %in% c("Ephemeroptera")),
                             family = poisson(link=log)) # glm model
anova(model_Ephemeroptera_lin,test="Chisq")
model_Ephemeroptera_qua<-glm(CountSum~elevation+I(elevation^2), 
                             data=macroselev |>dplyr::filter(Order %in% c("Ephemeroptera")),
                             family = poisson(link=log)) # Quadratic model
anova(model_Ephemeroptera_qua,model_Ephemeroptera_lin,test="Chisq")

#calculate what is the best model for Plecoptera
model_Plecoptera_lin<-glm(CountSum~elevation, 
                          data=macroselev |>dplyr::filter(Order %in% c("Plecoptera")),
                          family = poisson(link=log)) # glm model
anova(model_Plecoptera_lin,test="Chisq")
model_Plecoptera_qua<-glm(CountSum~elevation+I(elevation^2), 
                          data=macroselev |>dplyr::filter(Order %in% c("Plecoptera")),
                          family = poisson(link=log)) # Quadratic model
anova(model_Plecoptera_qua,model_Plecoptera_lin,test="Chisq")

#calculate what is the best model for Trichoptera
model_Trichoptera_lin<-glm(CountSum~elevation, 
                           data=macroselev |>dplyr::filter(Order %in% c("Trichoptera")),
                           family = poisson(link=log)) # glm model
anova(model_Trichoptera_lin,test="Chisq")

model_Trichoptera_qua<-glm(CountSum~elevation+I(elevation^2), 
                           data=macroselev |>dplyr::filter(Order %in% c("Trichoptera")),
                           family = poisson(link=log)) # Quadratic model
anova(model_Trichoptera_qua,model_Trichoptera_lin,test="Chisq")

# plot the best model for EPT
ggsave("plots/Fig_best EPT model_Mara River.png", width = 6, height = 4, dpi=300, units = "in")
