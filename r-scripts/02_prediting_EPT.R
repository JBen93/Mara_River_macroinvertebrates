# Multiple regression in a glm context:Quadratic effect of elevation on EPT distribution
# Q: How does the EPTA distribution along the Mara River gradient change between years and elevation?

# clear everything in memory (of R and restore libraries)
remove(list=ls())

#restore libraries to the environment if shared with other collaborators 

renv::restore()

# load the tidyverse libraries since it was in the R environment

library(tidyverse)

# database used  (remove hashtag)
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
  ylab("count") + xlab("Elevation (m)") +
  facet_wrap(~Order,ncol=1,scales="free") # what is in rows ~what is in columns

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
ggsave("plots/Fig_best EPT model_Mara River.pdf", width = 6, height = 4, dpi=300, units = "in")






