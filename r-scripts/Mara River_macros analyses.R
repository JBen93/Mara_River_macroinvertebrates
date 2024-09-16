# Multiple regression in a glm context: calculation and plotting
# Q: How does the distribution of the macroinvertebrate along the Mara River gradient change between years?

# clear everything in memory (of R and restore libraries)
remove(list=ls())

#restore libraries to the environment if shared with other collaborators 

renv::restore()

# load the tidyverse libraries since it was in the R environment

library(tidyverse)


# read the macros data abundance data, call the dataset macrodat
# filter to use only years (2021,2022,2023) only since 2024 sampling explored new sites which are not comparable to the previous years. 
# group by year and Location_Code
# calculate the sum of the number of EPT found per year and Location_ID
# do all the above in one pipeline
# only use 3 or less replicates

macrosdat <- readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1254679428&single=true&output=csv")|>
  dplyr::filter(Year %in% c(2021, 2022, 2023)) |>
  dplyr::filter(Order=="Ephemeroptera") |>
  dplyr::group_by(Year, Location_Code) |>
  dplyr::summarise(CountSum = sum(Count, na.rm = TRUE))
print(macrosdat)

# read the Fact Elevation data, filter and select right variables
elevdat<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=880897856&single=true&output=csv") |>
  dplyr::filter(Year %in% c(2021,202,2023) & !is.na(Location_Code))|>
  dplyr::select(Year,Location_Code,Elevation)   # select  only distance_m and elevation 
elevdat

# join  the elevation with  the orchestia data by year and TransectPoint_ID, and filter to retain only transect points between 200 and 1000 (where it lives)
# make year a factor
orchdat3<-left_join(orchdat, elevdat, by=c("year", "TransectPoint_ID")) |>
  dplyr::filter(TransectPoint_ID>=200 & TransectPoint_ID<=1000) |>
  dplyr::filter(!(year==2022 & TransectPoint_ID==260))


# explore how Orchestia abundance changes along the transect in a bar plot
orchdat3 |>
  ggplot2::ggplot(mapping=aes(x=factor(TransectPoint_ID), y=CountSum, group = year)) +
  geom_bar(stat="identity") + # value that is already in the table (always need to specify statistic)
  facet_grid(year~.) # what is in rows vs what is in columns (years is in rows nothing is in columns)

# explore how elevation changes along the transect in a barplot


# plot Orchestia (y) versus elevation_m (x) in ggplot as a scatterplot, with each year as a different color
p1<- orchdat3  |> 
  ggplot2::ggplot(mapping=aes(x=elevation_m, y=CountSum, color=factor(year))) +
  geom_point(size=3) 
p1

# calculate the optimal preferred elevation by Orchestia for each year using weighted.mean function
orchdat3 |> 
  dplyr::group_by(year) |>
  dplyr::summarise(wa_elevation=weighted.mean(elevation_m, w=CountSum, na.rm=T)) #w = weighing factor)

## explore response to elevation and year as a linear model, call this m1
# first only elevation (not yet year)
orchdat3
m1<-lm(CountSum~elevation_m, data=orchdat3)
print(m1)
# so b0=-34.80 (intercept) and b1= 76.37 (coefficient or slope)
# so the model is y=-34.80 + 76.37x

#add the linear model to the plot
orchdat3$pred1<-predict(m1) #calculate the predicted value of m1 for every observation
p2 <- p1 + geom_line(data=orchdat3, aes(y=pred1), color="black", linewidth=1.2) # add the new predicted line to the previous plot p1)
print(p2)
#does not look like a good linear model

# fitmodel  m2 by adding a quadratic term for elevation_m to check for an ecological optimum
# y=b0 + b1x1 + b2*I(x^2) #expect that b2 is negative (gives us an optimum (berg parabool))
m2<-lm(CountSum~elevation_m + I(elevation_m^2), data=orchdat3) # I() is a function that tells R to treat the term as a single variable so we can transform it
print(m2)

#add the linear model to the plot
# calculate the predicted value of m2 for every observation, add to the dataset as a variable as pred2
# add the new predicted line to the previous plot p2, store as object p3 and show it
orchdat3$pred2<- predict(m2)
p2 + geom_line(data=orchdat3, aes(y=pred2), color="red", linewidth=1.2)
#this line is a bit better but it goes below zero

# test if the new model m2   significantly explains more variation than the first model m1
anova(m1,m2)
#statistically significant, so m2 is better than m1

# predict the orchestia abundance at 1.5 m elevation


# add year (as a factor) to the model and fit it as model m3 (does year explain something that is not already explained by elevation)
# so y = b0 + b1x1 + b2x1^2 +b3x2 
# test if it is significant, and better than the previous one


#add the linear model to the plot
# calculate the predicted value of m2 for every observation, add to the dataset as a variable as pred2
# add the new predicted line to the previous plot p2, store as object p3 and show it
m3<-lm(CountSum~elevation_m + I(elevation_m^2) + factor(year), data=orchdat3) #add year as a factor
print(m3)
anova(m3) #elevation has an quadratic effect but also the year has an effect
#is there interaction (response of elevation to year)?
orchdat3$pred3<-predict(m3)
p1 + geom_line(data=orchdat3, aes(y=pred3, color=factor(year)), linewidth=1.2)

# does this significantly explain more variation than the model without year?
anova(m3,m2)

# yes this is a better model

# include the interaction between elevation_m and year (as a factor) in the model 
m4<-lm(CountSum~elevation_m + I(elevation_m^2) + factor(year) + elevation_m*factor(year), data=orchdat3)
print(m4)

#is there a significant interaction between elevation and year?
orchdat3$pred4<-predict(m4)
p1 + geom_line(data=orchdat3, aes(y=pred4, color=factor(year)), linewidth=1.2)
anova(m4,m3)


#no, the interaction is not significant, so the model without interaction is better

#add the linear model to the plot
# calculate the predicted value of m2 for every observation, add to the dataset as a variable as pred2
# add the new predicted line to the previous plot p2, store as object p3 and show it


#### m3 is the preferred model in this case (outcome of model selection)

# explore the consequences of a log transformation of y values
x<-c(0,1,2,3,4,5)
y<-c(1,10,100,1000,10000,100000)
dat<-data.frame(x,y)
dat
dat %>% ggplot(aes(x=x,y=y)) +
  geom_point(shape=16,size=5) +
  geom_line(size=1.2)
dat %>% ggplot(aes(x=x,y=log10(y))) +
  geom_point(shape=16,size=3) +
  geom_line(size=1.2)
dat %>% ggplot(aes(x=x,y=log(y))) + # use log with base number e
  geom_point(shape=16,size=3) +
  geom_line(size=1.2)
# explore the consequences of the log base number (10 or e)
log(0)
log10(1)
log10(10)
log(0)
log(1)
exp(1)
log(exp(1))

### develop, test for significance and plot different models of increasing complexity 
#  using  multiple regresssion, assuming a poisson distribution (so use a generalized linear model)
# Explore  how the abundance of Orchestia depends on elevation_m and year,  their potential interaction,
# and a potential ecological optimum of Orchestia with respect to elevation_m
# show the effect of elevation but now in a generalized linear model instead of linear model, using a log link function and a poisson distribution


#add the linear model to the plot
# calculate the predicted value of m2 for every observation, add to the dataset as a variable as pred2
# add the new predicted line to the previous plot p2, store as object p3 and show it



# now test and show  the effect of both elevation , elevation squared and year


#add the linear model to the plot
# calculate the predicted value of m2 for every observation, add to the dataset as a variable as pred2
# add the new predicted line to the previous plot p2, store as object p3 and show it


# better than the previous?


# add the interaction to the model: elevation + elevation ^2 + year + elevation*year
# now test and show  the effect of both elevation + year


#add the  model to the plot
# calculate the predicted value of m2 for every observation, add to the dataset as a variable as pred2
# add the new predicted line to the previous plot p2, store as object p3 and show it





