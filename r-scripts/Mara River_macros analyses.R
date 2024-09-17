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

macrosdat <- readr::read_csv ("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1254679428&single=true&output=csv")|>
  dplyr::filter(year %in% c(2021, 2022, 2023)) |>
  dplyr::filter(Order=="Ephemeroptera") |>
  dplyr::group_by(year, Location_ID) |>
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




# explore how Ephemeroptera abundance changes along the Mara River gradient in a bar plot
macroselev|>
  ggplot2::ggplot(mapping=aes(x=factor(elevation), y=CountSum, group = year)) + # x is the location, y is the count, group is the Year
  geom_bar(stat="identity") + # value that is already in the table (always need to specify statistic)
  facet_grid(year~.) # what is in rows vs what is in columns (years is in rows nothing is in columns)

# explore how elevation changes along the transect in a barplot
# explore how elevation changes along the transect in a barplot
macroselev|>
  ggplot2::ggplot(mapping=aes(x=factor(elevation),y=CountSum,
                              group = year)) +
  geom_bar(stat = "identity") +
  facet_grid(year~.)  # what is in rows ~what is in columns
#  facet_wrap(~year)  # what is in rows ~what is in columns


# plot macros (y) versus elevation (x) in ggplot as a scatterplot, with each year as a different color
p1<- macroselev|> 
  ggplot2::ggplot(mapping=aes(x=elevation, y=CountSum, color=factor(year))) +
  geom_point(size=3) 
p1

# calculate the optimal preferred elevation by macros for each year using weighted.mean function
macroselev|> 
  dplyr::group_by(year) |>
  summarize(wa_elevation=weighted.mean(x=elevation,
       w=CountSum,na.rm=T)) #w = weighing factor)

## explore response to elevation and year as a linear model, call this m1
# first only elevation (not yet year)
macroselev
m1<-lm(CountSum~elevation, data=macroselev) # linear model of CountSum as a function of Location_ID
print(m1)
# so b0=-34.80 (intercept) and b1= 76.37 (coefficient or slope)
# so the model is y=-34.80 + 76.37x

#add the linear model to the plot
macroselev$pred1<-predict(m1) #calculate the predicted value of m1 for every observation
p2 <- p1 + geom_line(data=macroselev, aes(y=pred1), color="black", linewidth=1.2) # add the new predicted line to the previous plot p1)
print(p2)
#does not look like a good linear model

# fitmodel  m2 by adding a quadratic term for elevation_m to check for an ecological optimum
# y=b0 + b1x1 + b2*I(x^2) #expect that b2 is negative (gives us an optimum (berg parabool))
m2<-lm(CountSum~elevation + I(elevation^2), data=macroselev) # I() is a function that tells R to treat the term as a single variable so we can transform it
print(m2)

#add the linear model to the plot
# calculate the predicted value of m2 for every observation, add to the dataset as a variable as pred2
# add the new predicted line to the previous plot p2, store as object p3 and show it
macroselev$pred2<- predict(m2)
p2 + geom_line(data=macroselev, aes(y=pred2), color="red", linewidth=1.2)
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
m3<-lm(CountSum~elevation + I(elevation^2) + factor(year), data=macroselev) #add year as a factor
print(m3)
anova(m3) #elevation has an quadratic effect but also the year has an effect
#is there interaction (response of elevation to year)?
macroselev$pred3<-predict(m3)
p1 + geom_line(data=macroselev, aes(y=pred3, color=factor(year)), linewidth=1.2)

# does this significantly explain more variation than the model without year?
anova(m3,m2)

# yes this is a better model

# include the interaction between elevation_m and year (as a factor) in the model 
m4<-lm(CountSum~elevation + I(elevation^2) + factor(year) + elevation*factor(year), data=macroselev)
print(m4)

#is there a significant interaction between elevation and year?
macroselev$pred4<-predict(m4)
p1 + geom_line(data=macroselev, aes(y=pred4, color=factor(year)), linewidth=1.2)
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
macroselev
p1

#add the linear model to the plot
# calculate the predicted value of m2 for every observation, add to the dataset as a variable as pred2
# add the new predicted line to the previous plot p2, store as object p3 and show it
m5<-glm(CountSum~elevation,
        family = poisson(log),
        data=macroselev) #add year as a factor


# now test and show  the effect of both elevation , elevation squared and year
anova(m5,test="Chisq")

#add the linear model to the plot
# calculate the predicted value of m2 for every observation, add to the dataset as a variable as pred2
# add the new predicted line to the previous plot p2, store as object p3 and show it

macroselev$pred5<-predict(m5,type="response")
p1  + geom_line(data=macroselev,aes(y=pred5),linewidth=1.2)

# account for the ecological optimum
m6<-glm(CountSum~elevation+I(elevation^2),
        family = poisson(log),
        data=macroselev) #add year as a factor
macroselev$pred6<-predict(m6,type="response")
p1  + geom_line(data=macroselev,aes(y=pred6),linewidth=1.2)


# better than the previous?


# add the interaction to the model: elevation + elevation ^2 + year + elevation*year
# now test and show  the effect of both elevation + year


#add the  model to the plot
# calculate the predicted value of m2 for every observation, add to the dataset as a variable as pred2
# add the new predicted line to the previous plot p2, store as object p3 and show it





