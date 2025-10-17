# Compare the historical and current macroinvertebrate community structure in the Mara River using NMDS, DCA, and PCoA to determine the temporal changes in the community structure.
# clear everything in memory (of R)
remove(list=ls())

renv::restore()

library(vegan) 
library(psych) 
library(tidyverse)
library(ggpubr)


#database source
#browseURL("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/edit?usp=sharing")
# read the factbenthos historical and current data
combinedmacros<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1573424697&single=true&output=csv") |> 
tibble::column_to_rownames(var="Location_ID") # #names as variable to use in nmds otherwise it will not work, only works with numbers. 

#The NMDS plot will be used to determine the temporal changes in the macroinvertebrates community structure in the Mara river.

# Calculate the Jaccard dissimilarity matrix since it's a presence/absence data
set.seed(123)
jmacros <- vegdist(combinedmacros, method = "jaccard")# Jaccard dissimilarity matrix
print(jmacros)

# Perform NMDS using the Jaccard dissimilarity matrix
nmds_jaccard <- metaMDS(jmacros, k = 2, trace =T)
stressplot(nmds_jaccard)
#plot out our results and see if there is a difference between the historical and current
#Identify the time period as groups:
treat=as.matrix(c(rep("Historical",4),rep("Current",4)))

#Plot out the points (Location_ID):
 
ordiplot(nmds_jaccard,type="n",xlim=c(-.5,.5),ylim=c(-.5,.5))

## species scores not available
orditorp(nmds_jaccard,display="sites",col=c(rep("black",4),rep("brown",4)),air=0.01,cex=1.25)
legend(-.8,.5, c("Historical (2008-2009)","Current (2021-2023)"), cex=0.8,
       col=c("black","brown"), pch=15:15)

#Add  a dotted eclipse around the points to show the groups:

ordiellipse(nmds_jaccard, treat, display="si",lty=3,kind = "sd",conf = 0.75, col="black", show.groups="Historical")

ordiellipse(nmds_jaccard, treat, display="si", lty=3,kind = "sd",conf = 0.75, col="brown", show.groups="Current")

#NMDS plot is not the best for this data, so we will use DCA to determine the temporal changes in the macroinvertebrates community structure in the Mara river.

#####################################################################################################
#using PCoA to determine the temporal changes in the macroinvertebrates community structure in the Mara river.
#####################################################################################################
set.seed(123)
jmacros <- vegdist(combinedmacros, method = "jaccard")
print(jmacros)
#You are going to use the cmdscale function in the stats package to run the PCoA:
cmd<-cmdscale(jmacros, k=5, eig=TRUE) 
cmd

#The “points” are the coordinates of each sites. They are the eigenvectors scaled by the square root of their eigenvalues (i.e. the standard deviation):
cmd$points

#make a PCoA table to look at the eigenvalues, and the proportional and cumulative variance:
eigenvalues<-cmd$eig[1:5]
propVar<-eigenvalues/sum(eigenvalues)
cumVar<-cumsum(propVar)
PCoA_Table<-cbind(eigenvalues,propVar,cumVar)
PCoA_Table

#How many axes should you keep? You can use a scree plot to help you decide.
#Scree plot:
plot(eigenvalues)
lines(lowess(eigenvalues))

#The scree plot shows that the first two axes are the most important.

#Now, plot the first two PCoA axes:
x<-cmd$points[,1]
y<-cmd$points[,2]
plot(x,y,xlab= "PCoA1(72.9%)", ylab="PCoA2 (13.5%)", xlim=range(x)*1.2,ylim=range(y)*1.2, type="n")
text(x,y,labels=rownames(cmd$points), cex=.9)
# Add the title:
#title(main="Macroinvertebrate taxa composition", cex.main=1)

#add dotted lines to the inside of the plot to show the axes:
#abline(v=0, lty=3)
#abline(h=0, lty=3)

#use different colors for the historical and current data points using ecllipes
ordiellipse(cmd, treat, display="sites",lty=3,kind = "sd",conf = 0.75, col="black", show.groups="Historical")
ordiellipse(cmd, treat, display="sites", lty=3,kind = "sd",conf = 0.60, col="brown", show.groups="Current")

#leged 
legend(-0.3,0.3, c("2008-2009","2021-2023"), cex=0.8,
       col=c("black","brown"), pch=15:15)
####################################################################
#plot the figure without the Curr in 2021-2023 data
library(vegan)

set.seed(123)

# 1) Distance + Euclidean correction (prevents negative eigenvalues)
# Use binary=TRUE for Jaccard on presence/absence data
D <- vegdist(combinedmacros, method = "jaccard", binary = TRUE)

# Use wcmdscale(add=TRUE) -> Cailliez correction, returns $eig and $points
cmd <- wcmdscale(D, k = 5, eig = TRUE, add = TRUE)

# 2) Proportion of variation on the first two axes
# Use only positive eigenvalues for the denominator
eig_all <- cmd$eig
prop <- round(100 * eig_all[1:2] / sum(eig_all[eig_all > 0]), 1)

# (Optional) Table of eigenvalues, proportional and cumulative variance
eigenvalues <- eig_all[1:5]
propVar <- eigenvalues / sum(eig_all[eig_all > 0])
cumVar  <- cumsum(propVar)
PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
PCoA_Table

# 3) Groups (use the original rownames), and clean labels only for plotting
rn <- rownames(cmd$points)
treat <- factor(ifelse(grepl("Curr$", rn), "2021-2023", "2008-2009"),
                levels = c("2008-2009", "2021-2023"))
clean_labels <- gsub("Curr$", "", rn)

# 4) Plot
x <- cmd$points[, 1]; y <- cmd$points[, 2]
plot(x, y,
     xlab = paste0("PCoA1 (", prop[1], "%)"),
     ylab = paste0("PCoA2 (", prop[2], "%)"),
     xlim = range(x) * 1.2, ylim = range(y) * 1.5, type = "n")
text(x, y, labels = clean_labels, cex = 0.9)

# 5) Ellipses (pass coordinates + grouping)
# 5) Ellipses with thicker borders
ordiellipse(cmd$points[, 1:2], groups = treat,
            kind = "sd", conf = 0.75,
            col = c("black", "brown"),
            lty = 3, lwd = 2,   # <- thicker border line
            draw = "lines")


# 6) Legend with thicker ellipse borders
legend("topleft", bty = "n",
       legend = c("2008-2009", "2021-2023"),
       col    = c("black", "brown"),
       lty    = 3, 
       lwd    = 2,    # <- match thickness of ordiellipse lines
       pt.cex = 1)

#####################################################################################################
#PerMANOVA test to determine if there is a significant difference between the historical and current macroinvertebrate community structure in the Mara river.
#####################################################################################################
#next step is to determine if the historical and current sites are significantly different from each other using a permutation test.
#You will use the adonis function in the vegan package to run the permutation test:
set.seed(11) #set seed for reproducibility
permmacros<-adonis2(jmacros ~ treat, permutations=999)
permmacros

#The p-value is 0.024*, which means that the historical and current sites are significantly different from each other.
#Number of permutations: 999
#F-statistic: 15.369
#R-squared: 0.719


#####################################################################################################
#betadisper test 
#use the betadisper function in the vegan package :to test the differences betweeen the current sites
#betadisper test is used to determine the homogeneity of variance between the groups.
#####################################################################################################
set.seed(11) #set seed for reproducibility
#use the combined data-historical and current (presence, absence data)
# Euclidean distances between sites
dis<- vegdist(combinedmacros, "euclidean")

# groups are the 2 different years (historical and current)
groups <- as.factor(c(rep("Historical",4),rep("Current",4)))
#multivariate dispersions
MVdisp <- betadisper(dis, groups)
MVdisp

#Perform parametric test to determine if the multivariate dispersions are significantly different between the historical and current sites.
disp_aov<-anova(MVdisp)
disp_aov

# Tukey's Honest Significant Differences (HSD) test to determine which groups are significantly different from each other.
MVdisp.HSD <- TukeyHSD(MVdisp)
MVdisp.HSD

## non-parametric test: Permutation test for F-statistic to determine if the multivariate dispersions are significantly different between the historical and current sites.
perm_MVdisp <- permutest(MVdisp, permutations = 99, pairwise = TRUE)
perm_MVdisp


#####################################################################################################
#SIMPER analysis- (Similarity Percentage analysis)
# ----- Setup -----
remove(list = ls())

library(vegan)
library(dplyr)
library(ggplot2)
library(readr)

# ----- 1) Load data -----
data <- read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1573424697&single=true&output=csv",
  show_col_types = FALSE
) %>%
  mutate(
    Period = ifelse(grepl("Curr", Location_ID), "2021-2023", "2008-2009"),
    Period = factor(Period, levels = c("2008-2009","2021-2023"))
  )

# Replace NA with 0 (presence–absence)
data[is.na(data)] <- 0

# ----- 2) Species matrix & grouping factor -----
species <- data %>% select(-Location_ID, -Period)
# ensure strictly binary (optional)
# species[] <- ifelse(species > 0, 1, 0)

group <- data$Period

# sanity check: all numeric
stopifnot(all(sapply(species, is.numeric)))

# ----- 3) SIMPER with Jaccard -----
set.seed(123)
sim <- simper(
  comm = as.matrix(species),   # <-- key change (comm, not x)
  group = group,
  permutations = 999,
  method = "jaccard",
  binary = TRUE,
  trace = FALSE
)

sumry <- summary(sim, ordered = TRUE)
sim

# get the 2008-2009 vs 2021-2023 contrast name robustly
contrast_name <- names(sumry)[grepl("2008-2009.*2021-2023|2021-2023.*2008-2009", names(sumry))][1]
if (is.na(contrast_name)) contrast_name <- names(sumry)[1]

contrib <- as.data.frame(sumry[[contrast_name]])
contrib$Taxon <- rownames(contrib)

# ----- 4) Percent contribution & ordering -----
contrib <- contrib %>%
  mutate(pct = 100 * average / sum(average, na.rm = TRUE),
         cum_pct = cumsum(pct)) %>%
  arrange(desc(pct))

# ----- 5) Plot top N taxa -----
N <- 10
top_contrib <- contrib %>% slice_head(n = N)

ggplot(top_contrib, aes(x = reorder(Taxon, pct), y = pct)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Taxon",
    y = "Contribution to Dissimilarity (%)",
    title = paste0("Top ", N, " SIMPER Contributors (", contrast_name, ")")
  ) +
  theme_bw()

# Optional: view table
print(top_contrib[, c("Taxon","pct","cum_pct","average","sd","ava","avb","p")], row.names = FALSE)


###########################################################################################
# Historical macroinvertebrate taxa composition
####################################################################################################
#NMDS analysis 
#database source 

#browseURL("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/edit?usp=sharing")

#load data filter 2008,2009 and also group by Location_ID, month, year,River_reach and Family
historicalmacros<-readr::read_csv("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/pub?gid=1151562191&single=true&output=csv") |> 
  dplyr::filter(year %in% c(2008, 2009))|>
  tibble::column_to_rownames(var="Sample_ID")# convert distance_m to the row names of the tibble
head(historicalmacros)
historicalmacros2<-historicalmacros %>% select(-c(Location_ID, month, year, Reach)) #remove the columns that are not needed for the analysis
  head(historicalmacros2)

#PerMANOVA test to determine if the Location_ID, month, year, and River_reach are significant factors in explaining the variation in the macroinvertebrate community structure.
set.seed(123) #set seed for reproducibility
permmacros <- adonis2(historicalmacros2 ~month+ Reach, 
                      data = historicalmacros, 
                      method = "euclidean", permutations = 999,
                      by = "margin")
permmacros

#the PerMANOVA test for the historical data shows that the Location_ID, month, year, and River_reach are not significant factors in explaining the variation in the macroinvertebrate community structure.

#plot the NMDS ordination of the historical macroinvertebrate community structure
histomacros<-vegdist(historicalmacros2, "bray") #we use bray curtis distance matrix since it's count data
histomacros

#You are going to use the metaMDS function in the vegan package.
#K =2 because we are interested in only two dimension (which is common for NMDS)
#Trymax=100 because we have a small data set
nmdshisto<-metaMDS(histomacros,k=2, trace=T, trymax=100)
nmdshisto

stressplot(nmdshisto) #plot the stress plot
#What do the stress value and the fit (R2) of the monotonic regression tell you about the NMDS plot?
#The stress value is  0.092  which implies fairly good fit. The stress value is a measure of how well the data are represented in the NMDS plot.
#The closer the stress value is to 0, the better the representation.
#The R2 value is 0.992, which is also very good. The R2 value is a measure of how well the data are represented in the NMDS plot.
#The closer the R2 value is to 1, the better the representation.

#First, we need to create a grouping variable for the historical groups. #We will use the ca package to do this. Identify the river reach, month, year and site as groups: #River reach as a grouping variable
#First, we need to extract the site scores from the historical NMDS object.
#We will use the scores function in the vegan package to do this.
#We will also convert the site scores to a data frame.
histodata.scores <- as.data.frame(scores(nmdshisto)) #Using the scores function from vegan to extract the site scores and convert to a data.frame
histodata.scores$sites <- rownames(historicalmacros) # create a column of River reach from the original data frame macros
histodata.scores$Reach<-(historicalmacros$Reach) # create a column of 
head(histodata.scores)  #look at the data


#Plot NMDS for the current macroinvertebrate community structure

ggplot(histodata.scores, aes(x = NMDS1, y = NMDS2, col = Reach)) + 
  geom_point() +
  stat_ellipse(linetype = "dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-2.5, 2.5) +
  ylim(-2, 3) +
  labs() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 2, y = 3, label = paste("Stress = ", round(nmdshisto$stress, 3)))


####################################################################################################
#Current macroinvertebrate taxa composition
####################################################################################################
#NMDS alysis 
#database source 
#browseURL("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/edit?usp=sharing")

#load data filter 2021, 2022, 2023, and also group by Location_ID, month, year,River_reach and Family
currentmacros<-readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=798918621&single=true&output=csv")|> 
dplyr::filter(year %in% c(2021, 2022, 2023))|>
tibble::column_to_rownames(var="sample_ID")# convert distance_m to the row names of the tibble 
head(currentmacros)
currentmacros2<-currentmacros %>% select(-c(observation_ID, month, year, Reach)) #remove the columns that are not needed for the analysis
head(currentmacros2)

#PerMANOVA test to determine if the Location_ID, month, year, and River_reach are significant factors in explaining the variation in the macroinvertebrate community structure.
set.seed(123) #set seed for reproducibility
permcurrent <- adonis2(currentmacros2 ~ month + year + Reach, 
                      data = currentmacros, 
                      method = "euclidean", permutations = 999,
                      by = "margin")
permcurrent
###################################################################################################
#NMDS plot for the current macroinvertebrate community structure
#Since we have the count data, we will use the bray-curtis distance matrix
bmacros<-vegdist(currentmacros2, "bray") #we use bray curtis distance matrix since it's count data
bmacros

#You are going to use the metaMDS function in the vegan package.
#K =2 because we are interested in only two dimension (which is common for NMDS)
#Trymax=100 because we have a small data set
nmdsmacros<-metaMDS(bmacros,k=2, trace=T, trymax=100)
nmdsmacros

stressplot(nmdsmacros) #plot the stress plot
#What do the stress value and the fit (R2) of the monotonic regression tell you about the NMDS plot?
#The stress value is  0.166, which implies fairly good fit. The stress value is a measure of how well the data are represented in the NMDS plot.
#The closer the stress value is to 0, the better the representation.
#The R2 value is 0.972, which is also very good. The R2 value is a measure of how well the data are represented in the NMDS plot.
#The closer the R2 value is to 1, the better the representation.

#First, we need to create a grouping variable for the spatial-temporal macroinvertebrates groups. #We will use the ca package to do this. Identify the river reach, month, year and site as groups: #River reach as a grouping variable
#First, we need to extract the site scores from the NMDS object.
#We will use the scores function in the vegan package to do this.
#We will also convert the site scores to a data frame.
data.scores <- as.data.frame(scores(nmdsmacros)) #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$sites <- rownames(currentmacros) # create a column of River reach from the original data frame macros
data.scores$Reach <-(currentmacros$Reach) # create a column of 
head(data.scores)  #look at the data

#Plot NMDS for the current macroinvertebrate community structure

ggplot(data.scores, aes(x = NMDS1, y = NMDS2, col = Reach)) + 
  geom_point() +
  stat_ellipse(linetype = "dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-2.5, 2.5) +
  ylim(-2, 2) +
  labs() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 1.5, y = 2, label = paste("Stress = ", round(nmdsmacros$stress, 3)))

############################################################################
#the combined NMDS plot for the historical and current macroinvertebrate community structure
# Clear the environment
remove(list = ls())

# Load necessary libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(ggrepel)
# Load historical macroinvertebrate data
historicalmacros <- read_csv("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/pub?gid=1151562191&single=true&output=csv") %>% 
  filter(year %in% c(2008, 2009)) %>%
  column_to_rownames(var = "Sample_ID")

# Prepare historical data for NMDS
historicalmacros2 <- historicalmacros %>% select(-c(Location_ID, month, year, Reach))
histomacros <- vegdist(historicalmacros2, "bray")
nmdshisto <- metaMDS(histomacros, k = 2, trace = TRUE, trymax = 100)

# Extract site scores for historical NMDS
histodata.scores <- as.data.frame(scores(nmdshisto))
histodata.scores$sites <- rownames(historicalmacros)
histodata.scores$Reach <- historicalmacros$Reach
histodata.scores$Time_Period <- "2008-2009"  # Add time period label

# Load current macroinvertebrate data
currentmacros <- read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=798918621&single=true&output=csv") %>%
  filter(year %in% c(2021, 2022, 2023)) %>%
  column_to_rownames(var = "sample_ID")

# Prepare current data for NMDS
currentmacros2 <- currentmacros %>% select(-c(observation_ID, month, year, Reach))
bmacros <- vegdist(currentmacros2, "bray")
nmdsmacros <- metaMDS(bmacros, k = 2, trace = TRUE, trymax = 100)

# Extract site scores for current NMDS
data.scores <- as.data.frame(scores(nmdsmacros))
data.scores$sites <- rownames(currentmacros)
data.scores$Reach <- currentmacros$Reach
data.scores$Time_Period <- "2021-2023"  # Add time period label

# Combine historical and current NMDS data
combined_nmds <- bind_rows(histodata.scores, data.scores)

# Create a new data frame for stress values
stress_values <- data.frame(
  Time_Period = c("2008-2009", "2021-2023"),
  NMDS1 = -2,  # Position on x-axis
  NMDS2 = 2.5,  # Position on y-axis
  Stress = c("Stress = 0.092", "Stress = 0.166")
)

# Plot with stress values
ggplot(combined_nmds, aes(x = NMDS1, y = NMDS2, col = Reach)) + 
  geom_point(size = 3, alpha = 0.7, na.rm = TRUE) +  # Remove missing points
  stat_ellipse(linetype = "dashed", na.rm = TRUE) +  # Remove missing ellipses
  facet_wrap(~Time_Period, ncol = 2) +  # Split by time period
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),  # Simplifies facet labels
    strip.text = element_text(size = 12, face = "bold"),  # Adjust facet label size
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust margins
  ) +
  xlim(-3, 3) +  # Adjusted to include all data
  ylim(-3, 3) +  # Adjusted to include all data
  labs(x = "NMDS1", y = "NMDS2", col = "Reach") +
  # Add stress values using geom_text()
  geom_text(data = stress_values, aes(x = NMDS1, y = NMDS2, label = Stress), 
            inherit.aes = FALSE, size = 4, color = "black")
############################################################################
