# Create data for 2009
data_2009 <- data.frame(
  Location = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9"),
  Agriculture_Percentage = c(50.43, 55.45, 52.89, 65.43, 65.19, 64.07, 63.21, 58.60, 36.80),
  Forest_Percentage = c(31.55, 32.92, 37.85, 22.10, 21.87, 20.67, 20.24, 18.78, 10.61),
  Shrub_Woodland_Percentage = c(14.28, 9.53, 7.70, 9.28, 9.44, 9.84, 9.82, 9.01, 15.77),
  Grassland_Percentage = c(3.74, 2.10, 1.56, 3.19, 3.49, 5.40, 6.73, 13.60, 36.82)
)

# Create data for 2021-2023
data_2021_2023 <- data.frame(
  Location = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9"),
  Agriculture_Percentage = c(65.40, 67.61, 60.24, 72.83, 72.35, 71.86, 72.47, 66.66, 42.91),
  Forest_Percentage = c(31.20, 30.02, 39.04, 19.55, 19.06, 17.50, 16.22, 14.38, 7.84),
  Shrub_Woodland_Percentage = c(2.58, 1.82, 0.62, 2.99, 3.20, 3.37, 3.34, 3.31, 7.81),
  Grassland_Percentage = c(0.82, 0.56, 0.11, 4.63, 5.39, 7.26, 7.97, 15.65, 41.44)
)

# Merge the data for comparison
comparison <- merge(data_2009, data_2021_2023, by = "Location", suffixes = c("_2009", "_2021_2023"))

# Calculate the differences
comparison$Agriculture_Change <- comparison$Agriculture_Percentage_2021_2023 - comparison$Agriculture_Percentage_2009
comparison$Forest_Change <- comparison$Forest_Percentage_2021_2023 - comparison$Forest_Percentage_2009
comparison$Shrub_Woodland_Change <- comparison$Shrub_Woodland_Percentage_2021_2023 - comparison$Shrub_Woodland_Percentage_2009
comparison$Grassland_Change <- comparison$Grassland_Percentage_2021_2023 - comparison$Grassland_Percentage_2009

# View the comparison
print(comparison)

######################################
## --- Setup ----
library(dplyr)
library(tidyr)
library(car)        # for Levene's test (optional)
library(broom)      # tidy outputs

# Your data frames: data_2009 and data_2021_2023 (as provided)

## --- Tidy to long format and bind ----
to_long <- function(df, year_lab){
  df %>%
    pivot_longer(-Location, names_to = "CoverType", values_to = "Percent") %>%
    mutate(YearPeriod = year_lab)
}

long_2009       <- to_long(data_2009, "2009")
long_2021_2023  <- to_long(data_2021_2023, "2021_2023")

long_all <- bind_rows(long_2009, long_2021_2023) %>%
  mutate(
    Location   = factor(Location),
    YearPeriod = factor(YearPeriod, levels = c("2009","2021_2023")),
    CoverType  = factor(CoverType,
                        levels = c("Agriculture_Percentage","Forest_Percentage",
                                   "Shrub_Woodland_Percentage","Grassland_Percentage"),
                        labels = c("Agriculture","Forest","Shrub/Woodland","Grassland"))
  )

## --- (A) Repeated-measures ANOVA: one model per cover type ----
run_rm_aov <- function(cover){
  dat <- filter(long_all, CoverType == cover)
  # Repeated-measures ANOVA: Location as subject (error term),
  # YearPeriod is the within-subject factor (2 levels)
  fit <- aov(Percent ~ YearPeriod + Error(Location/YearPeriod), data = dat)
  list(cover = cover, anova = summary(fit))
}

results_aov <- lapply(levels(long_all$CoverType), run_rm_aov)

# Print ANOVA tables
for (res in results_aov) {
  cat("\n====", res$cover, "====\n")
  print(res$anova)
}

## --- (B) Paired t-tests (equivalent here) with effect size ----
paired_tests <- long_all %>%
  select(Location, CoverType, YearPeriod, Percent) %>%
  pivot_wider(names_from = YearPeriod, values_from = Percent) %>%
  group_by(CoverType) %>%
  summarise(
    t_test = list(t.test(`2021_2023`, `2009`, paired = TRUE))
  )

# Tidy print
for (i in seq_len(nrow(paired_tests))) {
  cover <- paired_tests$CoverType[i]
  tt    <- paired_tests$t_test[[i]]
  # Cohen's d for paired samples
  diffs <- with(tt$data, x - y)  # x=2021_2023, y=2009 inside t.test env
  d     <- mean(diffs) / sd(diffs)
  cat(sprintf("\n==== %s (paired t-test) ====\n", cover))
  cat(sprintf("mean change = %.2f\n", mean(diffs)))
  cat(sprintf("t(%d) = %.3f, p = %.4f\n", tt$parameter, tt$statistic, tt$p.value))
  cat(sprintf("Cohen's d (paired) = %.3f\n", d))
}

## --- Assumption checks (quick) ----
# 1) Normality of within-site differences (more appropriate than RM-ANOVA residuals here)
diffs_df <- long_all %>%
  select(Location, CoverType, YearPeriod, Percent) %>%
  pivot_wider(names_from = YearPeriod, values_from = Percent) %>%
  mutate(Diff = `2021_2023` - `2009`)

by(diffs_df$Diff, diffs_df$CoverType, shapiro.test)  # p > .05 supports normality

# 2) (Optional) Homogeneity is not required for paired tests, but if you still want
#    to look at residual spreads from a simple model:
#    Levene test on absolute residuals per period (illustrative)
leveneTest(Percent ~ YearPeriod, data = filter(long_all, CoverType == "Agriculture"))

