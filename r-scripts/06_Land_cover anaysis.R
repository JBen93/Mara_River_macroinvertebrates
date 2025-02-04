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
