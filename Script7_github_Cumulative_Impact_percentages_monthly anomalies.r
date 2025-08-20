# ------------------------------------------------------------#
#   Comparing methodological choices for environmental cumulative impacts analysis: The Black Sea as a case study
#   SCRIPT 6                                            
#   CUMULATIVE IMPACTS
#   Based on percentages of monthly anomalies
#   BLACK SEA
#   SURFACE DATA (Yearly matrix as input)
#   Authors: Lucia Espasandin, Fran Ramirez and Elena Lloret-Lloret
#   Last update: 20/08/2025                                                         
#------------------------------------------------------------ #

# AIM OF THE SCRIPT

# Calculate a cumulative impact analysis based on monthly data.
# We use a set of physical and a set of biogeochemical variables over a relatively large time range
# and calculate the percentage of anomaly change over time, per pixel.

# Here we use SURFACE environmental variables as an example, 
# but any of the previous pre-process outputs can be used instead (vertical integration weighted mean,vertical integration weighted sum, vertical integration 1m interpolation)



#------------------------------#
# 0. Set up ---- 
#------------------------------#

#------------------------------#
## 0. Load packages ---- 
#------------------------------#

# install.packages(c(
#   "raster", "fields", "plyr", "circular", "CircStats", "arrayhelpers",
#   "gstat", "ggplot2", "sp", "ncdf4", "reshape2", "terra", "sf",
#   "rstudioapi", "lubridate", "dplyr", "reshape", "nc", "tidyverse",
#   "R.matlab", "mapview", "rnaturalearth", "rnaturalearthdata", 
#   "abind", "plotly"
# ))

library(raster)
library(fields)
library(plyr)
library(circular)
library(CircStats)
library(arrayhelpers)
library(gstat)
library(ggplot2)
library(sp)
library(ncdf4)
library(reshape2)
library(terra)
library(sf)
library(rstudioapi)
library(lubridate)
library(dplyr)
library(reshape)
library(nc)
library(tidyverse)
library(R.matlab)
library(mapview)
library(rnaturalearth)
library(rnaturalearthdata)
library(abind)
library(plotly)

# ------------------------------ #
## 0.1. Set working directory -----
# ------------------------------- #
#this sets the working directory directly where your script is:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ----------------------------------------------------- #
## 0.2. Define parameters that need to be changed -----
# ---------------------------------------------------- #

# Determine the start and end year at the beginning of each script
# your data needs to have complete years (12 months)!!

start_year<- 1998
end_year<- 2021

years_all<- seq(start_year, end_year, 1)

# var time for linear models
length(years_all)
years_all 
years_all <- as.Date(paste0(years_all, "-01-01"))
years_all
times_all <- seq(1, length(years_all), by=1) #we create a time object to ease the script
times_all

# Define longitude and latitude ranges for the analysis
# If you don't need a subset, you still need to define the min and max of your entire dataset
# These values must be EXACTLY the SAME as in the nc. (otherwise it won't work)

# Lon lat for the Physical variables
lon_range_phy <- c(27.37000, 41.96259)  #all black sea
lat_range_phy <- c(40.86000, 46.80444)  #all black_sea

#Lon lat for the Biological variables
lon_range_bio <- c(27.250, 42.000) #all black sea
lat_range_bio <- c(40.500, 47.000) #all black sea

# Create colour palettes for layer use
blue_to_red_palette <- colorRampPalette(colors = c("blue", "beige", "red"))
magenta_palette <- colorRampPalette(colors = c('#3288bd','#32A8BD' , '#7EC9D6', '#93DACC','#fdae61','#9e0142'))
magenta_to_blue_palette <- colorRampPalette(colors = c('#3288bd','#32A8BD' ,'#fdae61','#9e0142'))

# --------------------------------------------- #
# 1. READ DATASETS and EXTRACT PARAMETERS ----
# --------------------------------------------- #

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --------------------------------------------- #
## 1.1. Load areas shapefile ----
# --------------------------------------------- #

# Define and load your shapefile area
# areas <- st_read("./Data/Shapefiles/Only_Black_Sea_CLEAN_good.shp")
# plot(areas)
# class(areas)
# crs(areas)
# str(areas)
# areas
# mapview(areas)

# Download administrative boundaries for your country
world <- ne_countries(scale = "medium", returnclass = "sf")
mapview(world)

World<- st_read("./Data/Shapefiles/World_Countries_Generalized/World_Countries_Generalized.shp")
mapview(World)       
plot(World)
class(World)
crs(World)
str(World)
World
mapview(World)
World_transformed <- st_transform(World, crs = 4326)  # EPSG: 4326 for WGS 84

EEZ<- st_read("./Data/Shapefiles/World_EEZ_v12_20231025/eez_boundaries_v12.shp")
mapview(EEZ) 
plot(EEZ)
class(EEZ)
crs(EEZ)
str(EEZ)
EEZ_transformed <- st_transform(EEZ, crs = 4326)  # EPSG: 4326 for WGS 84

#------------------------------#
## 1.2 ENVIRONMENTAL PRODUCTS ---- 
#------------------------------#

#------------------------------#
### 1.2.1. LOAD PHYSICAL FILES ----
#------------------------------# 
#Load
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

thetao <-readRDS( "./Data/CI_monthly_anomalies/surface_temp_matrix_Black_Sea_1993_2021.RData")
str(thetao)
summary(thetao)
dim(thetao)

# Sub-setting the third dimension (years 1998 to 2021)
thetao_subset <- thetao[, , 61:348]
dim(thetao_subset)

# Calculate and visualize the spatial mean of all the dataset
thetao_mean1998_2021 <- apply(thetao_subset,c(1,2),"mean")
image.plot(thetao_mean1998_2021, main="Mean surface Temp 98-21")

# Sea water potential temperature at sea floor: bottomT [°C]
bottomT<-readRDS( "./Data/CI_monthly_anomalies/temp_bot_matrix_Black_Sea_1993_2021.RData")
str(bottomT)
summary(bottomT)
dim(bottomT)

# Sub-setting the third dimension (years 1998 to 2021)
bottomT_subset <- bottomT[, , 61:348]
dim(bottomT_subset)

# Calculate and visualize the spatial mean of all the dataset
bottomT_mean1998_2021 <- apply(bottomT_subset,c(1,2),"mean")
image.plot(bottomT_mean1998_2021, main="Mean bottom Temp 98-21")

# Sea water salinity: so [10-3]
so<-readRDS( "./Data/CI_monthly_anomalies/surface_salinity_matrix_Black_Sea_1993_2021.RData")
str(so)
summary(so)
dim(so)

# Sub-setting the third dimension (years 1998 to 2021)
so_subset <- so[, , 61:348]
dim(so_subset)

# Calculate and visualize the spatial mean of all the dataset
so_mean1998_2021 <- apply(so_subset,c(1,2),"mean")
image.plot(so_mean1998_2021, main="Mean Surface Salinity 98-21")

#-------------------------------------------#
### 1.2.2. LOAD BIOGEOCHEMICAL FILES ----
#------------------------------------------#
# Chlorophyll
chl<-readRDS( "./Data/CI_monthly_anomalies/surface_chl_matrix_Black_Sea_1992_2022.RData")
str(chl)
summary(chl)
dim(chl) 

# Sub-setting the third dimension (years 1998 to 2021)
chl_subset <- chl[, , 73:360]
dim(chl_subset)

#Calculate and visualize the spatial mean of all the dataset
chl_mean1998_2021 <- apply(chl_subset,c(1,2),"mean")
image.plot(chl_mean1998_2021, main="Mean Surface Chl 1998-2021")

# NPP
nppv<-readRDS( "./Data/CI_monthly_anomalies/surface_npp_matrix_Black_Sea_1992_2022.RData")
str(nppv)
summary(nppv)
dim(nppv)

# Sub-setting the third dimension (years 1998 to 2021)
nppv_subset <- nppv[, , 73:360]
dim(nppv_subset)

# Calculate and visualize the spatial mean of all the dataset
nppv_mean1998_2021 <- apply(nppv_subset,c(1,2),"mean")
image.plot(nppv_mean1998_2021, main="Mean Surface nppv 1998-2021")

# Oxygen
o2<-readRDS( "./Data/CI_monthly_anomalies/surface_oxygen_matrix_Black_Sea_1992_2022.RData")
str(o2)
summary(o2)
dim(o2)

# Sub-setting the third dimension (years 1998 to 2021)
o2_subset <- o2[, , 73:360]
dim(o2_subset)

# Calculate and visualize the spatial mean of all the dataset
o2_mean1998_2021 <- apply(o2_subset,c(1,2),"mean")
image.plot(o2_mean1998_2021, main="Mean Surface oxygen 1998-2021")


#------------------------------#
# 2. CALCULATE CLIMATOLOGIES---- 
#------------------------------#	

## PHYSICAL VARIABLES

# Temperature climatological means
# This generates a vector all.months with month abbreviations ("Jan", "Feb", etc.) for each time point in thetao_subset
all.months <- rep(month.abb, length.out=dim(thetao_subset)[3])  # "Jan", "Feb", ..., for each time step
sel <- which(all.months == "Jan") # It identifies all time points that correspond to January and selects the data for those time points.
sel.inp.thetao <- thetao_subset[ , , sel]  # Extract data for January across all years (for every grid point)
clim_jan <- apply(sel.inp.thetao, c(1, 2), mean, na.rm=T)  # Average across the time dimension for January

# Initialize an empty array to store the climatology for each month (12 month
thetao_clim <- array(NA, dim = c(dim(thetao_subset)[1], dim(thetao_subset)[2], 12))
# This will hold the climatology for each month [longitude, latitude, month]
#i="Jan"
# Loop over each month
for(i in unique(all.months)){
  sel <- which(all.months == i)  # Select indices for the current month
  sel.inp <- thetao_subset[ , , sel]  # Extract data for that month across all years
  clim_temp <- apply(sel.inp, c(1, 2), mean, na.rm=T)  # Compute the mean across time
  # Store the climatology for the current month in the new array
  month_index <- match(i, month.abb)  # Get the index for the current month (1 for Jan, 2 for Feb, etc.)
  thetao_clim[ , , month_index] <- clim_temp
}

summary(thetao_clim)
str(thetao_clim)

# Verify output
# Extract climatology for January (or any other month)
image.plot(thetao_clim[,,1])  # January climatology

# Plot
# Calculate the mean ( of all the area) across all years for each month
thetao_climatology_mean <- apply(thetao_clim, c(3), mean, na.rm = TRUE)  # Averaging over the first two dimensions (116 and 128)

# Prepare individual yearly data
# Create a matrix to hold yearly averages
years <- seq(1998, 2021)
thetao_monthly_yearly_avg <- matrix(NA, nrow = length(years), ncol = 12)

# Loop through each year and calculate monthly averages
for (i in 1:length(years)) {
  # Calculate the starting and ending index for the current year in thetao_subset
  start_index <- (i - 1) * 12 + 1
  end_index <- i * 12  # This should give us 12 months for each year
  # Extract the monthly data for the current year
  thetao_monthly_data <- thetao_subset[, , start_index:end_index]  # Extract 3D array for this year
  
  # Calculate the mean for each month across all longitudes and latitudes
  for (month in 1:12) {
    thetao_monthly_yearly_avg[i, month] <- mean(thetao_monthly_data[, , month], na.rm = TRUE)
  }
}

# Convert to data frame for ggplot
thetao_monthly_yearly_avg_df <- melt(thetao_monthly_yearly_avg)
colnames(thetao_monthly_yearly_avg_df) <- c("Year", "Month", "Value")
thetao_monthly_yearly_avg_df$Year <- years  # Directly assign years here

# Prepare climatology data for plotting
thetao_climatology_df <- data.frame(Month = 1:12, Value = thetao_climatology_mean)

# Plot
ggplot() +
  geom_line(data = thetao_climatology_df, aes(x = Month, y = Value), color = '#d53e4f', size = 1.5) +
  geom_line(data = thetao_monthly_yearly_avg_df, aes(x = Month, y = Value, group = Year, color = as.factor(Year)), alpha = 0.3) +
  labs(title = "Monthly Surface Temperature Climatology (1998-2021)",
       x = "Month", y = "Temperature (°C)",
       color = "Year") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_classic()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/Monthly_surface_temperature_climatology_1998-2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()


#Bottom Temperature climatological means
# This generates a vector all.months with month abbreviations ("Jan", "Feb", etc.) for each time point in bottomT_subset
all.months <- rep(month.abb, length.out=dim(bottomT_subset)[3])  # "Jan", "Feb", ..., for each time step
sel <- which(all.months == "Jan") #It identifies all time points that correspond to January and selects the data for those time points.
sel.inp.bottomT <- bottomT_subset[ , , sel]  # Extract data for January across all years (for every grid point)
clim_jan <- apply(sel.inp.bottomT, c(1, 2), mean, na.rm=T)  # Average across the time dimension for January

# Initialize an empty array to store the climatology for each month (12 month
bottomT_clim <- array(NA, dim = c(dim(bottomT_subset)[1], dim(bottomT_subset)[2], 12))
# This will hold the climatology for each month [longitude, latitude, month]

# Loop over each month
for(i in unique(all.months)){
  sel <- which(all.months == i)  # Select indices for the current month
  sel.inp <- bottomT_subset[ , , sel]  # Extract data for that month across all years
  clim_temp <- apply(sel.inp, c(1, 2), mean, na.rm=T)  # Compute the mean across time
  # Store the climatology for the current month in the new array
  month_index <- match(i, month.abb)  # Get the index for the current month (1 for Jan, 2 for Feb, etc.)
  bottomT_clim[ , , month_index] <- clim_temp
}

summary(bottomT_clim)
str(bottomT_clim)

# Plot
# Calculate the mean across all years for each month
bottomT_climatology_mean <- apply(bottomT_clim, c(3), mean, na.rm = TRUE)  # Averaging over the first two dimensions (116 and 128)

# Prepare individual yearly data
# Create a matrix to hold yearly averages
years <- seq(1998, 2021)
bottomT_monthly_yearly_avg <- matrix(NA, nrow = length(years), ncol = 12)

# Loop through each year and calculate monthly averages
for (i in 1:length(years)) {
  # Calculate the starting and ending index for the current year in bottomT_subset
  start_index <- (i - 1) * 12 + 1
  end_index <- i * 12  # This should give us 12 months for each year
  # Extract the monthly data for the current year
  bottomT_monthly_data <- bottomT_subset[, , start_index:end_index]  # Extract 3D array for this year
  
  # Calculate the mean for each month across all longitudes and latitudes
  for (month in 1:12) {
    bottomT_monthly_yearly_avg[i, month] <- mean(bottomT_monthly_data[, , month], na.rm = TRUE)
  }
}

# Convert to data frame for ggplot
bottomT_monthly_yearly_avg_df <- melt(bottomT_monthly_yearly_avg)
colnames(bottomT_monthly_yearly_avg_df) <- c("Year", "Month", "Value")
bottomT_monthly_yearly_avg_df$Year <- years  # Directly assign years here

# Prepare climatology data for plotting
bottomT_climatology_df <- data.frame(Month = 1:12, Value = bottomT_climatology_mean)

# Plot
ggplot() +
  geom_line(data = bottomT_climatology_df, aes(x = Month, y = Value), color = "blue", size = 1.2, linetype = "dashed") +
  geom_line(data = bottomT_monthly_yearly_avg_df, aes(x = Month, y = Value, group = Year, color = as.factor(Year)), alpha = 0.5) +
  labs(title = "Monthly Bottom Temperature Climatology and Yearly Averages Black Sea (1998-2021)",
       x = "Month", y = "Temperature (°C)",
       color = "Year") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_minimal()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/Monthly_bottom_temperature_climatology_1998-2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

#Salinity climatological means
#This generates a vector all.months with month abbreviations ("Jan", "Feb", etc.) for each time point in so_subset
all.months <- rep(month.abb, length.out=dim(so_subset)[3])  # "Jan", "Feb", ..., for each time step
sel <- which(all.months == "Jan") #It identifies all time points that correspond to January and selects the data for those time points.
sel.inp.so <- so_subset[ , , sel]  # Extract data for January across all years (for every grid point)
clim_jan <- apply(sel.inp.so, c(1, 2), mean, na.rm=T)  # Average across the time dimension for January

# Initialize an empty array to store the climatology for each month (12 month
so_clim <- array(NA, dim = c(dim(so_subset)[1], dim(so_subset)[2], 12))
# This will hold the climatology for each month [longitude, latitude, month]

# Loop over each month
for(i in unique(all.months)){
  sel <- which(all.months == i)  # Select indices for the current month
  sel.inp <- so_subset[ , , sel]  # Extract data for that month across all years
  clim_temp <- apply(sel.inp, c(1, 2), mean, na.rm=T)  # Compute the mean across time
  # Store the climatology for the current month in the new array
  month_index <- match(i, month.abb)  # Get the index for the current month (1 for Jan, 2 for Feb, etc.)
  so_clim[ , , month_index] <- clim_temp
}

summary(so_clim)
str(so_clim)

# Plot
# Calculate the mean across all years for each month
so_climatology_mean <- apply(so_clim, c(3), mean, na.rm = TRUE)  # Averaging over the first two dimensions (116 and 128)

# Prepare individual yearly data
# Create a matrix to hold yearly averages
years <- seq(1998, 2021)
so_monthly_yearly_avg <- matrix(NA, nrow = length(years), ncol = 12)

# Loop through each year and calculate monthly averages
for (i in 1:length(years)) {
  # Calculate the starting and ending index for the current year in so_subset
  start_index <- (i - 1) * 12 + 1
  end_index <- i * 12  # This should give us 12 months for each year
  # Extract the monthly data for the current year
  so_monthly_data <- so_subset[, , start_index:end_index]  # Extract 3D array for this year
  
  # Calculate the mean for each month across all longitudes and latitudes
  for (month in 1:12) {
    so_monthly_yearly_avg[i, month] <- mean(so_monthly_data[, , month], na.rm = TRUE)
  }
}

# Convert to data frame for ggplot
so_monthly_yearly_avg_df <- melt(so_monthly_yearly_avg)
colnames(so_monthly_yearly_avg_df) <- c("Year", "Month", "Value")
so_monthly_yearly_avg_df$Year <- years  # Directly assign years here

# Prepare climatology data for plotting
so_climatology_df <- data.frame(Month = 1:12, Value = so_climatology_mean)

# Plot
ggplot() +
  geom_line(data = so_climatology_df, aes(x = Month, y = Value), color = "blue", size = 1.2, linetype = "dashed") +
  geom_line(data = so_monthly_yearly_avg_df, aes(x = Month, y = Value, group = Year, color = as.factor(Year)), alpha = 0.5) +
  labs(title = "Monthly surface Salinity Climatology and Yearly Averages Black Sea (1998-2021)",
       x = "Month", y = "Salinity (PSU)",
       color = "Year") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_minimal()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/Monthly_surface_salinity_climatology_1998-2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

## BIOGEOCHEMICAL VARIABLES

#  Chlorophyll climatological means
#This generates a vector all.months with month abbreviations ("Jan", "Feb", etc.) for each time point in chl_subset
all.months <- rep(month.abb, length.out=dim(chl_subset)[3])  # "Jan", "Feb", ..., for each time step
sel <- which(all.months == "Jan") #It identifies all time points that correspond to January and selects the data for those time points.
sel.inp.chl <- chl_subset[ , , sel]  # Extract data for January across all years (for every grid point)
clim_jan <- apply(sel.inp.chl, c(1, 2), mean, na.rm=T)  # Average across the time dimension for January


# Initialize an empty array to store the climatology for each month (12 month
chl_clim <- array(NA, dim = c(dim(chl_subset)[1], dim(chl_subset)[2], 12))
# This will hold the climatology for each month [longitude, latitude, month]

# Loop over each month
for(i in unique(all.months)){
  sel <- which(all.months == i)  # Select indices for the current month
  sel.inp <- chl_subset[ , , sel]  # Extract data for that month across all years
  clim_temp <- apply(sel.inp, c(1, 2), mean, na.rm=T)  # Compute the mean across time
  # Store the climatology for the current month in the new array
  month_index <- match(i, month.abb)  # Get the index for the current month (1 for Jan, 2 for Feb, etc.)
  chl_clim[ , , month_index] <- clim_temp
}

summary(chl_clim)
str(chl_clim)

# Plot
# Calculate the mean across all years for each month
chl_climatology_mean <- apply(chl_clim, c(3), mean, na.rm = TRUE)  # Averaging over the first two dimensions (116 and 128)

# Prepare individual yearly data
# Create a matrix to hold yearly averages
years <- seq(1998, 2021)
chl_monthly_yearly_avg <- matrix(NA, nrow = length(years), ncol = 12)

# Loop through each year and calculate monthly averages
for (i in 1:length(years)) {
  # Calculate the starting and ending index for the current year in chl_subset
  start_index <- (i - 1) * 12 + 1
  end_index <- i * 12  # This should give us 12 months for each year
  # Extract the monthly data for the current year
  chl_monthly_data <- chl_subset[, , start_index:end_index]  # Extract 3D array for this year
  
  # Calculate the mean for each month across all longitudes and latitudes
  for (month in 1:12) {
    chl_monthly_yearly_avg[i, month] <- mean(chl_monthly_data[, , month], na.rm = TRUE)
  }
}

# Convert to data frame for ggplot
chl_monthly_yearly_avg_df <- melt(chl_monthly_yearly_avg)
colnames(chl_monthly_yearly_avg_df) <- c("Year", "Month", "Value")
chl_monthly_yearly_avg_df$Year <- years  # Directly assign years here

# Prepare climatology data for plotting
chl_climatology_df <- data.frame(Month = 1:12, Value = chl_climatology_mean)

# Plot
ggplot() +
  geom_line(data = chl_climatology_df, aes(x = Month, y = Value), color = "blue", size = 1.2, linetype = "dashed") +
  geom_line(data = chl_monthly_yearly_avg_df, aes(x = Month, y = Value, group = Year, color = as.factor(Year)), alpha = 0.5) +
  labs(title = "Monthly surface Chlorophyll Climatology and Yearly Averages Balck Sea (1998-2021)",
       x = "Month", y = "Chl (mg/m3)",
       color = "Year") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_minimal()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/Monthly_surface_chl_climatology_1998-2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()


# NPP climatological means
# This generates a vector all.months with month abbreviations ("Jan", "Feb", etc.) for each time point in nppv_subset
all.months <- rep(month.abb, length.out=dim(nppv_subset)[3])  # "Jan", "Feb", ..., for each time step
sel <- which(all.months == "Jan") #It identifies all time points that correspond to January and selects the data for those time points.
sel.inp.nppv <- nppv_subset[ , , sel]  # Extract data for January across all years (for every grid point)
clim_jan <- apply(sel.inp.nppv, c(1, 2), mean, na.rm=T)  # Average across the time dimension for January

# Initialize an empty array to store the climatology for each month (12 month
nppv_clim <- array(NA, dim = c(dim(nppv_subset)[1], dim(nppv_subset)[2], 12))
# This will hold the climatology for each month [longitude, latitude, month]

# Loop over each month
for(i in unique(all.months)){
  sel <- which(all.months == i)  # Select indices for the current month
  sel.inp <- nppv_subset[ , , sel]  # Extract data for that month across all years
  clim_temp <- apply(sel.inp, c(1, 2), mean, na.rm=T)  # Compute the mean across time
  # Store the climatology for the current month in the new array
  month_index <- match(i, month.abb)  # Get the index for the current month (1 for Jan, 2 for Feb, etc.)
  nppv_clim[ , , month_index] <- clim_temp
}

summary(nppv_clim)
str(nppv_clim)

# Plot
# Calculate the mean across all years for each month
nppv_climatology_mean <- apply(nppv_clim, c(3), mean, na.rm = TRUE)  # Averaging over the first two dimensions (116 and 128)

# Prepare individual yearly data
# Create a matrix to hold yearly averages
years <- seq(1998, 2021)
nppv_monthly_yearly_avg <- matrix(NA, nrow = length(years), ncol = 12)

# Loop through each year and calculate monthly averages
for (i in 1:length(years)) {
  # Calculate the starting and ending index for the current year in nppv_subset
  start_index <- (i - 1) * 12 + 1
  end_index <- i * 12  # This should give us 12 months for each year
  # Extract the monthly data for the current year
  nppv_monthly_data <- nppv_subset[, , start_index:end_index]  # Extract 3D array for this year
  
  # Calculate the mean for each month across all longitudes and latitudes
  for (month in 1:12) {
    nppv_monthly_yearly_avg[i, month] <- mean(nppv_monthly_data[, , month], na.rm = TRUE)
  }
}

# Convert to data frame for ggplot
nppv_monthly_yearly_avg_df <- melt(nppv_monthly_yearly_avg)
colnames(nppv_monthly_yearly_avg_df) <- c("Year", "Month", "Value")
nppv_monthly_yearly_avg_df$Year <- years  # Directly assign years here

# Prepare climatology data for plotting
nppv_climatology_df <- data.frame(Month = 1:12, Value = nppv_climatology_mean)

# Plot
ggplot() +
  geom_line(data = nppv_climatology_df, aes(x = Month, y = Value), color = "blue", size = 1.2, linetype = "dashed") +
  geom_line(data = nppv_monthly_yearly_avg_df, aes(x = Month, y = Value, group = Year, color = as.factor(Year)), alpha = 0.5) +
  labs(title = "Monthly surface NPP Climatology and Yearly Averages Black Sea (1998-2021)",
       x = "Month", y = "NPP",
       color = "Year") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_minimal()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/Monthly_surface_NPP_climatology_1998-2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

# Oxygen climatological means
# This generates a vector all.months with month abbreviations ("Jan", "Feb", etc.) for each time point in o2_subset
all.months <- rep(month.abb, length.out=dim(o2_subset)[3])  # "Jan", "Feb", ..., for each time step
sel <- which(all.months == "Jan") #It identifies all time points that correspond to January and selects the data for those time points.
sel.inp.o2 <- o2_subset[ , , sel]  # Extract data for January across all years (for every grid point)
clim_jan <- apply(sel.inp.o2, c(1, 2), mean, na.rm=T)  # Average across the time dimension for January

# Initialize an empty array to store the climatology for each month (12 month
o2_clim <- array(NA, dim = c(dim(o2_subset)[1], dim(o2_subset)[2], 12))
# This will hold the climatology for each month [longitude, latitude, month]

# Loop over each month
for(i in unique(all.months)){
  sel <- which(all.months == i)  # Select indices for the current month
  sel.inp <- o2_subset[ , , sel]  # Extract data for that month across all years
  clim_temp <- apply(sel.inp, c(1, 2), mean, na.rm=T)  # Compute the mean across time
  # Store the climatology for the current month in the new array
  month_index <- match(i, month.abb)  # Get the index for the current month (1 for Jan, 2 for Feb, etc.)
  o2_clim[ , , month_index] <- clim_temp
}

summary(o2_clim)
str(o2_clim)

# Plot
# Calculate the mean across all years for each month
o2_climatology_mean <- apply(o2_clim, c(3), mean, na.rm = TRUE)  # Averaging over the first two dimensions (116 and 128)

# Prepare individual yearly data
# Create a matrix to hold yearly averages
years <- seq(1998, 2021)
o2_monthly_yearly_avg <- matrix(NA, nrow = length(years), ncol = 12)

# Loop through each year and calculate monthly averages
for (i in 1:length(years)) {
  # Calculate the starting and ending index for the current year in o2_subset
  start_index <- (i - 1) * 12 + 1
  end_index <- i * 12  # This should give us 12 months for each year
  # Extract the monthly data for the current year
  o2_monthly_data <- o2_subset[, , start_index:end_index]  # Extract 3D array for this year
  
  # Calculate the mean for each month across all longitudes and latitudes
  for (month in 1:12) {
    o2_monthly_yearly_avg[i, month] <- mean(o2_monthly_data[, , month], na.rm = TRUE)
  }
}

# Convert to data frame for ggplot
o2_monthly_yearly_avg_df <- melt(o2_monthly_yearly_avg)
colnames(o2_monthly_yearly_avg_df) <- c("Year", "Month", "Value")
o2_monthly_yearly_avg_df$Year <- years  # Directly assign years here

# Prepare climatology data for plotting
o2_climatology_df <- data.frame(Month = 1:12, Value = o2_climatology_mean)

# Plot
ggplot() +
  geom_line(data = o2_climatology_df, aes(x = Month, y = Value), color = "blue", size = 1.2, linetype = "dashed") +
  geom_line(data = o2_monthly_yearly_avg_df, aes(x = Month, y = Value, group = Year, color = as.factor(Year)), alpha = 0.5) +
  labs(title = "Monthly surface Climatology and Yearly Averages Black Sea (1998-2021)",
       x = "Month", y = "Oxygen",
       color = "Year") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_minimal()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/Monthly_surface_O2_climatology_1998-2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

#-------------------------------------#
# 2.1. CALCULATE THE ANOMALIES %  ---- 
#-------------------------------------#	
# Compute anomalies as PERCENTAGES

#PHYSICAL VARIABLES

#Temperature
# Empty array to store the anomalies
thetao_anom <- array(NA, dim = dim(thetao_subset))  # Same dimensions as thetao_subset

for(i in 1:dim(thetao_subset)[3]) {
  print(i)
  # Get the current month corresponding to the time step i
  current_month <- all.months[i]
  # Find the corresponding climatology for the current month
  month_index <- match(current_month, month.abb)  # Get index for current month (1 for Jan, 2 for Feb, ..., 12 for Dec)
  # Calculate the anomaly as a percentage relative to the climatology
  thetao_anom[ , , i] <- (thetao_subset[ , , i] - thetao_clim[ , , month_index]) / thetao_clim[ , , month_index] * 100
}

# Verify output
# Extract anomalies for January 1993
image.plot(thetao_anom[,,1])  # Assuming the first time step is January 1993

# Plot the anomalies
# Assuming thetao_anom is [lon, lat, time], calculate the mean anomalies for each time step
thetao_anom_avg <- apply(thetao_anom, 3, mean, na.rm = TRUE)  # Averaging spatially for each time step

# Create a data frame for anomalies
thetao_anomalies_df <- data.frame(
  Time = seq.Date(from = as.Date("1998-01-01"), by = "month", length.out = dim(thetao_anom)[3]),  # Time series from Jan 1993
  Anomalies = thetao_anom_avg,  # The averaged anomalies
  Month = factor(rep(month.abb, length.out = dim(thetao_anom)[3]), levels = month.abb)  # Repeat month labels for 288 months
)
# Plotting the anomalies
ggplot() +
  geom_point(data = thetao_anomalies_df, aes(x = Month, y = Anomalies, group = Time, color = Time), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
  labs(title = "Surface Temperature Monthly Climatology and Anomalies Black Sea",
       x = "Month", 
       y = "Temperature / Anomalies (%)",
       color = "Year") +
  theme_minimal() +
  theme(legend.position = "bottom")

#and
ggplot(thetao_anomalies_df, aes(x = Time, y = Anomalies)) +
  geom_line(color = '#d53e4f', size= 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black",linetype="dashed") +  # Add a linear trend line
  labs(title = "Average monthly surface Temperature Anomalies",
       x = "Month",
       y = "Anomaly (%)") +
  theme_classic()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/Average_monthly_surface_Temperature_Anomalies_2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

# Bottom Temperature
# Empty array to store the anomalies
bottomT_anom <- array(NA, dim = dim(bottomT_subset))  # Same dimensions as bottomT_subset

for(i in 1:dim(bottomT_subset)[3]) {
  print(i)
  # Get the current month corresponding to the time step i
  current_month <- all.months[i]
  # Find the corresponding climatology for the current month
  month_index <- match(current_month, month.abb)  # Get index for current month (1 for Jan, 2 for Feb, ..., 12 for Dec)
  # Calculate the anomaly as a percentage relative to the climatology
  bottomT_anom[ , , i] <- (bottomT_subset[ , , i] - bottomT_clim[ , , month_index]) / bottomT_clim[ , , month_index] * 100
}

# Verify output
# Extract anomalies for January 1998
image.plot(bottomT_anom[,,1])  # Assuming the first time step is January 1998

# Plot the anomalies
# Assuming bottomT_anom is [lon, lat, time], calculate the mean anomalies for each time step
bottomT_anom_avg <- apply(bottomT_anom, 3, mean, na.rm = TRUE)  # Averaging spatially for each time step

# Create a data frame for anomalies
bottomT_anomalies_df <- data.frame(
  Time = seq.Date(from = as.Date("1998-01-01"), by = "month", length.out = dim(bottomT_anom)[3]),  # Time series from Jan 1993
  Anomalies = bottomT_anom_avg,  # The averaged anomalies
  Month = factor(rep(month.abb, length.out = dim(bottomT_anom)[3]), levels = month.abb)  # Repeat month labels for 288 months
)

# Plotting the anomalies
ggplot() +
  geom_point(data = bottomT_anomalies_df, aes(x = Month, y = Anomalies, group = Time, color = Time), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
  labs(title = "Bottom Temperature Monthly Climatology and Anomalies Black Sea",
       x = "Month", 
       y = "Temperature / Anomalies (%)",
       color = "Year") +
  theme_minimal() +
  theme(legend.position = "bottom")

#and
ggplot(bottomT_anomalies_df, aes(x = Time, y = Anomalies)) +
  geom_line(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a linear trend line
  labs(title = "Average Monthly Bottom Temperature Anomalies Black Sea (1998-2021)",
       x = "Month",
       y = "Anomaly (%)") +
  theme_minimal()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/Average_monthly_bottom_Temperature_Anomalies_2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

# Salinity
# Empty array to store the anomalies
so_anom <- array(NA, dim = dim(so_subset))  # Same dimensions as so_subset

for(i in 1:dim(so_subset)[3]) {
  print(i)
  # Get the current month corresponding to the time step i
  current_month <- all.months[i]
  # Find the corresponding climatology for the current month
  month_index <- match(current_month, month.abb)  # Get index for current month (1 for Jan, 2 for Feb, ..., 12 for Dec)
  # Calculate the anomaly as a percentage relative to the climatology
  so_anom[ , , i] <- (so_subset[ , , i] - so_clim[ , , month_index]) / so_clim[ , , month_index] * 100
}

# Verify output
# Extract anomalies for January 1998
image.plot(so_anom[,,1])  # Assuming the first time step is January 1998

# Plot the anomalies
# Assuming so_anom is [lon, lat, time], calculate the mean anomalies for each time step
so_anom_avg <- apply(so_anom, 3, mean, na.rm = TRUE)  # Averaging spatially for each time step

# Create a data frame for anomalies
so_anomalies_df <- data.frame(
  Time = seq.Date(from = as.Date("1998-01-01"), by = "month", length.out = dim(so_anom)[3]),  # Time series from Jan 1993
  Anomalies = so_anom_avg,  # The averaged anomalies
  Month = factor(rep(month.abb, length.out = dim(so_anom)[3]), levels = month.abb)  # Repeat month labels for 288 months
)

# Plotting the anomalies
ggplot() +
  geom_point(data = so_anomalies_df, aes(x = Month, y = Anomalies, group = Time, color = Time), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
  labs(title = "Surface Salinity Monthly Climatology and Anomalies Black Sea",
       x = "Month", 
       y = "Salinity / Anomalies (%)",
       color = "Year") +
  theme_minimal() +
  theme(legend.position = "bottom")

#and
ggplot(so_anomalies_df, aes(x = Time, y = Anomalies)) +
  geom_line(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a linear trend line
  labs(title = "Average Monthly surface Salinity Anomalies Black Sea (1998-2021)",
       x = "Month",
       y = "Anomaly (%)") +
  theme_minimal()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/Average_monthly_surface_salinity_Anomalies_2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

## BIOGEOCHEMICAL VARIABLES

# Chlorophyll
# Empty array to store the anomalies
chl_anom <- array(NA, dim = dim(chl_subset))  # Same dimensions as chl_subset

for(i in 1:dim(chl_subset)[3]) {
  print(i)
  # Get the current month corresponding to the time step i
  current_month <- all.months[i]
  # Find the corresponding climatology for the current month
  month_index <- match(current_month, month.abb)  # Get index for current month (1 for Jan, 2 for Feb, ..., 12 for Dec)
  # Calculate the anomaly as a percentage relative to the climatology
  chl_anom[ , , i] <- (chl_subset[ , , i] - chl_clim[ , , month_index]) / chl_clim[ , , month_index] * 100
}

# Verify output
# Extract anomalies for January 1998
image.plot(chl_anom[,,1])  # Assuming the first time step is January 1998

# Plot the anomalies
# Assuming chl_anom is [lon, lat, time], calculate the mean anomalies for each time step
chl_anom_avg <- apply(chl_anom, 3, mean, na.rm = TRUE)  # Averaging spatially for each time step

# Create a data frame for anomalies
chl_anomalies_df <- data.frame(
  Time = seq.Date(from = as.Date("1998-01-01"), by = "month", length.out = dim(chl_anom)[3]),  # Time series from Jan 1993
  Anomalies = chl_anom_avg,  # The averaged anomalies
  Month = factor(rep(month.abb, length.out = dim(chl_anom)[3]), levels = month.abb)  # Repeat month labels for 288 months
)

# Plotting the anomalies
ggplot() +
  geom_point(data = chl_anomalies_df, aes(x = Month, y = Anomalies, group = Time, color = Time), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
  labs(title = "Surface Chlorophyll Monthly Climatology and Anomalies Black Sea",
       x = "Month", 
       y = "Chl / Anomalies (%)",
       color = "Year") +
  theme_minimal() +
  theme(legend.position = "bottom")

#and
ggplot(chl_anomalies_df, aes(x = Time, y = Anomalies)) +
  geom_line(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a linear trend line
  labs(title = "Average Monthly Surface Chlorophyll Anomalies Black Sea (1998-2021)",
       x = "Month",
       y = "Anomaly (%)") +
  theme_minimal()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/Average_monthly_surface_chl_Anomalies_2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

# NPP
# Empty array to store the anomalies
nppv_anom <- array(NA, dim = dim(nppv_subset))  # Same dimensions as nppv_subset

for(i in 1:dim(nppv_subset)[3]) {
  print(i)
  # Get the current month corresponding to the time step i
  current_month <- all.months[i]
  # Find the corresponding climatology for the current month
  month_index <- match(current_month, month.abb)  # Get index for current month (1 for Jan, 2 for Feb, ..., 12 for Dec)
  # Calculate the anomaly as a percentage relative to the climatology
  nppv_anom[ , , i] <- (nppv_subset[ , , i] - nppv_clim[ , , month_index]) / nppv_clim[ , , month_index] * 100
}

# Verify output
# Extract anomalies for January 1998
image.plot(nppv_anom[,,1])  # Assuming the first time step is January 1998

# Plot the anomalies
# Assuming nppv_anom is [lon, lat, time], calculate the mean anomalies for each time step
nppv_anom_avg <- apply(nppv_anom, 3, mean, na.rm = TRUE)  # Averaging spatially for each time step

# Create a data frame for anomalies
nppv_anomalies_df <- data.frame(
  Time = seq.Date(from = as.Date("1998-01-01"), by = "month", length.out = dim(nppv_anom)[3]),  # Time series from Jan 1993
  Anomalies = nppv_anom_avg,  # The averaged anomalies
  Month = factor(rep(month.abb, length.out = dim(nppv_anom)[3]), levels = month.abb)  # Repeat month labels for 288 months
)

# Plotting the anomalies
ggplot() +
  geom_point(data = nppv_anomalies_df, aes(x = Month, y = Anomalies, group = Time, color = Time), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
  labs(title = "Surface NPPV Monthly Climatology and Anomalies Black Sea",
       x = "Month", 
       y = "NPP / Anomalies (%)",
       color = "Year") +
  theme_minimal() +
  theme(legend.position = "bottom")

#and
ggplot(nppv_anomalies_df, aes(x = Time, y = Anomalies)) +
  geom_line(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a linear trend line
  labs(title = "Average Monthly surface NPP Anomalies Black Sea (1998-2021)",
       x = "Month",
       y = "Anomaly (%)") +
  theme_minimal()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/Average_monthly_surface_NPP_Anomalies_2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

# Oxygen
# Empty array to store the anomalies
o2_anom <- array(NA, dim = dim(o2_subset))  # Same dimensions as o2_subset

for(i in 1:dim(o2_subset)[3]) {
  print(i)
  # Get the current month corresponding to the time step i
  current_month <- all.months[i]
  # Find the corresponding climatology for the current month
  month_index <- match(current_month, month.abb)  # Get index for current month (1 for Jan, 2 for Feb, ..., 12 for Dec)
  # Calculate the anomaly as a percentage relative to the climatology
  o2_anom[ , , i] <- (o2_subset[ , , i] - o2_clim[ , , month_index]) / o2_clim[ , , month_index] * 100
}

# Verify output
# Extract anomalies for January 1998
image.plot(o2_anom[,,1])  # Assuming the first time step is January 1998

# Plot the anomalies
# Assuming o2_anom is [lon, lat, time], calculate the mean anomalies for each time step
o2_anom_avg <- apply(o2_anom, 3, mean, na.rm = TRUE)  # Averaging spatially for each time step

# Create a data frame for anomalies
o2_anomalies_df <- data.frame(
  Time = seq.Date(from = as.Date("1998-01-01"), by = "month", length.out = dim(o2_anom)[3]),  # Time series from Jan 1993
  Anomalies = o2_anom_avg,  # The averaged anomalies
  Month = factor(rep(month.abb, length.out = dim(o2_anom)[3]), levels = month.abb)  # Repeat month labels for 288 months
)

# Plotting the anomalies
ggplot() +
  geom_point(data = o2_anomalies_df, aes(x = Month, y = Anomalies, group = Time, color = Time), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
  labs(title = "Surface Oxygen Monthly Climatology and Anomalies Black Sea",
       x = "Month", 
       y = "Oxygen/ Anomalies (%)",
       color = "Year") +
  theme_minimal() +
  theme(legend.position = "bottom")

#and
ggplot(o2_anomalies_df, aes(x = Time, y = Anomalies)) +
  geom_line(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a linear trend line
  labs(title = "Average Monthly surface Oxygen Anomalies Black Sea(1998-2021)",
       x = "Month",
       y = "Anomaly (%)") +
  theme_minimal()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/Average_monthly_surface_O2_Anomalies_2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

#------------------------------#
# 3. CALCULATE THE SLOPE    ---- 
#------------------------------#	

#-----------------------------------#
## 3.1. Slope for PHYSICAL FILES ----
#-----------------------------------#

# Define the function to calculate the annual rate of change
f.rata <- function(x) {
  if (sum(!is.na(x)) < 100) {
    rata <- NA  # If less than 100 valid values, don't compute
  } else {
    ind <- 1:length(x)
    mdl <- lm(x ~ ind, na.action = na.exclude)  # Linear model
    slp <- mdl$coefficients[2]  # Get slope
    intr <- mdl$coefficients[1]  # Get intercept
    #rata <- slp * 12 # alternative added my me
    
    # Compute rate as evolution per YEAR
    first.intercept <- intr + slp  # Intercept at the first index
    last.ind <- ind[length(ind)]  # Last index
    fin.intercept <- intr + slp * last.ind  # Intercept at the last index
    rata <- (fin.intercept - first.intercept) / (last.ind-1) * 12  # Annual rate #(value at last time−value at first time)/total times *12
    }
  return(rata)
}
# rata <- rata <- slp * 12
#----------------------------------#
### 3.1.1. Apply the function ----
#----------------------------------#
# Applying fun2 (slope calculation) 
# In apply(): c(1,2) indicates that the calculation is performed for both rows and columns
# It calculates the regression of the value over the 3rd dimension which is time
thetao_mean1998_2021_slopes <- apply(thetao_anom,c(1,2),"f.rata")
bottomT_mean1998_2021_slopes <- apply(bottomT_anom,c(1,2),"f.rata")
so_mean1998_2021_slopes <- apply(so_anom,c(1,2),"f.rata")

#------------------------------------#  
### 3.1.2. Visualization of arrays ----
#-------------------------------------#
# You can visualize arrays with image.plot()
par(mfrow=c(2,2))
image.plot(thetao_mean1998_2021_slopes, col = magenta_to_blue_palette(100), main="Surface Temp 98-21")
image.plot(bottomT_mean1998_2021_slopes,col = magenta_to_blue_palette(100), main="Bottom Temp 98-21")
image.plot(so_mean1998_2021_slopes,col = magenta_to_blue_palette(100), main="Surface Salinity 98-21")
  
# See the data distribution (helps detect outliers)
thetao_values <- as.vector(thetao_mean1998_2021_slopes)
summary(thetao_values)
bottomT_values <- as.vector(bottomT_mean1998_2021_slopes)
summary(bottomT_values)
so_values <- as.vector(so_mean1998_2021_slopes)
summary(so_values)

# Plot the histogram of temperature
par(mfrow=c(2,2))
hist(thetao_values, 
     breaks = 50,  # Number of bins
     main = "Distribution of Matrix Values so", 
     xlab = "Values", 
     col = "blue", 
     border = "black")

p10_thetao <-quantile(thetao_values, probs = 0.02, na.rm = TRUE)
p90_thetao <- quantile(thetao_values, probs = 0.98, na.rm = TRUE)
# Plot ablines
abline(v = p10_thetao, col = "red", lwd = 2, lty = 2)
abline(v = p90_thetao, col = "green", lwd = 
         2, lty = 2)

# Plot the histogram of bottom Temperature
hist(bottomT_values, 
     breaks = 50,  # Number of bins
     main = "Distribution of Matrix Values so", 
     xlab = "Values", 
     col = "blue", 
     border = "black")

p10_bottomT <-quantile(bottomT_values, probs = 0.02, na.rm = TRUE)
p90_bottomT <- quantile(bottomT_values, probs = 0.98, na.rm = TRUE)
# Plot ablines
abline(v = p10_bottomT, col = "red", lwd = 2, lty = 2)
abline(v = p90_bottomT, col = "green", lwd = 
         2, lty = 2)

# Plot the histogram of salinity
hist(so_values, 
     breaks = 50,  # Number of bins
     main = "Distribution of Matrix Values so", 
     xlab = "Values", 
     col = "blue", 
     border = "black")

p10_so <-quantile(so_values, probs = 0.02, na.rm = TRUE)
p90_so <- quantile(so_values, probs = 0.98, na.rm = TRUE)
# Plot ablines
abline(v = p10_so, col = "red", lwd = 2, lty = 2)
abline(v = p90_so, col = "green", lwd = 
         2, lty = 2)


#-----------------------------------#
### 3.1.3. Transform to raster ----
#-----------------------------------#
 
# Create the bounding box of the area of my interest 
# Must have the same dimensions (lat, lon) as the matrix in question
bb_phy<- extent(min(lon_range_phy),max(lon_range_phy),min(lat_range_phy),max(lat_range_phy))  

#Transform each variable to raster
# Temperature
thetao_mean1998_2021_slopes_raster <- raster(thetao_mean1998_2021_slopes) # transform the matrix into a raster
thetao_mean1998_2021_slopes_raster <- t(flip(thetao_mean1998_2021_slopes_raster, direction="x")) #to rotate the raster
projection(thetao_mean1998_2021_slopes_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
thetao_mean1998_2021_slopes_raster<-setExtent(thetao_mean1998_2021_slopes_raster,bb_phy,keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
thetao_mean1998_2021_slopes_raster

# Bottom temperature
bottomT_mean1998_2021_slopes_raster <- raster(bottomT_mean1998_2021_slopes) # transform the matrix into a raster
bottomT_mean1998_2021_slopes_raster <- t(flip(bottomT_mean1998_2021_slopes_raster, direction="x")) #to rotate the raster
projection(bottomT_mean1998_2021_slopes_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
bottomT_mean1998_2021_slopes_raster<-setExtent(bottomT_mean1998_2021_slopes_raster,bb_phy,keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
bottomT_mean1998_2021_slopes_raster

# Salinity	
so_mean1998_2021_slopes_raster <- raster(so_mean1998_2021_slopes) # transform the matrix into a raster
so_mean1998_2021_slopes_raster <- t(flip(so_mean1998_2021_slopes_raster, direction="x")) #to rotate the raster
projection(so_mean1998_2021_slopes_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
so_mean1998_2021_slopes_raster<-setExtent(so_mean1998_2021_slopes_raster,bb_phy,keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
so_mean1998_2021_slopes_raster


#---------------------------------------#
### 3.1.4. Visualization of rasters ----
#---------------------------------------#
#x11();
par(mfrow=c(2,2))
plot(thetao_mean1998_2021_slopes_raster, col = magenta_to_blue_palette(100), main="Raster Surface Temp 98-21")
plot(bottomT_mean1998_2021_slopes_raster, col = magenta_to_blue_palette(100), main="Raster Bottom Temp 98-21")
plot(so_mean1998_2021_slopes_raster, col = magenta_to_blue_palette(100), main="Raster Surface Salinity 98-21")

# Stack
stacked_raster_physical <- stack(thetao_mean1998_2021_slopes_raster,
                                 bottomT_mean1998_2021_slopes_raster, 
                                 so_mean1998_2021_slopes_raster)

# Rename
names(stacked_raster_physical) = c("thetao", "bottomT", "so")

# Plot
par(cex.axis = 1.5)  # Increase the axis text size
plot(stacked_raster_physical, 
     col = magenta_to_blue_palette(100),  
     axes = TRUE)


# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/1.a.Physical_trend.tiff", width = 800, height = 600, res = 300, units= "px")
#dev.off()

#---------------------------------------#
### 3.1.5. Save the rasters ----
#---------------------------------------#
#writeRaster(thetao_mean1998_2021_slopes_raster, "./Outputs/CI_analysis/CI_monthly_anomalies/Rasters/surface_thetao_mean1998_2021_slopes_raster.tif", sep="", overwrite=TRUE) #geotif
#writeRaster(bottomT_mean1998_2021_slopes_raster, "./Outputs/CI_analysis/CI_monthly_anomalies/Rasters/bottomT_mean1998_2021_slopes_raster.tif", sep="", overwrite=TRUE) #geotif
#writeRaster(so_mean1998_2021_slopes_raster, "./Outputs/CI_analysis/CI_monthly_anomalies/Rasters/surface_so_mean1998_2021_slopes_raster.tif", sep="", overwrite=TRUE) #geotif

#---------------------------------------#
## 3.2. SLOPE  biogeochemical data ----
#---------------------------------------#

# Define the function to calculate the annual rate of change
f.rata <- function(x) {
  if (sum(!is.na(x)) < 100) {
    rata <- NA  # If less than 100 valid values, don't compute
  } else {
    ind <- 1:length(x)
    mdl <- lm(x ~ ind, na.action = na.exclude)  # Linear model
    slp <- mdl$coefficients[2]  # Get slope
    intr <- mdl$coefficients[1]  # Get intercept
    
    # Compute rate as evolution per YEAR
    first.intercept <- intr + slp  # Intercept at the first index
    last.ind <- ind[length(ind)]  # Last index
    fin.intercept <- intr + slp * last.ind  # Intercept at the last index
    rata <- (fin.intercept - first.intercept) / last.ind * 12  # Annual rate
  }
  return(rata)
}

#---------------------------------------#
### 3.2.1. Apply the function  ----
#---------------------------------------#
# Applying fun2 (slope calculation) 
# In apply(): c(1,2) indicates that the calculation is performed for both rows and columns
chl_mean1998_2021_slopes <-apply(chl_anom,c(1,2),"f.rata")
nppv_mean1998_2021_slopes <-apply(nppv_anom,c(1,2),"f.rata")
o2_mean1998_2021_slopes <-apply(o2_anom,c(1,2),"f.rata")

#---------------------------------------#
### 3.2.2. Visualization of arrays ----
#---------------------------------------#
# You can visualize arrays with image.plot()
image.plot(chl_mean1998_2021_slopes,col = magenta_to_blue_palette(100), main=" Surface Chl 1998-2021")
image.plot(nppv_mean1998_2021_slopes, col = magenta_to_blue_palette(100), main="Surface NPP 1998-2021")
image.plot(o2_mean1998_2021_slopes, col = magenta_to_blue_palette(100), main="Surface Oxygen 1998-2021")

# See the data distribution (helps detect outliers)
chl_values <- as.vector(chl_mean1998_2021_slopes)
summary(chl_values)

npp_values <- as.vector(nppv_mean1998_2021_slopes)
summary(npp_values)

o2_values <- as.vector(o2_mean1998_2021_slopes)
summary(o2_values)

par(mfrow=c(2,2))
# Plot the histogram of Chl
hist(chl_values, 
     breaks = 50,  # Number of bins
     main = "Distribution of Matrix Values NPP", 
     xlab = "Values", 
     col = "blue", 
     border = "black")

p10_chl <-quantile(chl_values, probs = 0.02, na.rm = TRUE)
p90_chl <- quantile(chl_values, probs = 0.98, na.rm = TRUE)
# Plot histogram
abline(v = p10_chl, col = "red", lwd = 2, lty = 2)
abline(v = p90_chl, col = "green", lwd = 
         2, lty = 2)

# Plot the histogram of NPP
hist(npp_values, 
     breaks = 50,  # Number of bins
     main = "Distribution of Matrix Values NPP", 
     xlab = "Values", 
     col = "blue", 
     border = "black")

p10_npp <-quantile(npp_values, probs = 0.02, na.rm = TRUE)
p90_npp <- quantile(npp_values, probs = 0.98, na.rm = TRUE)
# Plot histogram
abline(v = p10_npp, col = "red", lwd = 2, lty = 2)
abline(v = p90_npp, col = "green", lwd = 
         2, lty = 2)

# Plot the histogram of Oxygen
hist(o2_values, 
     breaks = 50,  # Number of bins
     main = "Distribution of Matrix Values NPP", 
     xlab = "Values", 
     col = "blue", 
     border = "black")

p10_o2 <-quantile(o2_values, probs = 0.02, na.rm = TRUE)
p90_o2 <- quantile(o2_values, probs = 0.98, na.rm = TRUE)
# Plot histogram
abline(v = p10_o2, col = "red", lwd = 2, lty = 2)
abline(v = p90_o2, col = "green", lwd = 
         2, lty = 2)

#---------------------------------------#
### 3.2.3. Transform to raster ----
#---------------------------------------#
# Create the bounding box of the area of my interest
bb_bio<-extent(min(lon_range_bio),max(lon_range_bio),min(lat_range_bio),max(lat_range_bio))

# Transform each variable to raster
#Chlorophyll
chl_mean1998_2021_slopes_raster <- raster(chl_mean1998_2021_slopes) # transform the matrix into a raster
chl_mean1998_2021_slopes_raster <- t(flip(chl_mean1998_2021_slopes_raster, direction="x")) #to rotate the raster
projection(chl_mean1998_2021_slopes_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
chl_mean1998_2021_slopes_raster<-setExtent(chl_mean1998_2021_slopes_raster,bb_bio,keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
chl_mean1998_2021_slopes_raster

#NPP
nppv_mean1998_2021_slopes_raster <- raster(nppv_mean1998_2021_slopes) # transform the matrix into a raster
nppv_mean1998_2021_slopes_raster <- t(flip(nppv_mean1998_2021_slopes_raster, direction="x")) #to rotate the raster
projection(nppv_mean1998_2021_slopes_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
nppv_mean1998_2021_slopes_raster<-setExtent(nppv_mean1998_2021_slopes_raster,bb_bio,keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
nppv_mean1998_2021_slopes_raster

#Oxygen
o2_mean1998_2021_slopes_raster <- raster(o2_mean1998_2021_slopes) # transform the matrix into a raster
o2_mean1998_2021_slopes_raster <- t(flip(o2_mean1998_2021_slopes_raster, direction="x")) #to rotate the raster
projection(o2_mean1998_2021_slopes_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
o2_mean1998_2021_slopes_raster<-setExtent(o2_mean1998_2021_slopes_raster,bb_bio,keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
o2_mean1998_2021_slopes_raster

#---------------------------------------#
### 3.2.4. Visualization of rasters ----
#---------------------------------------#
# You can visualize rasters with plot()
plot(chl_mean1998_2021_slopes_raster, col = blue_to_red_palette(100), main="Raster Surfae Chl 1998-2021")
plot(nppv_mean1998_2021_slopes_raster, col = blue_to_red_palette(100), main="Raster Surface NPP 1998-2021")
plot(o2_mean1998_2021_slopes_raster, col = blue_to_red_palette(100), main="Raster Surface Oxygen 1998-2021")

#Stack
stacked_raster_biogeo <- stack(chl_mean1998_2021_slopes_raster, 
  	                         nppv_mean1998_2021_slopes_raster,
                           o2_mean1998_2021_slopes_raster)

#Rename
names(stacked_raster_biogeo) = c("chl", "nppv", "o2")

# Plot
par(cex.axis = 1.5)  # Increase the axis text size
plot(stacked_raster_biogeo, 
     col = magenta_to_blue_palette(100),  
     axes = TRUE)

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/1.b.BGC_trend.tiff",width = 800, height = 600, res = 300, units= "px")
#dev.off()

#---------------------------------------#	
### 3.2.5. Save the rasters ----
#---------------------------------------#
#writeRaster(chl_mean1998_2021_slopes_raster, "./Outputs/CI_analysis/CI_monthly_anomalies/Rasters/surface_chl_mean1998_2021_slopes_raster.tif", sep="", overwrite=TRUE) #geotif
#writeRaster(nppv_mean1998_2021_slopes_raster, "./Outputs/CI_analysis/CI_monthly_anomalies/Rasters/surface_nppv_mean1998_2021_slopes_raster.tif", sep="", overwrite=TRUE) #geotif
#writeRaster(o2_mean1998_2021_slopes_raster,  "./Outputs/CI_analysis/CI_monthly_anomalies/Rasters/surface_o2_mean1998_2022_slopes_raster.tif", sep="", overwrite=TRUE) #geotif

#------------------------------#
# 4. CUMULATIVE IMPACTS    ---- 
#------------------------------#

#---------------------------------------#
## 4.1. CUMULATIVE IMPACTS Product 1 ----
#---------------------------------------#

#---------------------------------------#
### 4.1.1. Calculate absolute values ----
#---------------------------------------#
# We consider as chang2 all the slope values (positive or negative) 
abs_thetao_mean1998_2021_slopes_raster <- abs(thetao_mean1998_2021_slopes_raster)
abs_bottomT_mean1998_2021_slopes_raster <- abs(bottomT_mean1998_2021_slopes_raster)
abs_so_mean1998_2021_slopes_raster <- abs(so_mean1998_2021_slopes_raster)

#----------------------------------------------------#
### 4.1.2. Normalization Physical Impacts ----
#-------------------------------------------------------#
# We use the 10th and 90th percentiles based normalization. To reduce the impact of extreme values

# Temperature
par(mfrow=c(2,2))
# Calculate quantiles
p10_value_thetao <-quantile(values(abs_thetao_mean1998_2021_slopes_raster), probs = 0.02, na.rm = TRUE)
p90_value_thetao <- quantile(values(abs_thetao_mean1998_2021_slopes_raster), probs = 0.98, na.rm = TRUE)
# Plot histogram
hist(values(abs_thetao_mean1998_2021_slopes_raster), 
     breaks = 50, main = "Histogram of Raster Values",
     xlab = "Values", col = "lightblue", border = "black")
# Ablines
abline(v = p10_value_thetao, col = "red", lwd = 2, lty = 2)
abline(v = p90_value_thetao, col = "green", lwd = 
         2, lty = 2)
# Normalize using 10th and 90th percentiles (percentile based normalization)
normalize_thetao_mean1998_2021_slopes_raster <- (abs_thetao_mean1998_2021_slopes_raster - p10_value_thetao) / (p90_value_thetao - p10_value_thetao)
plot(normalize_thetao_mean1998_2021_slopes_raster)
# Duplicate to avoid overwrite
norm_thetao_mean1998_2021_slopes_raster<-normalize_thetao_mean1998_2021_slopes_raster
# Correct to 0-1
norm_thetao_mean1998_2021_slopes_raster[norm_thetao_mean1998_2021_slopes_raster < 0] <- 0
norm_thetao_mean1998_2021_slopes_raster[norm_thetao_mean1998_2021_slopes_raster > 1] <- 1
# Plot
plot(norm_thetao_mean1998_2021_slopes_raster, col= magenta_to_blue_palette(100), main="Temperature Trend 1998-2021")

# Bottom temperature
par(mfrow=c(2,2))
# Calculate quantiles
p10_value_bottomT <-quantile(values(abs_bottomT_mean1998_2021_slopes_raster), probs = 0.02, na.rm = TRUE)
p90_value_bottomT <- quantile(values(abs_bottomT_mean1998_2021_slopes_raster), probs = 0.98, na.rm = TRUE)
# Plot histogram
hist(values(abs_bottomT_mean1998_2021_slopes_raster), 
     breaks = 50, main = "Histogram of Raster Values Bottom T",
     xlab = "Values", col = "lightblue", border = "black")
# Ablines
abline(v = p10_value_bottomT, col = "red", lwd = 2, lty = 2)
abline(v = p90_value_bottomT, col = "green", lwd = 
         2, lty = 2)
# Normalize using 10th and 90th percentiles (percentile based normalization)
normalize_bottomT_mean1998_2021_slopes_raster <- (abs_bottomT_mean1998_2021_slopes_raster - p10_value_bottomT) / (p90_value_bottomT - p10_value_bottomT)
plot(normalize_bottomT_mean1998_2021_slopes_raster)
# Duplicate to avoid overwrite
norm_bottomT_mean1998_2021_slopes_raster<-normalize_bottomT_mean1998_2021_slopes_raster
# Correct to 0-1
norm_bottomT_mean1998_2021_slopes_raster[norm_bottomT_mean1998_2021_slopes_raster < 0] <- 0
norm_bottomT_mean1998_2021_slopes_raster[norm_bottomT_mean1998_2021_slopes_raster > 1] <- 1
# Plot
plot(norm_bottomT_mean1998_2021_slopes_raster, col= magenta_to_blue_palette(100), main="Bottom Temperature Trend 1998-2021")

# Salinity
par(mfrow=c(2,2))
# Calculate quantiles
p10_value_so <-quantile(values(abs_so_mean1998_2021_slopes_raster), probs = 0.02, na.rm = TRUE)
p90_value_so <- quantile(values(abs_so_mean1998_2021_slopes_raster), probs = 0.98, na.rm = TRUE)
# Plot histogram
hist(values(abs_so_mean1998_2021_slopes_raster), 
     breaks = 50, main = "Histogram of Raster Values Salinity",
     xlab = "Values", col = "lightblue", border = "black")
# Ablines
abline(v = p10_value_so, col = "red", lwd = 2, lty = 2)
abline(v = p90_value_so, col = "green", lwd = 
         2, lty = 2)
# Normalize using 10th and 90th percentiles (percentile based normalization)
normalize_so_mean1998_2021_slopes_raster <- (abs_so_mean1998_2021_slopes_raster - p10_value_so) / (p90_value_so - p10_value_so)
plot(normalize_so_mean1998_2021_slopes_raster)
# Duplicate to avoid overwrite
norm_so_mean1998_2021_slopes_raster<-normalize_so_mean1998_2021_slopes_raster
# Correct to 0-1
norm_so_mean1998_2021_slopes_raster[norm_so_mean1998_2021_slopes_raster < 0] <- 0
norm_so_mean1998_2021_slopes_raster[norm_so_mean1998_2021_slopes_raster > 1] <- 1
# Plot
plot(norm_so_mean1998_2021_slopes_raster, col= magenta_to_blue_palette(100), main="Salinity Trend 1998-2021")

#------------------------------------------------------------------#
### 4.1.3. Create a collection of physical rasters (normalized) ----
#-------------------------------------------------------------------#

# Stack
stacked_raster <- stack(norm_thetao_mean1998_2021_slopes_raster,
                        norm_bottomT_mean1998_2021_slopes_raster, 
                        norm_so_mean1998_2021_slopes_raster)

# Rename
names(stacked_raster) = c("thetao","bottomT", "so" )

# Plot
par(cex.axis = 1.5)
plot(stacked_raster, 
     col = magenta_to_blue_palette(100),  
     axes = TRUE)

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/2.a.Physical_trend_normalizes.tiff", width = 800, height = 600, res = 300, units= "px")
#dev.off()

#---------------------------------------#
## 4.2. CUMULATIVE IMPACTS Product 2 ----
#---------------------------------------#   

#---------------------------------------#
### 4.2.1. Calculate absolute values ----
#---------------------------------------#
# We consider as "change" all the slope values (positive or negative) 
abs_chl_mean1998_2021_slopes_raster <- abs(chl_mean1998_2021_slopes_raster)
abs_nppv_mean1998_2021_slopes_raster <- abs(nppv_mean1998_2021_slopes_raster)
abs_o2_mean1998_2021_slopes_raster <- abs(o2_mean1998_2021_slopes_raster)

#--------------------------------------------------#
### 4.2.2. Normalization Biogeochemical Impacts ----
#--------------------------------------------------#
#Chl
par(mfrow=c(2,2))
# Calculate quantiles
p10_value_chl <-quantile(values(abs_chl_mean1998_2021_slopes_raster), probs = 0.02, na.rm = TRUE)
p90_value_chl <- quantile(values(abs_chl_mean1998_2021_slopes_raster), probs = 0.98, na.rm = TRUE)
# Plot histogram
hist(values(abs_chl_mean1998_2021_slopes_raster), 
     breaks = 50, main = "Histogram of Raster Values Chl",
     xlab = "Values", col = "lightblue", border = "black")
# Ablines
abline(v = p10_value_chl, col = "red", lwd = 2, lty = 2)
abline(v = p90_value_chl, col = "green", lwd = 
         2, lty = 2)
# Normalize using 10th and 90th percentiles (percentile based normalization)
normalize_chl_mean1998_2021_slopes_raster <- (abs_chl_mean1998_2021_slopes_raster - p10_value_chl) / (p90_value_chl - p10_value_chl)
plot(normalize_chl_mean1998_2021_slopes_raster)
# Duplicate to avoid overwrite
norm_chl_mean1998_2021_slopes_raster<-normalize_chl_mean1998_2021_slopes_raster
# Correct to 0-1
norm_chl_mean1998_2021_slopes_raster[norm_chl_mean1998_2021_slopes_raster < 0] <- 0
norm_chl_mean1998_2021_slopes_raster[norm_chl_mean1998_2021_slopes_raster > 1] <- 1
# Plot
plot(norm_chl_mean1998_2021_slopes_raster, col= magenta_to_blue_palette(100), main="Chlorophyll Trend 1998-2021")

#NPP
par(mfrow=c(2,2))
# Calculate quantiles
p10_value_npp <-quantile(values(abs_nppv_mean1998_2021_slopes_raster), probs = 0.02, na.rm = TRUE)
p90_value_npp <- quantile(values(abs_nppv_mean1998_2021_slopes_raster), probs = 0.98, na.rm = TRUE)
# Plot histogram
hist(values(abs_nppv_mean1998_2021_slopes_raster), 
     breaks = 50, main = "Histogram of Raster Values NPP",
     xlab = "Values", col = "lightblue", border = "black")
# Ablines
abline(v = p10_value_npp, col = "red", lwd = 2, lty = 2)
abline(v = p90_value_npp, col = "green", lwd = 
         2, lty = 2)
# Normalize using 10th and 90th percentiles (percentile based normalization)
normalize_nppv_mean1998_2021_slopes_raster <- (abs_nppv_mean1998_2021_slopes_raster - p10_value_npp) / (p90_value_npp - p10_value_npp)
plot(normalize_nppv_mean1998_2021_slopes_raster)
# Duplicate to avoid overwrite
norm_nppv_mean1998_2021_slopes_raster<-normalize_nppv_mean1998_2021_slopes_raster
# Correct to 0-1
norm_nppv_mean1998_2021_slopes_raster[norm_nppv_mean1998_2021_slopes_raster < 0] <- 0
norm_nppv_mean1998_2021_slopes_raster[norm_nppv_mean1998_2021_slopes_raster > 1] <- 1
# Plot
plot(norm_nppv_mean1998_2021_slopes_raster, col= magenta_to_blue_palette(100), main="NPP Trend 1998-2021")

#Oxygen
par(mfrow=c(2,2))
# Calculate quantiles
p10_value_o2 <-quantile(values(abs_o2_mean1998_2021_slopes_raster), probs = 0.02, na.rm = TRUE)
p90_value_o2 <- quantile(values(abs_o2_mean1998_2021_slopes_raster), probs = 0.98, na.rm = TRUE)
# Plot histogram
hist(values(abs_o2_mean1998_2021_slopes_raster), 
     breaks = 50, main = "Histogram of Raster Values Oxygen",
     xlab = "Values", col = "lightblue", border = "black")
# Ablines
abline(v = p10_value_o2, col = "red", lwd = 2, lty = 2)
abline(v = p90_value_o2, col = "green", lwd = 
         2, lty = 2)
# Normalize using 10th and 90th percentiles (percentile based normalization)
normalize_o2_mean1998_2021_slopes_raster <- (abs_o2_mean1998_2021_slopes_raster - p10_value_o2) / (p90_value_o2 - p10_value_o2)
plot(normalize_o2_mean1998_2021_slopes_raster)
# Duplicate to avoid overwrite
norm_o2_mean1998_2021_slopes_raster<-normalize_o2_mean1998_2021_slopes_raster
# Correct to 0-1
norm_o2_mean1998_2021_slopes_raster[norm_o2_mean1998_2021_slopes_raster < 0] <- 0
norm_o2_mean1998_2021_slopes_raster[norm_o2_mean1998_2021_slopes_raster > 1] <- 1
# Plot
plot(norm_o2_mean1998_2021_slopes_raster,col= magenta_to_blue_palette(100), main="Oxygen Trend 1998-2021")

#------------------------------------------------------------------------#
### 4.2.3. Create a collection of biogeochemical rasters (normalized) ----
#--------------------------------------------------------------------------#
# Stack
 stacked_raster2 <- stack(norm_chl_mean1998_2021_slopes_raster, 
                          norm_nppv_mean1998_2021_slopes_raster,
                          norm_o2_mean1998_2021_slopes_raster)

# Rename
names(stacked_raster2) = c("chl", "nppv", "o2")

# Plot
par(cex.axis = 1.5)
plot(stacked_raster2, 
     col = magenta_to_blue_palette(100),  
     axes = TRUE)

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/2.b.BGC_trend_normalizes.tiff",  width = 800, height = 600, res = 300, units= "px")
#dev.off()


# ----------------------------------------------------- #
## 4.3. CUMULATIVE IMPACTS Both ----
# ----------------------------------------------------- #

# ----------------------------------------------------- #
### 4.3.1. check resolution and resample if needed ----
# ----------------------------------------------------- #

# Create new one to avoid overwriting them
stand_BGC_impact<-stacked_raster2 
stand_phy_impact<-stacked_raster

# Check both resolution to see if resample is needed
# Always resample to the coarser resolution
stand_BGC_impact #(x,y)
stand_phy_impact #(x,y)

# Resample
stand_BGC_impact_res <- resample(stand_BGC_impact, stand_phy_impact, method="bilinear", na.rm=TRUE)

# Checkresampled object
stand_BGC_impact_res
summary(stand_BGC_impact_res)
plot(stand_BGC_impact_res,col=blue_to_red_palette(100), main= "All BGC impacts" )

# Ensure the min and max still [0-1]. Resampling can sometimes alter this and will affect the following steps.
# This function changes the min and max values
standardize_layer <- function(layer) {
  min_value <- cellStats(layer, min, na.rm = TRUE)
  max_value <- cellStats(layer, max, na.rm = TRUE)
  scaled_layer <- (layer - min_value) / (max_value - min_value)
  return(scaled_layer)
}

# Apply standardize_layer to each layer of the RasterBrick
stand_BGC_impact_final <- stack(lapply(1:nlayers(stand_BGC_impact_res), function(i) {
  standardize_layer(raster(stand_BGC_impact_res, layer=i))
}))

# Print the standardized values
print(stand_BGC_impact_final)
plot(stand_BGC_impact_final, col=blue_to_red_palette(100), main= "All bgc impacts rescaled" )


# ----------------------------------------------------- #
### 4.3.2. create a unique layer of Physical impacts----
# ----------------------------------------------------- #

# Add all the layers of one raster stack
CI_phy <- sum(stand_phy_impact) # sum of all the RasterStacks

# Plot
par(mfrow=c(2,2))
plot(CI_phy, col=magenta_to_blue_palette(100), main= "Physical cumulative impacts")

# Save raster
#writeRaster(CI_phy, "./Outputs/CI_analysis/CI_monthly_anomalies/Rasters/surface_CI_sum_physical_1998-2021_norm2_98.tif", overwrite=TRUE) #geotif

# ----------------------------------------------------- #
### 4.3.3. create a unique layer of BGC impacts----
# ----------------------------------------------------- #
# Add all the layers of one raster stack
CI_BGC <- sum(stand_BGC_impact_final) # sum of all the RasterStacks

# Plot
plot(CI_BGC, col=magenta_to_blue_palette(100), main= "Biogeochemical Cumulative Impacts")

#Save Raster
#writeRaster(CI_BGC, "./Outputs/CI_analysis/CI_monthly_anomalies/Rasters/surface_CI_sum_BGC_1998-2021_norm_2_98.tif", overwrite=TRUE) #geotif

# If we want to save the image (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/3.Physical_sum_&_BGC_sum.tiff", width = 800, height = 600, res = 300, units= "px")
#dev.off()


# ---------------------------------------------------------------- #
### 4.3.4. create a unique layer of Physical and BGC impacts----
# ---------------------------------------------------------- #

# Add the two stacks of variables phy and bio together
CI_sum <- CI_phy + CI_BGC

# Plot
par(mfrow=c(1,1))
plot(CI_sum,col=magenta_to_blue_palette(100), main= "All Cumulative Impacts" )

#Save raster
#writeRaster(CI_sum,  "./Outputs/CI_analysis/CI_monthly_anomalies/Rasters/surface_CI_sum_1998-2021_norm_2_98.tif", sep="", overwrite=TRUE) #geotif

# If we want to save the image (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/4.All_sum..tiff", width = 800, height = 600, res = 300, units= "px")
#dev.off()


# ----------------------------------------------------- #
## 4.4. RECLASSIFY RASTERS INTO ITS QUARTILES ----
# ----------------------------------------------------- #

# ------------------------------------------------------------ #
### 4.4.1. Reclassify Product 1 into its quantiles ----
# ------------------------------------------------------------ #  

# Calculate the quantiles for the physical variables
quantiles1 <- quantile(values(CI_phy), probs = c(0.25, 0.5, 0.75), na.rm=TRUE)
  
# Initialize a new map to avoid overwriting
phy_qtile_hab <-  CI_phy 

# Data frame with the min, Q1, Q2, Q3 and max. 
rcl_phy <- c(0, quantiles1[1], 1, quantiles1[1],
             quantiles1[2], 2, quantiles1[2], 
             quantiles1[3], 3, quantiles1[3], 
             cellStats(CI_phy, stat=max), 4)
# Transform into a matrix
rcl_phy <- matrix(rcl_phy, ncol=3, byrow=TRUE) #reclass matrix
rcl_phy

# Reclassification of the original raster (CI_phy)into a 1 to 4 numerations
phy_qtile_hab <-reclassify(CI_phy, rcl_phy, na.rm=TRUE) 

# Plot
par(mfrow=c(2,2))
plot(phy_qtile_hab, col=magenta_to_blue_palette(100), main="Physical quartiles")

  
# ----------------------------------------------------- #
### 4.4.2. Reclassify Product 2 into its quartiles ----
# ----------------------------------------------------- #
# Calculate the quantiles for the physical variables
quantiles2 <- quantile(values(CI_BGC), probs = c(0.25, 0.5, 0.75), na.rm=TRUE)
    
# Initialize a new map 
bgc_qtile_hab <- CI_BGC
    
# Data frame with the min, Q1, Q2, Q3 and max. 
rcl_BGC <- c(0, quantiles2[1], 1, quantiles2[1], 
             quantiles2[2], 2, quantiles2[2], 
             quantiles2[3], 3, quantiles2[3], 
             cellStats(CI_BGC, stat=max), 4)
# Transform into a matrix
rcl_BGC <- matrix(rcl_BGC, ncol=3, byrow=TRUE) #reclass matrix
rcl_BGC

# Reclassification of the original raster (CI_phy)into a 1 to 4 numerations
bgc_qtile_hab <- reclassify(CI_BGC, rcl_BGC, na.rm=TRUE) 

# Plot
plot(bgc_qtile_hab, col=magenta_to_blue_palette(100), main="BGC quartiles")

# If we want to save the image (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/5.Physical_&_BGC_quantiles.tiff", width = 800, height = 600, res = 300, units= "px")
#dev.off()


# ----------------------------------------------------- #
## 4.5. FINAL CUMULATIVE IMPACT MAP----
# ----------------------------------------------------- #

# ----------------------------------------- #
### 4.5.1. Visualization bivariate plot ----
# ------------------------------------------ #
    
# Palettes
p <- c('#9e0142', '#d53e4f', '#f46d43', '#fdae61')
    
c3 <- c('#3288bd','#32A8BD' , '#7EC9D6', '#93DACC')

# Create a dataframe
colors_df <- crossing(x = 0:4, y = 0:4) %>%
      mutate(bi_class = paste(x, y, sep = '-')) %>%
      mutate(bi_fill = case_when(bi_class == '0-0' ~ '#ffffff',
                                 bi_class %in% c('0-1') ~ colorspace::lighten(p[4], .8),
                                 bi_class %in% c('0-2') ~ colorspace::lighten(p[4], .6),
                                 bi_class %in% c('0-3') ~ colorspace::lighten(p[4], .2),
                                 bi_class %in% c('1-0') ~ colorspace::lighten(c3[4], .8),
                                 bi_class %in% c('2-0') ~ colorspace::lighten(c3[4], .6),
                                 bi_class %in% c('3-0') ~ colorspace::lighten(c3[4], .2),
    
                                 bi_class %in% c('1-1') ~ '#f0f0f0',
                                # bi_class %in% c(, '2-1') ~ '#d9d9d9',
                                 bi_class %in% c('2-2') ~ '#bdbdbd',
                                 #bi_class %in% c(, '3-2') ~ '#969696',
                                 bi_class %in% c('3-3') ~ '#737373',
    
                                 bi_class %in% c('0-4') ~ p[4],
                                 bi_class %in% c('1-4') ~ p[3],
                                 bi_class %in% c('1-3') ~ colorspace::lighten(p[3], .2),
                                 bi_class %in% c('1-2') ~ colorspace::lighten(p[3], .4),
                                 bi_class %in% c('2-4') ~ p[2],
                                 bi_class %in% c('2-3') ~ colorspace::lighten(p[2], .4),
                                 bi_class %in% c('3-4') ~ p[1],
    
                                 bi_class %in% c('4-0') ~ c3[4],
                                 bi_class %in% c('4-1') ~ c3[3],
                                 bi_class %in% c('3-1') ~ colorspace::lighten(c3[3], .2),
                                 bi_class %in% c('2-1') ~ colorspace::lighten(c3[3], .4),
                                 bi_class %in% c('4-2') ~ c3[2],
                                 bi_class %in% c('3-2') ~ colorspace::lighten(c3[2], .4),
                                 bi_class %in% c('4-3') ~ c3[1],
                                 bi_class == '4-4' ~ 'black', ### black
                                 TRUE ~ 'blue'))

# Create a ggplot object (for the palette legend)
color_plot <- ggplot(colors_df, aes(x = x, y = y, fill = bi_fill, label = bi_class)) +
  geom_tile(color = "black") +  # Add tiles with black borders
  geom_text(size = 4) +  # Add text labels
  scale_fill_identity() +  # Use fill colors as they are
  labs(x = NULL, y = NULL) +  # Remove axis labels
  theme_void()  # Remove background, gridlines, and axis ticks

color_plot

# Create the plot of the final palette legend
panel_lgd_raw <- ggplot(colors_df, aes(x, y, fill = bi_class)) +
  geom_tile(show.legend = FALSE, color = 'white', linewidth = .25) +
  scale_fill_manual(breaks = colors_df$bi_class, values = colors_df$bi_fill) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))+
      theme_void()
panel_lgd_raw

# Convert to data frames (x=long, y=lat)
phy_df <- as.data.frame(phy_qtile_hab, xy = TRUE, na.rm=TRUE) %>% 
  setNames(c('x', 'y', 'phy'))
summary(phy_qtile_hab)
unique(phy_df$phy)

bgc_df <- as.data.frame(bgc_qtile_hab, xy = TRUE, na.rm=TRUE) %>% # he añadido lo de na.rm=TRUE
  setNames(c('x', 'y', 'bgc'))
summary(bgc_qtile_hab)
unique(bgc_df$bgc) 

# Merge PHY and BGC data frames
phy_bgc_df <- merge(phy_df, bgc_df, by = c('x', 'y'), na.rm=TRUE)
phy_bgc_df <- mutate(phy_bgc_df, bi_class = paste(phy, bgc, sep = '-'))
phy_bgc_df

# Plot the data
# with axis
biplot_hab2 <- ggplot() +
  geom_raster(data = phy_bgc_df, aes(x, y, fill = bi_class)) +
  scale_fill_manual(breaks = colors_df$bi_class, values = colors_df$bi_fill, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +  # Remove axis labels, but keep axes visible
  theme_bw() +  # Use theme_bw for visible axis lines
  theme(
    panel.border = element_blank(),  # Remove the border around the plot area
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove panel background
    axis.line = element_line(color = "black"),  # Keep axis lines visible
    axis.title = element_blank(),  # Ensure axis titles are blank
    axis.text = element_text(size = 12)  # Customize axis tick labels
  )

biplot_hab2

# Assemble the figure (map + legend)
panel_a_map <- cowplot::get_panel(biplot_hab2)
panel_a_lgd <- cowplot::get_panel(panel_lgd_raw)
lgd_p <- c(x = 0.75, y = .8, h = .15, w = .15) #right corner
#lgd_p <- c(x = 0.05, y = .85, h = .15, w = .15) #left corner

figS2_panel <- cowplot::ggdraw() +
  cowplot::draw_plot(biplot_hab2, x = 0, y = 0, height = 1, width = 1) +  # Use the entire `biplot_hab2` object
  ### Draw legends:
  cowplot::draw_plot(panel_lgd_raw, x = lgd_p['x'], y = lgd_p['y'], height = lgd_p['h'], width = lgd_p['w']) +
  ### Legend labels:
  cowplot::draw_label('Cum Phy',
                      x = lgd_p['x'], y = lgd_p['y'] - .01,
                      hjust = 0, vjust = 1, angle = 0,
                      size = 12, color = 'black') +
  cowplot::draw_label('Cum BGC',
                      x = lgd_p['x'] - .01, y = lgd_p['y'],
                      hjust = 0, vjust = 0, angle = 90,
                      size = 12, color = 'black')

figS2_panel   

# If we want to save the image (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/6.a.CI_based_on_monthly_anomalies_surface.tiff", width = 800, height = 600, res = 300, units= "px")
#dev.off()

# ----------------------------------------------------- #
### 4.5.2. Rasterize Cumulative impacts map ----
# ----------------------------------------------------- #

# Ensure bi_class is a factor and convert to integer codes
phy_bgc_df$bi_class <- as.factor(phy_bgc_df$bi_class)
#the function raster() cant rasterize a filed in the form of 1-1/1-2/4-4, etc so we need to create a single integer for eahc code
phy_bgc_df$bi_class_int <- as.integer(phy_bgc_df$bi_class)
str(phy_bgc_df)

# Check which integer has been assigned to each number code
unique_combinations <- unique(phy_bgc_df[c("bi_class", "bi_class_int")])
unique_combinations <- unique_combinations[order(unique_combinations$bi_class), ]
print(unique_combinations)

# create an empty raster first with the same extent and resolution as the original rasters
CI_r <- raster(extent(phy_qtile_hab), nrow = nrow(phy_qtile_hab), ncol = ncol(phy_qtile_hab))
crs(CI_r) <- crs(phy_qtile_hab)
CI_r

# Convert the data frame to a SpatialPointsDataFrame
phy_bgc_spatial<-phy_bgc_df
coordinates(phy_bgc_spatial) <- ~x + y

# Rasterize the bi_class_int column from the phy_bgc_df
CI_raster <- rasterize(phy_bgc_spatial, CI_r, field="bi_class_int")

# ----------------------------------------------------- #
#### 4.5.2.1. Visaulization as Raster Cumulative impacts ----
# ----------------------------------------------------- #

# Create a progressive color palette from blue to white to red
# Define the initial palettes
p <- c('#9e0142', '#d53e4f', '#f46d43', '#fdae61')
c3 <- c('#3288bd','#32A8BD' , '#7EC9D6', '#93DACC')

# Create the `colors_df` with the correct order
colors_df <- crossing(x = 0:4, y = 0:4) %>%
  mutate(bi_class = paste(x, y, sep = '-')) %>%
  mutate(bi_fill = case_when(
    bi_class == '0-0' ~ '#ffffff',
    bi_class %in% c('0-1') ~ colorspace::lighten(p[4], .8),
    bi_class %in% c('0-2') ~ colorspace::lighten(p[4], .6),
    bi_class %in% c('0-3') ~ colorspace::lighten(p[4], .2),
    bi_class %in% c('1-0') ~ colorspace::lighten(c3[4], .8),
    bi_class %in% c('2-0') ~ colorspace::lighten(c3[4], .6),
    bi_class %in% c('3-0') ~ colorspace::lighten(c3[4], .2),
    bi_class %in% c('1-1') ~ '#f0f0f0',
    bi_class %in% c('2-2') ~ '#bdbdbd',
    bi_class %in% c('3-3') ~ '#737373',
    bi_class %in% c('0-4') ~ p[4],
    bi_class %in% c('1-4') ~ p[3],
    bi_class %in% c('1-3') ~ colorspace::lighten(p[3], .2),
    bi_class %in% c('1-2') ~ colorspace::lighten(p[3], .4),
    bi_class %in% c('2-4') ~ p[2],
    bi_class %in% c('2-3') ~ colorspace::lighten(p[2], .4),
    bi_class %in% c('3-4') ~ p[1],
    bi_class %in% c('4-0') ~ c3[4],
    bi_class %in% c('4-1') ~ c3[3],
    bi_class %in% c('3-1') ~ colorspace::lighten(c3[3], .2),
    bi_class %in% c('2-1') ~ colorspace::lighten(c3[3], .4),
    bi_class %in% c('4-2') ~ c3[2],
    bi_class %in% c('3-2') ~ colorspace::lighten(c3[2], .4),
    bi_class %in% c('4-3') ~ c3[1],
    bi_class == '4-4' ~ 'black', 
    TRUE ~ 'blue'
  ))

# Extract the colors in the exact order of bi_class
color_palette <- setNames(colors_df$bi_fill, colors_df$bi_class)

# 1. Create the base plot without axis labels
par(mfrow=c(1,1))
plot(CI_raster, 
     main = "Cumulative impact Black Sea waters (surface)", 
     col = color_palette[levels(factor(phy_bgc_df$bi_class))], 
     asp = NA,
     xaxt = 'n',  # Suppress default x-axis labels
     yaxt = 'n')  # Suppress default y-axis labels

# Add custom x and y axis labels
axis(1, col.axis = "black", cex.axis = 1.2)  # Custom x-axis labels
axis(2, col.axis = "black", cex.axis = 1.2, las = 1)  


# 2. Add the shapefile (`areas`) on top of the raster
plot(World_transformed,
     col = "lightgrey",    # Use `NA` to keep the polygon transparent
     border = "black",  # Define border color for better visibility
     add = TRUE)  # Overlay on top of the raster

# # 3. Add EEZ boundaries (optional)
# plot(EEZ_transformed,
#      col = adjustcolor("red", alpha.f = 0.5), # Keep the polygon fill transparent
#      border = NA,   # Use a distinct color for land boundaries
#      lwd = 2,     # Increase line width for emphasis
#      add = TRUE)  # Overlay on the existing plot

# 4. Add a matching legend
legend("topright", 
       legend = levels(factor(phy_bgc_df$bi_class)), 
       fill = color_palette[levels(factor(phy_bgc_df$bi_class))], 
       title = "Bi Class", 
       cex = 0.8)

# 5. Draw a black box around the plot
box(col = "black", lwd = 1)  # Adjust the line width `lwd` as necessary

# If we want to save the image (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/6.b.CI_based_on_monthly_anomalies_surface_countries.tiff", width = 800, height = 600, res = 300, units= "px")
#dev.off()

# ----------------------------------------------------- #
#### 4.5.2.2. Save final CI raster ----
# ----------------------------------------------------- #

#save raster
#writeRaster(CI_raster, "./Outputs/CI_analysis/CI_monthly_anomalies/CI_based_on_monthly_anomalies_surface_raster.tif", overwrite=TRUE) #geotif

# ----------------------------------------------------- #
### 4.5.3.Transform to Polygon----
# ----------------------------------------------------- #

# Convert the raster to a SpatRaster object
CI_terra <- rast(CI_raster)

# Convert the raster to polygons
CI_polygons <- as.polygons(CI_terra)

# ----------------------------------------------------- #
#### 4.5.3.1. Plot as Polygon----
# ----------------------------------------------------- #

# Plot the polygons using the layer attribute for colors
plot(CI_polygons, col =color_palette[levels(factor(phy_bgc_df$bi_class))], main = "Polygons representing raster levels", )
spplot(CI_polygons, col.regions = color_palette[levels(factor(phy_bgc_df$bi_class))], main = "Polygons representing raster levels")

# ----------------------------------------------------- #
#### 4.5.3.2. Save as Polygon----
# ----------------------------------------------------- #

# Write the polygons to a shapefile
# Convert the polygons to an sf object first
sf_polygons <- st_as_sf(CI_polygons)

# Write the shapefile
#st_write(sf_polygons, "./Outputs/CI_analysis/CI_monthly_anomalies/Polygon/CI_based_on_monthly_anomalies_surface.shp")

#------------------------------#
# 5. CALCULATE THE P VALUE  ---- 
#------------------------------#	

# Create the function to extract the overall ANOVA p-value out of a linear model object

#Alternative p-value function
# lmp <- function (modelobject) 	{
# 								if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
# 								f <- summary(modelobject)$fstatistic
# 								p <- pf(as.numeric(f[1]),as.numeric(f[2]),as.numeric(f[3]),lower.tail=F)
# 								attributes(p) <- NULL
# 								return(p)
# 								}
# 
# fun1 <- function(x) 			{ 
# 								if (is.na(x[1])) {
# 								NA 
# 								} else 	{ 
# 										m = lm(x ~ times_phy) 
# 										lmp(m)
# 										} 
# 								}

#
f.p.val <- function(x){
  if(sum(!is.na(x)) < 100) {
    p.val <- NA
  }
  if(sum(!is.na(x)) >= 100) {
    ind <- 1:length(x)
    mdl <- lm(x ~ ind, na.action=na.exclude)
    # get p-value
    p.val <- summary(mdl)$coefficients[2,4]
  }
  return(p.val)
}

# ----------------------------------------------------- #
## 5.1. PVALUE Product 1 ----
# ----------------------------------------------------- #

# ----------------------------------------------------- #
### 5.1.1. Apply the function ----
# ----------------------------------------------------- #
thetao_mean1998_2021_pvalue <-apply(thetao_anom,c(1,2),"f.p.val")
bottomT_mean1998_2021_pvalue <-apply(bottomT_anom,c(1,2),"f.p.val")
so_mean1998_2021_pvalue <-apply(so_anom,c(1,2),"f.p.val")

# ----------------------------------------------------- #
### 6.1.2. Visualization of significant trends ----
# ----------------------------------------------------- #
par(mfrow=c(2, 2))
breakpoints <- c(0,0.05,1)
image.plot(thetao_mean1998_2021_pvalue,breaks=breakpoints,col=blue_to_red_palette(2), main="MEAN temperature")
image.plot(bottomT_mean1998_2021_pvalue,breaks=breakpoints,col=blue_to_red_palette(2), main="MEAN bottom temp")
image.plot(so_mean1998_2021_pvalue,breaks=breakpoints,col=blue_to_red_palette(2), main="MEAN salinity")

# ----------------------------------------------------- #
### 5.1.3. Transform to raster ----
# ----------------------------------------------------- #  
# Create the bounding box of the area of my interest 
bb_phy<-extent(min(lon_range_phy),max(lon_range_phy),min(lat_range_phy),max(lat_range_phy))

# Transform each variable to raster
# Temperature
thetao_mean1998_2021_pvalue_raster <- raster(thetao_mean1998_2021_pvalue) # transform the matrix into a raster
thetao_mean1998_2021_pvalue_raster <- t(flip(thetao_mean1998_2021_pvalue_raster, direction="x")) #to rotate the raster
projection(thetao_mean1998_2021_pvalue_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
thetao_mean1998_2021_pvalue_raster<-setExtent(thetao_mean1998_2021_pvalue_raster,bb_phy,keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
thetao_mean1998_2021_pvalue_raster

# Bottom Temperature
bottomT_mean1998_2021_pvalue_raster <- raster(bottomT_mean1998_2021_pvalue) # transform the matrix into a raster
bottomT_mean1998_2021_pvalue_raster <- t(flip(bottomT_mean1998_2021_pvalue_raster, direction="x")) #to rotate the raster
projection(bottomT_mean1998_2021_pvalue_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
bottomT_mean1998_2021_pvalue_raster<-setExtent(bottomT_mean1998_2021_pvalue_raster,bb_phy,keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
bottomT_mean1998_2021_pvalue_raster

# Salinity
so_mean1998_2021_pvalue_raster <- raster(so_mean1998_2021_pvalue) # transform the matrix into a raster
so_mean1998_2021_pvalue_raster <- t(flip(so_mean1998_2021_pvalue_raster, direction="x")) #to rotate the raster
projection(so_mean1998_2021_pvalue_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
so_mean1998_2021_pvalue_raster<-setExtent(so_mean1998_2021_pvalue_raster,bb_phy,keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
so_mean1998_2021_pvalue_raster

# ----------------------------------------------------- #
### 5.1.4. Transform to shapefiles ----
# ----------------------------------------------------- #
# Temperature
# Create the reclassified table
m <- c(cellStats(thetao_mean1998_2021_pvalue_raster, min), 0.05, 0, 0.05 , cellStats(thetao_mean1998_2021_pvalue_raster, max),  1)
rcl <- matrix(m, ncol=3, byrow=TRUE) #reclass matrix
# Reclassification of the raster
thetao_mean1998_2021_pvalue_shp <-reclassify(thetao_mean1998_2021_pvalue_raster, rcl, na.rm=TRUE)
# Transform to polygon
thetao_mean1998_2021_pvalue_shp <- rasterToPolygons(thetao_mean1998_2021_pvalue_shp, fun=function(x){x==1}, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)

# Bottom Temperature
# Create the reclassified table
m <- c(cellStats(bottomT_mean1998_2021_pvalue_raster, min), 0.05, 0, 0.05 , cellStats(bottomT_mean1998_2021_pvalue_raster, max),  1)
rcl <- matrix(m, ncol=3, byrow=TRUE) #reclass matrix
# Reclassification of the raster
bottomT_mean1998_2021_pvalue_shp<-reclassify(bottomT_mean1998_2021_pvalue_raster, rcl, na.rm=TRUE)
# Transform to polygon
bottomT_mean1998_2021_pvalue_shp <- rasterToPolygons(bottomT_mean1998_2021_pvalue_shp, fun=function(x){x==1}, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)

# Salnity
# Create the reclassified table
m <- c(cellStats(so_mean1998_2021_pvalue_raster, min), 0.05, 0, 0.05 , cellStats(so_mean1998_2021_pvalue_raster, max),  1)  
rcl <- matrix(m, ncol=3, byrow=TRUE) #reclass matrix
# Reclassification of the raster
so_mean1998_2021_pvalue_shp <-reclassify(so_mean1998_2021_pvalue_raster, rcl, na.rm=TRUE)
# Transform to polygon
so_mean1998_2021_pvalue_shp <- rasterToPolygons(so_mean1998_2021_pvalue_shp, fun=function(x){x==1}, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)


# ----------------------------------------------------- #  
## 5.2. PVALUE Product 2 ----
# ----------------------------------------------------- #  	

# -------------------------------- #
### 5.2.1. Apply the function ----
# -------------------------------- # 	
chl_mean1998_2021_pvalue <-apply(chl_anom,c(1,2),"f.p.val")
nppv_mean1998_2021_pvalue <-apply(nppv_anom,c(1,2),"f.p.val")
o2_mean1998_2021_pvalue <-apply(o2_anom,c(1,2),"f.p.val")

# ----------------------------------------------------- #
### 5.2.2. Visualization of significant trends ----
# ----------------------------------------------------- #
par(mfrow=c(2,2))
breakpoints <- c(0,0.05,1)
image.plot(chl_mean1998_2021_pvalue,breaks=breakpoints,col=blue_to_red_palette(2), main="MEAN chl 1998-2021")
image.plot(nppv_mean1998_2021_pvalue,breaks=breakpoints,col=blue_to_red_palette(2), main="MEAN nppv 1998-2021")
image.plot(o2_mean1998_2021_pvalue,breaks=breakpoints,col=blue_to_red_palette(2), main="MEAN O2")

# ----------------------------------------------------- #
### 5.2.3. Transform to raster ----
# ----------------------------------------------------- #	
# Create the bounding box of the area of my interest
bb_bio<-extent(min(lon_range_bio),max(lon_range_bio),min(lat_range_bio),max(lat_range_bio))
  	
# Transform each variable to raster
# Chl
chl_mean1998_2021_pvalue_raster <- raster(chl_mean1998_2021_pvalue) # transform the matrix into a raster
chl_mean1998_2021_pvalue_raster <- t(flip(chl_mean1998_2021_pvalue_raster, direction="x")) #to rotate the raster
projection(chl_mean1998_2021_pvalue_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
chl_mean1998_2021_pvalue_raster<-setExtent(chl_mean1998_2021_pvalue_raster,bb_bio, keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
chl_mean1998_2021_pvalue_raster

# NPP
nppv_mean1998_2021_pvalue_raster <- raster(nppv_mean1998_2021_pvalue) # transform the matrix into a raster
nppv_mean1998_2021_pvalue_raster <- t(flip(nppv_mean1998_2021_pvalue_raster, direction="x")) #to rotate the raster
projection(nppv_mean1998_2021_pvalue_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
nppv_mean1998_2021_pvalue_raster<-setExtent(nppv_mean1998_2021_pvalue_raster,bb_bio, keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
nppv_mean1998_2021_pvalue_raster

# Oxygen
o2_mean1998_2021_pvalue_raster <- raster(o2_mean1998_2021_pvalue) # transform the matrix into a raster
o2_mean1998_2021_pvalue_raster <- t(flip(o2_mean1998_2021_pvalue_raster, direction="x")) #to rotate the raster
projection(o2_mean1998_2021_pvalue_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
o2_mean1998_2021_pvalue_raster<-setExtent(o2_mean1998_2021_pvalue_raster,bb_bio,keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
o2_mean1998_2021_pvalue_raster

# ----------------------------------------------------- #
### 5.2.4. Transform to shapefiles ----
# ----------------------------------------------------- #
  	
# Chl 
# Create the reclassified table
m <- c(cellStats(chl_mean1998_2021_pvalue_raster, "min"), 0.05, 0, 0.05 , cellStats(chl_mean1998_2021_pvalue_raster, max),  1)
rcl <- matrix(m, ncol=3, byrow=TRUE) #reclass matrix
rcl
# Reclassification of the raster
chl_mean1998_2021_pvalue_shp <- reclassify(chl_mean1998_2021_pvalue_raster, rcl, na.rm=TRUE)
chl_mean1998_2021_pvalue_shp
# Transform to polygon
chl_mean1998_2021_pvalue_shp <- rasterToPolygons(chl_mean1998_2021_pvalue_shp, fun=function(x){x==1}, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)

# NPP
# Create the reclassified table
m <- c(cellStats(nppv_mean1998_2021_pvalue_raster, "min"), 0.05, 0, 0.05 , cellStats(nppv_mean1998_2021_pvalue_raster, max),  1)  
rcl <- matrix(m, ncol=3, byrow=TRUE) #reclass matrix
# Reclassification of the raster
nppv_mean1998_2021_pvalue_shp <- reclassify(nppv_mean1998_2021_pvalue_raster, rcl, na.rm=TRUE)
# Transform to polygon
nppv_mean1998_2021_pvalue_shp <- rasterToPolygons(nppv_mean1998_2021_pvalue_shp, fun=function(x){x==1}, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)

# Oxygen
# Create the reclassified table
m <- c(cellStats(o2_mean1998_2021_pvalue_raster, "min"), 0.05, 0, 0.05 , cellStats(o2_mean1998_2021_pvalue_raster, max),  1)  
rcl <- matrix(m, ncol=3, byrow=TRUE) #reclass matrix
# Reclassification of the raster
o2_mean1998_2021_pvalue_shp <- reclassify(o2_mean1998_2021_pvalue_raster, rcl, na.rm=TRUE)
# Transform to polygon
o2_mean1998_2021_pvalue_shp <- rasterToPolygons(o2_mean1998_2021_pvalue_shp, fun=function(x){x==1}, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)

# ----------------------------------------------------- #  
## 5.3. Plot maps with p-values overlaid----
# ----------------------------------------------------- #  

par(mfrow=c(2,2))
par(cex.axis = 1.5)
transparent_grey <- rgb(0.5, 0.5, 0.5, alpha = 0.5)  # Grey color with 50% transparency

# Temperature
plot(norm_thetao_mean1998_2021_slopes_raster, 
     main="Temperature Trend 1998-2021",
     col = magenta_to_blue_palette(100),  
     axes = TRUE)

# Plot the shapefile with transparency (non-sig values)
plot(thetao_mean1998_2021_pvalue_shp, add = TRUE,     
     lwd = 0.5,                
     col = transparent_grey,   
     lty = 1                   
)

# Bottom Temperature
plot(norm_bottomT_mean1998_2021_slopes_raster, 
     main="Bottom Temperature Trend 1998-2021",
     col = magenta_to_blue_palette(100),  
     axes = TRUE)

# Plot the shapefile with transparency (non-sig values)
plot(bottomT_mean1998_2021_pvalue_shp, add = TRUE, 
     lwd = 0.5,                
     col = transparent_grey,   
     lty = 1                  
)

# Salinity
plot(norm_so_mean1998_2021_slopes_raster, 
     main="Salinity Trend 1998-2021",
     col = magenta_to_blue_palette(100),  
     axes = TRUE)

# Plot the shapefile with transparency (non-sig values)
plot(so_mean1998_2021_pvalue_shp, add = TRUE, 
     lwd = 0.5,               
     col = transparent_grey,  
     lty = 1                  
)

# If we want to save the image (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/7.a.Physical_trend_normalized_&_pvalues.tiff", width = 800, height = 600, res = 300, units= "px")
#dev.off()

# Chlorophyll
par(mfrow=c(2,2))

plot(norm_chl_mean1998_2021_slopes_raster, 
     main="Chl Trend 1998-2021",
     col = magenta_to_blue_palette(100),  
     axes = TRUE)
# Plot the shapefile with transparency (non-sig values)
plot(chl_mean1998_2021_pvalue_shp, add = TRUE, 
     lwd = 0.5,                
     col = transparent_grey,   
     lty = 1                   
)

# NPP
plot(norm_nppv_mean1998_2021_slopes_raster, 
     main="NPP Trend 1998-2021",
     col = magenta_to_blue_palette(100),  
     axes = TRUE)
# Plot the shapefile with transparency (non-sig values)
plot(nppv_mean1998_2021_pvalue_shp, add = TRUE, 
     lwd = 0.5,                
     col = transparent_grey,  
     lty = 1                   
)

# Oxygen
plot(norm_o2_mean1998_2021_slopes_raster, 
     main="Oxygen Trend 1998-2021",
     col = magenta_to_blue_palette(100),  
     axes = TRUE)
# Plot the shapefile with transparency (non-sig values)
plot(o2_mean1998_2021_pvalue_shp, add = TRUE, 
     lwd = 0.5,                
     col = transparent_grey,   
     lty = 1    )           

# If we want to save the image (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_monthly_anomalies/7.b.BGC_trend_normalized_&_pvalues.tiff", width = 800, height = 600, res = 300, units= "px")
#dev.off()


# ----------------------------------------------------- #
#END ----
# ----------------------------------------------------- #   	
