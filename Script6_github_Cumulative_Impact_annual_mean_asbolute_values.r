#------------------------------------------------------------ #
#   Comparing methodological choices for environmental cumulative impacts analysis: The Black Sea as a case study
#   SCRIPT 6                                            
#   CUMULATIVE IMPACTS
#   Based on annual mean absolute values
#   BLACK SEA
#   SURFACE DATA (Yearly matrix as input)
#   Authors: Lucia Espasandin, Fran Ramirez and Elena Lloret-Lloret
#   Last update: 20/08/2025                                                         
#------------------------------------------------------------ #

# AIM OF THE SCRIPT

# Calculate a cumulative impact analysis based on yearly averages.
# We use a set of physical and a set of biogeochemical variables over a relatively large time range
# and calculate the slopes of change over time, per pixel.
# 

# Here we use SURFACE environmental variables as an example, 
# but any of the previous pre-process outputs can be used instead ( vertical integration weighted mean,vertical integration weighted sum, vertical integration 1m interpolation)


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
# This sets the working directory directly where your script is:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ----------------------------------------------------- #
## 0.2. Define parameters that need to be changed -----
# ---------------------------------------------------- #

# Determine the start and end year at the beginning of each script
# Your data needs to have complete years (12 months)!!

start_year<- 1998
end_year<- 2021

years_all<- seq(start_year, end_year, 1)

# var time for linear models
length(years_all) # these products have 29 years data
years_all # from 1993 to 2021
years_all <- as.Date(paste0(years_all, "-01-01"))
years_all
times_all <- seq(1, length(years_all), by=1) #we create a time object to ease the script
times_all

# Define longitude and latitude ranges for the analysis
# If you don't need a subset, you still need to define the min and max of your entire dataset
# These values must be EXACTLY the SAME as in the nc. (otherwise it won't work)

#Lon lat for the Physical variables
lon_range_phy <- c(27.37000, 41.96259)  #all black sea
lat_range_phy <- c(40.86000, 46.80444)  #all black_sea

#Lon lat for the Biological variables
lon_range_bio <- c(27.250, 42.000) #all black sea
lat_range_bio <- c(40.500, 47.000) #all black sea

#create colour palettes for layer use
blue_to_red_palette <- colorRampPalette(colors = c("blue", "beige", "red"))
white_to_red_palette <- colorRampPalette(colors = c("beige", "red"))
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

#-----------------------------------#
## 1.2 ENVIRONMENTAL PRODUCTS ---- 
#----------------------------------#

#----------------------------------#
### 1.2.1. LOAD PHYSICAL FILES ----
#----------------------------------# 
# Load
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# For this analysis we load the yearly matrix
# It goes from 1993-2021
thetao <-readRDS( "./Data/CI_annual_mean_asbolute_values/surface_temp_matrix_YEARLY_Black_Sea_1993_2021.RData")
str(thetao)
summary(thetao)
dim(thetao)

# Sub-setting the third dimension (years 1998 to 2021)
thetao_subset <- thetao[, , 6:29]
dim(thetao_subset)

# Calculate and visualize the spatial mean of all the dataset
thetao_mean1998_2021 <- apply(thetao_subset,c(1,2),"mean")
image.plot(thetao_mean1998_2021, main="Mean surface Temp 98-21")

# Sea water potential temperature at sea floor: bottomT [°C]
bottomT<-readRDS( "./Data/CI_annual_mean_asbolute_values/temp_bot_matrix_YEARLY_Black_Sea_1993_2021.RData")
str(bottomT)
summary(bottomT)
dim(bottomT)

# Sub-setting the third dimension (years 1998 to 2021)
bottomT_subset <- bottomT[, , 6:29]
dim(bottomT_subset)

# Calculate and visualize the spatial mean of all the dataset
bottomT_mean1998_2021 <- apply(bottomT_subset,c(1,2),"mean")
image.plot(bottomT_mean1998_2021, main="Mean bottom Temp 98-21")

# Sea water salinity: so [10-3]
so<-readRDS( "./Data/CI_annual_mean_asbolute_values/surface_salinity_matrix_YEARLY_Black_Sea_1993_2021.RData")
str(so)
summary(so)
dim(so)

# Sub-setting the third dimension (years 1998 to 2021)
so_subset <- so[, , 6:29]
dim(so_subset)

# Calculate and visualize the spatial mean of all the dataset
so_mean1998_2021 <- apply(so_subset,c(1,2),"mean")
image.plot(so_mean1998_2021, main="Mean Salinity 98-21")

#-------------------------------------------#
### 1.2.2. LOAD BIOGEOCHEMICAL FILES ----
#------------------------------------------#
# Chlorophyll
chl<-readRDS( "./Data/CI_annual_mean_asbolute_values/surface_chl_matrix_YEARLY_Black_Sea_1992_2022.RData")
str(chl)
summary(chl)
dim(chl) 
# Sub-setting the third dimension (years 1998 to 2021)
chl_subset <- chl[, , 7:30]
dim(chl_subset)

# Calculate and visualize the spatial mean of all the dataset
chl_mean1998_2021 <- apply(chl_subset,c(1,2),"mean")
image.plot(chl_mean1998_2021, main="Mean Surface Chl 1998-2021")

#NPP
nppv<-readRDS( "./Data/CI_annual_mean_asbolute_values/surface_npp_matrix_YEARLY_Black_Sea_1992_2022_2.RData")
str(nppv)
summary(nppv)
dim(nppv)

# Sub-setting the third dimension (years 1998 to 2021)
nppv_subset <- nppv[, , 7:30]
dim(nppv_subset)

# Calculate and visualize the spatial mean of all the dataset
nppv_mean1998_2021 <- apply(nppv_subset,c(1,2),"mean")
image.plot(nppv_mean1998_2021, main="Mean Surface nppv 1998-2021")

#Oxygen
o2<-readRDS( "./Data/CI_annual_mean_asbolute_values/surface_oxygen_matrix_YEARLY_Black_Sea_1992_2022.RData")
str(o2)
summary(o2)
dim(o2)

# Sub-setting the third dimension (years 1998 to 2021)
o2_subset <- o2[, , 7:30]
dim(o2_subset)

#Calculate and visualize the spatial mean of all the dataset
o2_mean1998_2021 <- apply(o2_subset,c(1,2),"mean")
image.plot(o2_mean1998_2021, main="Mean Surface oxygen 1998-2021")



#------------------------------#
# 2. PLOT OVERALL TRENDS---- 
#------------------------------#	

#--------------------------------#
## 2.1. Physical variables ----
#----------------------------------#

# Temperature
# Calculate the mean (of all the area) across all years for each year
thetao_spatial_mean <- apply(thetao_subset, c(3), mean, na.rm = TRUE)  # Averaging over the first two dimensions (116 and 128)

# Prepare data to plot 
years <- seq(1998, 2021)
# Create a dataframe
thetao_df <- data.frame(Years = years, Temperature = thetao_spatial_mean)
head(thetao_df)

# Plot
ggplot() +
  geom_line(data=thetao_df, aes(x = Years, y = Temperature), color = '#d53e4f', size = 1.5) +
  labs(title = "Surface Temperature (1998-2021)",
       x = "Years",
       y = "Temperature (ºC)") +
  theme_minimal()


# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Yearly_surface_temperature_1998-2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

# Bottom Temperature
# Calculate the mean (of all the area) across all years for each year
bottomT_spatial_mean <- apply(bottomT_subset, c(3), mean, na.rm = TRUE)  # Averaging over the first two dimensions (116 and 128)

# Prepare data to plot 
years <- seq(1998, 2021)
# Create a dataframe
bottomT_df <- data.frame(Years = years, Temperature = bottomT_spatial_mean)
head(bottomT_df)

# Plot
ggplot() +
  geom_line(data=bottomT_df, aes(x = Years, y = Temperature), color = '#d53e4f', size = 1.5) +
  labs(title = "Bottom Temperature (1998-2021)",
       x = "Years",
       y = "Temperature (ºC)") +
  theme_minimal()


# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Yearly_bottom_temperature_1998-2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

# Salinity
# Calculate the mean (of all the area) across all years for each year
so_spatial_mean <- apply(so_subset, c(3), mean, na.rm = TRUE)  # Averaging over the first two dimensions (116 and 128)

# Prepare data to plot 
years <- seq(1998, 2021)
# Create a dataframe
so_df <- data.frame(Years = years, Salinity = so_spatial_mean)
head(so_df)

# Plot
ggplot() +
  geom_line(data=so_df, aes(x = Years, y = Salinity), color = '#d53e4f', size = 1.5) +
  labs(title = "Surface Salinity (1998-2021)",
       x = "Years",
       y = "Salinity (PSU)") +
  theme_minimal()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Yearly_surface_salinity_1998-2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

#-------------------------------------#
## 2.2. Biogeochemical variables ----
#------------------------------------#

# Chlorophyll-a
# Calculate the mean (of all the area) across all years for each year
chl_spatial_mean <- apply(chl_subset, c(3), mean, na.rm = TRUE)  # Averaging over the first two dimensions (116 and 128)

# Prepare data to plot 
years <- seq(1998, 2021)
# Create a dataframe
chl_df <- data.frame(Years = years, Chl = chl_spatial_mean)
head(chl_df)

# Plot
ggplot() +
  geom_line(data=chl_df, aes(x = Years, y = Chl), color = '#d53e4f', size = 1.5) +
  labs(title = "Surface Chlrophyll-a (1998-2021)",
       x = "Years",
       y = "Chlorophyll") +
  theme_minimal()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Yearly_surface_chl_1998-2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()


# NPP
# Calculate the mean (of all the area) across all years for each year
nppv_spatial_mean <- apply(nppv_subset, c(3), mean, na.rm = TRUE)  # Averaging over the first two dimensions (116 and 128)

# Prepare data to plot 
years <- seq(1998, 2021)
# Create a dataframe
nppv_df <- data.frame(Years = years, NPP = nppv_spatial_mean)
head(nppv_df)

# Plot
ggplot() +
  geom_line(data=nppv_df, aes(x = Years, y = NPP), color = '#d53e4f', size = 1.5) +
  labs(title = "Surface NPP (1998-2021)",
       x = "Years",
       y = "NPP") +
  theme_minimal()

# If we want to save the plot (as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Yearly_surface_NPP_1998-2021.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()



#------------------------------#
# 3. CALCULATE THE SLOPE    ---- 
#------------------------------#	

#------------------------------------#
## 3.1. Slope for PHYSICAL FILES ----
#------------------------------------#

# Create the function that calculates the rate of change per time (slope)
fun2 <- function(x) 			{ 
								if (is.na(x[1])) {
								NA 
								} else { 
								m = lm(x ~ times_all) 
								m$coefficients[2]
								} 
								} 


#---------------------------------#
### 3.1.1. Apply the function ----
#---------------------------------#
# Applying fun2 (slope calculation) 
# In apply(): c(1,2) indicates that the calculation is performed for both rows and columns
# It calculates the regression of the value over the 3rd dimension which is time
thetao_mean1998_2021_slopes <- apply(thetao_subset,c(1,2),"fun2")
bottomT_mean1998_2021_slopes <- apply(bottomT_subset,c(1,2),"fun2")  
so_mean1998_2021_slopes <- apply(so_subset,c(1,2),"fun2")

#--------------------------------------#  
### 3.1.2. Visualization of arrays ----
#--------------------------------------#
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
     breaks = 50, 
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
     breaks = 50, 
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
     breaks = 50, 
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
#Temperature
thetao_mean1998_2021_slopes_raster <- raster(thetao_mean1998_2021_slopes) # transform the matrix into a raster
thetao_mean1998_2021_slopes_raster <- t(flip(thetao_mean1998_2021_slopes_raster, direction="x")) #to rotate the raster
projection(thetao_mean1998_2021_slopes_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
thetao_mean1998_2021_slopes_raster<-setExtent(thetao_mean1998_2021_slopes_raster,bb_phy,keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
thetao_mean1998_2021_slopes_raster

#Bottom temperature
bottomT_mean1998_2021_slopes_raster <- raster(bottomT_mean1998_2021_slopes) # transform the matrix into a raster
bottomT_mean1998_2021_slopes_raster <- t(flip(bottomT_mean1998_2021_slopes_raster, direction="x")) #to rotate the raster
projection(bottomT_mean1998_2021_slopes_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
bottomT_mean1998_2021_slopes_raster<-setExtent(bottomT_mean1998_2021_slopes_raster,bb_phy,keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
bottomT_mean1998_2021_slopes_raster
	
#Salinity
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

#Plot
par(cex.axis = 1.5)  # Increase the axis text size
plot(stacked_raster_physical, 
     col = magenta_to_blue_palette(100),  
     axes = TRUE)

# If we want to save the image Physical( as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/1.a.Physical_trend.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

#-------------------------------#
### 3.1.5. Save the rasters ----
#------------------------------#

#writeRaster(thetao_mean1998_2021_slopes_raster, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Rasters/surface_thetao_mean1998_2021_slopes_raster.tif", sep="", overwrite=TRUE) #geotif
#writeRaster(bottomT_mean1998_2021_slopes_raster, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Rasters/bottomT_mean1998_2021_slopes_raster.tif", sep="", overwrite=TRUE) #geotif
#writeRaster(so_mean1998_2021_slopes_raster, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Rasters/surface_so_mean1998_2021_slopes_raster.tif", sep="", overwrite=TRUE) #geotif

#---------------------------------------#
## 3.2. SLOPE  biogeochemical data ----
#---------------------------------------#

# Create the function that calculates the rate of change per time (slope)
# In this case fun2 and fun3 are identical, but if the physical and biogeochemical datasets have different times, the function have to be adadted respectively
fun3 <- function(x) 			{ 
  if (is.na(x[1])) {
    NA 
  } else { 
    m = lm(x ~ times_all) 
    m$coefficients[2]
  } 
	} 

#---------------------------------------#
### 3.2.1. Apply the function  ----
#---------------------------------------#
# Applying fun3 (slope calculation) 
# In apply(): c(1,2) indicates that the calculation is performed for both rows and columns
chl_mean1998_2021_slopes <-apply(chl_subset,c(1,2),"fun3")
nppv_mean1998_2021_slopes <-apply(nppv_subset,c(1,2),"fun3")
o2_mean1998_2021_slopes <-apply(o2_subset,c(1,2),"fun3")

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

# Oxygen
o2_mean1998_2021_slopes_raster <- raster(o2_mean1998_2021_slopes) # transform the matrix into a raster
o2_mean1998_2021_slopes_raster <- t(flip(o2_mean1998_2021_slopes_raster, direction="x")) #to rotate the raster
projection(o2_mean1998_2021_slopes_raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #define projection
o2_mean1998_2021_slopes_raster<-setExtent(o2_mean1998_2021_slopes_raster,bb_bio,keepres=FALSE, snap=FALSE) # apply the bounding box to the raster
o2_mean1998_2021_slopes_raster

#---------------------------------------#
### 3.2.4. Visualization of rasters ----
#---------------------------------------#
par(mfrow=c(2,2))
plot(chl_mean1998_2021_slopes_raster, col = magenta_to_blue_palette(100), main="Raster Surfae Chl 1998-2021")
plot(nppv_mean1998_2021_slopes_raster, col = magenta_to_blue_palette(100), main="Raster Surface NPP 1998-2021")
plot(o2_mean1998_2021_slopes_raster, col = magenta_to_blue_palette(100), main="Raster Surface Oxygen 1998-2021")

#Stack biogeochemical variables:
stacked_raster_biogeo <- stack(chl_mean1998_2021_slopes_raster, 
  	                         nppv_mean1998_2021_slopes_raster,
                           o2_mean1998_2021_slopes_raster)
# Rename
names(stacked_raster_biogeo) = c("chl", "nppv", "o2")

# Plot
par(cex.axis = 1.5)  # Increase the axis text size
plot(stacked_raster_biogeo, 
     col = magenta_to_blue_palette(100),  
     axes = TRUE)

# If we want to save the image BGC( as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/1.b.BGC_trend.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

#---------------------------------------#	
### 3.2.5. Save the rasters ----
#---------------------------------------#
#writeRaster(chl_mean1998_2021_slopes_raster, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Rasters/surface_chl_mean1998_2021_slopes_raster.tif", sep="", overwrite=TRUE) #geotif
#writeRaster(nppv_mean1998_2021_slopes_raster, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Rasters/surface_nppv_mean1998_2021_slopes_raster.tif", sep="", overwrite=TRUE) #geotif
#writeRaster(o2_mean1998_2021_slopes_raster,  "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Rasters/surface_o2_mean1998_2022_slopes_raster.tif", sep="", overwrite=TRUE) #geotif

#------------------------------#
# 4. CUMULATIVE IMPACTS    ---- 
#------------------------------#

#---------------------------------------#
## 4.1. CUMULATIVE IMPACTS Product 1 ----
#---------------------------------------#

#---------------------------------------#
### 4.1.1. Calculate absolute values ----
#---------------------------------------#
### We consider as change all the slope values (positive or negative) 
abs_thetao_mean1998_2021_slopes_raster <- abs(thetao_mean1998_2021_slopes_raster)
abs_bottomT_mean1998_2021_slopes_raster <- abs(bottomT_mean1998_2021_slopes_raster)
abs_so_mean1998_2021_slopes_raster <- abs(so_mean1998_2021_slopes_raster)

#----------------------------------------------------#
### 4.1.2. Normalization Physical Impacts ----
#-------------------------------------------------------#
# We use the 10th and 90th percentiles based normalization. To reduce the impact of extreme values

#Temperature
par(mfrow=c(2,2))
# Calculate quantiles
p10_value_thetao <-quantile(values(abs_thetao_mean1998_2021_slopes_raster), probs = 0.02, na.rm = TRUE)
p90_value_thetao <- quantile(values(abs_thetao_mean1998_2021_slopes_raster), probs = 0.98, na.rm = TRUE)
# Plot histogram
hist(values(abs_thetao_mean1998_2021_slopes_raster), 
     breaks = 50, main = "Histogram of Raster Values",
     xlab = "Values", col = "lightblue", border = "black")
#ablines
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
#abline
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

# If we want to save the image Physical( as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/2.a.Physical_trend_normalized.tiff", width = 10, height = 10, units = "in", res = 300)
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
#ablines
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
#Plot
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
#ablines
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
#Plot
plot(norm_nppv_mean1998_2021_slopes_raster, col= magenta_to_blue_palette(100), main="NPP Trend 1998-2021")

# Oxygen
par(mfrow=c(2,2))
# Calculate quantiles
p10_value_o2 <-quantile(values(abs_o2_mean1998_2021_slopes_raster), probs = 0.02, na.rm = TRUE)
p90_value_o2 <- quantile(values(abs_o2_mean1998_2021_slopes_raster), probs = 0.98, na.rm = TRUE)
# Plot histogram
hist(values(abs_o2_mean1998_2021_slopes_raster), 
     breaks = 50, main = "Histogram of Raster Values Oxygen",
     xlab = "Values", col = "lightblue", border = "black")
#ablines
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
#Plot
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

# Pot
par(cex.axis = 1.5)
plot(stacked_raster2, 
     col = magenta_to_blue_palette(100),  
     axes = TRUE)

# If we want to save the image BGC( as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/2.b.BGC_trend_normalized.tiff", width = 10, height = 10, units = "in", res = 300)
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

stand_BGC_impact #(x,y)
stand_phy_impact #(x,y)

# Resample
# Aways resample to the coarser resolution
stand_BGC_impact_res <- resample(stand_BGC_impact, stand_phy_impact, method="bilinear", na.rm=TRUE)

# Check resampled object
stand_BGC_impact_res
summary(stand_BGC_impact_res)
plot(stand_BGC_impact_res,col=white_to_red_palette(100), main= "All BGC impacts" )

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

# Print and plot the standardized values
print(stand_BGC_impact_final)
plot(stand_BGC_impact_final, col=white_to_red_palette(100), main= "All bgc impacts rescaled" )


# ----------------------------------------------------- #
### 4.3.2. create a unique layer of Physical impacts----
# ----------------------------------------------------- #

# Add all the layers of one raster stack
par(mfrow=c(2,2))
CI_phy <- sum(stand_phy_impact) # sum of all the RasterStacks

# Plot
plot(CI_phy, col=magenta_to_blue_palette(100), main= "Physical cumulative impacts")

# Save raster
#writeRaster(CI_phy, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Rasters/surface_CI_sum_physical_1998-2021_norm2_98.tif", overwrite=TRUE) #geotif

# ----------------------------------------------------- #
### 4.3.3. Create a unique layer of BGC impacts----
# ----------------------------------------------------- #
# Add all the layers of one raster stack
CI_BGC <- sum(stand_BGC_impact_final) # sum of all the RasterStacks

# Plot
plot(CI_BGC, col=magenta_to_blue_palette(100), main= "Biogeochemical Cumulative Impacts")

# Save Raster
#writeRaster(CI_BGC, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Rasters/surface_CI_sum_BGC_1998-2021_norm_2_98.tif", overwrite=TRUE) #geotif

# If we want to save the image with sum physical and sum BGC( as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/3.Physical_sum_&_BGC_sum.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()
# ---------------------------------------------------------------- #
### 4.3.4. Create a unique layer of Physical and BGC impacts----
# ---------------------------------------------------------- #

# Add the two stacks of variables phy and bio together
CI_sum <- CI_phy + CI_BGC

# Plot
par(mfrow=c(1,1))
plot(CI_sum,col=magenta_to_blue_palette(100), main= "All Cumulative Impacts" )

# Save raster
 #writeRaster(CI_sum,  "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Rasters/surface_CI_sum_1998-2021_norm_2_98.tif", sep="", overwrite=TRUE) #geotif

# If we want to save the image with sum of all impacts (physical and BGC together ( as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/4.All_sum.tiff", width = 10, height = 10, units = "in", res = 300)
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
rcl_phy <- matrix(rcl_phy, ncol=3, byrow=TRUE)
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
rcl_BGC <- c(0,quantiles2[1], 1, quantiles2[1], 
             quantiles2[2], 2, quantiles2[2], 
             quantiles2[3], 3, quantiles2[3], 
             cellStats(CI_BGC, stat=max), 4)
# Transform into a matrix
rcl_BGC <- matrix(rcl_BGC, ncol=3, byrow=TRUE) #reclass matrix
rcl_BGC

# Reclassification of the original raster (CI_phy)into a 1 to 4 numerations
bgc_qtile_hab <- reclassify(CI_BGC, rcl_BGC, na.rm=TRUE) 

#Plot
plot(bgc_qtile_hab, col=magenta_to_blue_palette(100), main="BGC quartiles")

# If we want to save the physical and BGC quantiles( as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/5.Physical_&_BGC_quantiles.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

# ----------------------------------------------------- #
## 4.5. FINAL CUMULATIVE IMPACT MAP----
# ----------------------------------------------------- #
# ----------------------------------------------------- #
### 4.5.1. Visualization bivariate plot ----
# ----------------------------------------------------- #
    
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

# Convert quantiles (physical and biogeochemical impacts) to data frames (x=long, y=lat)
phy_df <- as.data.frame(phy_qtile_hab, xy = TRUE, na.rm=TRUE) %>% 
  setNames(c('x', 'y', 'phy'))
summary(phy_qtile_hab)
unique(phy_df$phy)

bgc_df <- as.data.frame(bgc_qtile_hab, xy = TRUE, na.rm=TRUE) %>%
  setNames(c('x', 'y', 'bgc'))
summary(bgc_qtile_hab)
unique(bgc_df$bgc) 

# Merge PHY and BGC quantiles data frames
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
    
# Assamble the figure (map + legend)
panel_a_map <- cowplot::get_panel(biplot_hab2)
panel_a_lgd <- cowplot::get_panel(panel_lgd_raw)
lgd_p <- c(x = 0.75, y = .8, h = .15, w = .15) # position of the legend right corner
#lgd_p <- c(x = 0.05, y = .85, h = .15, w = .15) # position of the legend left corner

#
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

# # As .tiff
# dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/6.a.CI_based_on_mean_annual_values_surface.tiff", width = 800, height = 600, res = 300, units= "px")
# dev.off()

# ----------------------------------------------------- #
### 4.5.2. Rasterize Cumulative impacts map ----
# ----------------------------------------------------- #

# Ensure bi_class is a factor and convert to integer codes
phy_bgc_df$bi_class <- as.factor(phy_bgc_df$bi_class)
# The function raster() cant rasterize a filed in the form of 1-1/1-2/4-4, etc so we need to create a single integer for eache code
phy_bgc_df$bi_class_int <- as.integer(phy_bgc_df$bi_class)
str(phy_bgc_df)

# Check which integer has been assigned to each number code
unique_combinations <- unique(phy_bgc_df[c("bi_class", "bi_class_int")])
unique_combinations <- unique_combinations[order(unique_combinations$bi_class), ]
print(unique_combinations)

# Create an empty raster first with the same extent and resolution as the original rasters
CI_r <- raster(extent(phy_qtile_hab), nrow = nrow(phy_qtile_hab), ncol = ncol(phy_qtile_hab))
crs(CI_r) <- crs(phy_qtile_hab)
CI_r

# Convert the data frame to a SpatialPointsDataFrame
phy_bgc_spatial<-phy_bgc_df
coordinates(phy_bgc_spatial) <- ~x + y

# Rasterize the bi_class_int column from the phy_bgc_df
CI_raster <- rasterize(phy_bgc_spatial, CI_r, field="bi_class_int")
CI_raster

# ---------------------------------------------------------- #
#### 4.5.2.1. Visaulization as Raster Cumulative impacts ----
# ---------------------------------------------------------- #

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

# Copy the plot to a TIFF file
# If we want to save the image for ( as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/6.b.CI_based_on_mean_annual_values_surface_countries.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()

# ----------------------------------------------------- #
#### 4.5.2.2. Save final CI raster ----
# ----------------------------------------------------- #

# Save raster
#writeRaster(CI_raster, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Rasters/CI_based_on_mean_annual_values_surface.tif", overwrite=TRUE) #geotif

# ----------------------------------------------------- #
### 4.5.3.Transform to Polygon----
# ----------------------------------------------------- #

# In case this is needed for a posterior spatial analysis.
# Convert the raster to a SpatRaster object
CI_terra <- rast(CI_raster)

# Convert the raster to polygons
CI_polygons <- as.polygons(CI_terra)

# -------------------------------- #
#### 4.5.3.1. Plot as Polygon----
# ------------------------------- #

# Plot the polygons using the layer attribute for colors
plot(CI_polygons, col =color_palette[levels(factor(phy_bgc_df$bi_class))], main = "Polygons representing raster levels", )
spplot(CI_polygons, col.regions = color_palette[levels(factor(phy_bgc_df$bi_class))], main = "Polygons representing raster levels")

# ------------------------------- #
#### 4.5.3.2. Save as Polygon----
# ------------------------------- #

# Write the polygons to a shapefile
# Convert the polygons to an sf object first
sf_polygons <- st_as_sf(CI_polygons)

# Write the shapefile
#st_write(sf_polygons, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/Polygon/CI_based_on_mean_annual_values_surface.shp")

#------------------------------#
# 5. CALCULATE THE P VALUE  ---- 
#------------------------------#	

# Create the function to extract the overall ANOVA p-value out of a linear model object

#
f.p.val <- function(x){
  if(sum(!is.na(x)) < 10) {
    p.val <- NA
  }
  if(sum(!is.na(x)) >= 10) {
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
thetao_mean1998_2021_pvalue <-apply(thetao_subset,c(1,2),"f.p.val")
bottomT_mean1998_2021_pvalue <-apply(bottomT_subset,c(1,2),"f.p.val")
so_mean1998_2021_pvalue <-apply(so_subset,c(1,2),"f.p.val")


# ----------------------------------------------------- #
### 6.1.2. Visualization of significant trends ----
# ----------------------------------------------------- #
par(mfrow=c(2, 2))
breakpoints <- c(0,0.05,1)
image.plot(bottomT_mean1998_2021_pvalue,breaks=breakpoints,col=blue_to_red_palette(2), main="MEAN bottom temp")
image.plot(so_mean1998_2021_pvalue,breaks=breakpoints,col=blue_to_red_palette(2), main="MEAN salinity")
image.plot(thetao_mean1998_2021_pvalue,breaks=breakpoints,col=blue_to_red_palette(2), main="MEAN temperature")

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

# Salinity
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

# ----------------------------------------------------- #
### 5.2.1. Apply the function ----
# ----------------------------------------------------- # 	

chl_mean1998_2021_pvalue <-apply(chl_subset,c(1,2),"f.p.val")
nppv_mean1998_2021_pvalue <-apply(nppv_subset,c(1,2),"f.p.val")
o2_mean1998_2021_pvalue <-apply(o2_subset,c(1,2),"f.p.val")

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
  	
#Transform each variable to raster
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

#Oxygen
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
par(cex.axis = 1.5)  # Increase the axis text size
transparent_grey <- rgb(0.5, 0.5, 0.5, alpha = 0.5)  # Grey color with 50% transparency

#Temperature
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

#Bottom Temperature
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

#Salinity
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

# If we want to save the image Physical( as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/7.a.Physical_trend_normalized_&_pvalues.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()


#Chlorophyll
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

#NPP
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

#Oxygen
plot(norm_o2_mean1998_2021_slopes_raster, 
     main="Oxygen Trend 1998-2021",
     col = magenta_to_blue_palette(100),  
     axes = TRUE)
# Plot the shapefile with transparency (non-sig values)
plot(o2_mean1998_2021_pvalue_shp, add = TRUE, 
     lwd = 0.5,                
     col = transparent_grey,   
     lty = 1    )             

# If we want to save the image for BGC ( as.tif)
#dev.copy(tiff, "./Outputs/CI_analysis/CI_annual_mean_asbolute_values/7.b.BGC_trend_normalized_&_pvalues.tiff", width = 10, height = 10, units = "in", res = 300)
#dev.off()


# ----------------------------------------------------- #
#END ----
# ----------------------------------------------------- #   	
