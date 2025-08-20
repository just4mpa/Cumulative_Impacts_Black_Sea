#------------------------------------------------------------#
#   Comparing methodological choices for environmental cumulative impacts analysis: The Black Sea as a case study
#   Script 1
#   Pre-process: Prepare input data
#   Surface data
#   Authors: Elena Lloret-Lloret & Camila Artana
#   Last update: 20/08/2025
# ------------------------------------------------------------ # 

# AIM OF THE SCRIPT

# Pre-process netcfd files to obtain surface data.
# With this scripts we generate surface matrix per year, and  per month. 
# This files will be used as input data in the different version of the Cumulative Impacts analysis (based on anual means and on monthly anomalies, respectively). 
# Here, this script is applied for Chlorophyll-a but it can be used for all the variables we want to include in the cumulative impacts analysis.

# Accommodate all the necessary parameters to your area in the first section of the script
# Please keep in mind that your data needs to have complete years (12 months) for the script to work properly.
# You will also need the latitude and longitude information of your study area.
# Your data has to be organised in the following manner:
#           -A folder containing .nc files of the variable of interest per year with  a monthly resolution
#           (if you are going to run the code for several variables, have one folder per variable)


# --------------- #
# 0. SET UP -----
# --------------- #

# ---------------------- #
## 0.1. Load packages-----
# ---------------------- #
# 
# install.packages(c(
#   "raster", "terra", "fields", "plyr", "ggplot2", "sp", "sf",
#   "ncdf4", "lubridate", "dplyr", "reshape", "tidyverse", "rstudioapi",
#   "R.matlab", "mapview", "rnaturalearth", "rnaturalearthdata",
#   "abind", "pracma"
# ))

library(raster)
library(terra)
library(fields) 
library(plyr)
library(ggplot2)
library(sp)
library(sf)
library(ncdf4)
library(lubridate)
library(dplyr)
library(reshape)
library(nc)
library(ncdf4)
library(tidyverse)
library(rstudioapi)
library(R.matlab)
library(mapview)
library(rnaturalearth)
library(rnaturalearthdata)
library(abind)
require(pracma)

# ------------------------------ #
## 0.2. Set working directory -----
# ------------------------------- #
#this sets the working directory directly where your script is:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ----------------------------------------------------- #
## 0.3. Define parameters that need to be changed -----
# ---------------------------------------------------- #
#define max depth to extract data. In this case we use 150m
z<-150

#Define the desired order of months
desired_order <- month.name
months <- month.name

# Determine the start and end year at the beginning of each script
# Your data needs to have complete years (12 months)!!
start_year <- 1992
end_year <- 2022

# Define longitude and latitude ranges for the analysis
# If you don't need a subset, you still need to define the min and max of your entire dataset
# These values must be EXACTLY the SAME as in the nc. (otherwise it won't work)

# Select one or the other accordinf to the area you need
lon_range <- c(27.250, 42.000) #all black sea
lat_range <- c(40.500, 47.000) #all black sea

# --------------------------------------------- #
# 1. READ DATASETS and EXTRACT PARAMETERS ----
# --------------------------------------------- #

# --------------------------------------------- #
## 1.1. Load areas shapefile ----
# --------------------------------------------- #

# # Define and load your shapefile area
# areas <- st_read("./Data/Shapefiles/Only_Black_Sea_CLEAN_good.shp")
# plot(areas)
# class(areas)
# crs(areas)
# str(areas)
# areas
# mapview(areas) #visualize the shapefile in a world map

# --------------------------------------------- #
## 1.2. Load product data (files) ----
# --------------------------------------------- #
# --------------------------------- #
### 1.2.a. Load all the files ----
# --------------------------------- #
#Load the files
files <- dir("./Data/Pre-process_prepare_data/Chl_data", pattern=".nc")
files

setwd("./Data/Pre-process_prepare_data/Chl_data") # must set the working directory for this to work!

# ------------------------------------------ #
### 1.2.b. Extract dimensions of 1 file ----
# -------------------------------------- #
# It is essential that ALL the files have exactly the same dimensions (exact same values for lon, lat and depth)

# Open the first file to get the dimensions
name <- files[1]
aa <- nc_open(name)

#lat & lon
lat1 <- ncvar_get(aa, varid = "latitude")
str(lat1)
summary(lat1)
lon1 <- ncvar_get(aa, varid = "longitude")
str(lon1)
summary(lat1)

#times & dates
time <- ncvar_get(aa, varid = "time")
dates <- as.POSIXct(1 * time, origin = "1970-01-01", tz = "GMT") #check this in the documentation of your data

#depth & layer width
depth<-ncvar_get(aa, varid = "depth")
#Calculate the width of each depth layer
espesor_depth <-diff(depth)
#identify the depth layer for which we want to do a subset of your data (this value can be changed)
ind<-which.min(abs(depth-z)) #in this case z=150 as defined at the beginning of the script
ind
#create a subset of depth from layer 1 up to the selected depth, in this case 150m
depth2<-depth[1:ind]
depth2 #chek it
#extract the "width" of each depth layer
espesor_depth<-diff(depth)
str(espesor_depth)
#cut it to our depth of interest
espesor_depth2<-espesor_depth[1:ind]
str(espesor_depth2)
espesor_depth2

# ------------------------------------------------ #
### 1.2.c. Do the subset of the area (lon & lat) ----
# -------------------------------------------------- #
# Do a subset of lon and lat corresponding to the min and max subset we want
# We are going to overwrite them so that if you don't want to do a subset the rest of the code will run anyway ()

# Subset "lon" vector
lon <- lon1[lon1 >= lon_range[1] & lon1 <= lon_range[2]]
str(lon)
# Subset "lat" vector
lat<- lat1[lat1 >= lat_range[1] & lat1 <= lat_range[2]]
str(lat)
# Identify to which indices this vlaues correspond in the matrix
lon_indices <- which(lon1 %in% lon)
lat_indices <- which(lat1 %in% lat)

# ------------------------------------------------------- #
## 1.3. Create a unique array with all the year files ----
# -------------------------------------------------- #

# Define the dimensions of the array
num_files <- length(files)
total_time <- num_files * length(dates)  # both of these lines do the same
total_time #372

# Create an empty array with he dimensions fo all the files combined so that we can store it
chl_combined<-array(NA,dim=c(length(lon),length(lat), length(depth),total_time))
str(chl_combined)

#Create an empty vector to store the dates name if needed
dates_combined <- vector("list", total_time)

#Loop
# If you have  a lot of data, this can take some time
for (i in 1:num_files) {
  name <- files[i]
  # Print the current time being processed
  print(paste("File", i))
  aa <- nc_open(name) # Open the .nc file
  chl <- ncvar_get(aa, varid = "chl") # Extract data
  str(chl)
  # Extract dates and append to dates_combined
  time <- ncvar_get(aa, varid = "time")
  dates <- as.POSIXct(1 * time, origin = "1970-01-01", tz = "GMT")
  dates_combined[i] <- list(dates)
  print(dates_combined[i]) #you can keep track if the files are chronologically in order
  # Subset the matrix using these indices
  chl2<- chl[lon_indices, lat_indices, , ]
  chl3<-chl2[,,1:ind,]
  # Assign data for the current file to the corresponding slice in chl_combined
  if (i == 1) {
    chl_combined <- chl3
  } else {
    chl_combined <- abind::abind(chl_combined, chl3, along = 4)
  }
  
  nc_close(aa) # Close the NetCDF file
}

#check the new combined array
str(chl_combined)
class(chl_combined)
dim(chl_combined)
summary(chl_combined) # if the matrix is very big, this might take some time
image.plot(chl_combined[,,1,12]) #plot a random slice of the matrix to check it

# Accessing a specific slice of the array if needed
#slice <- chl_combined[,, ,1]

# SAVE (we recommend saving this matrix that you don't need to run it everytime you use the variable)
# Go back to the working directory where the script is so that the results are correctly saved in outputs
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#saveRDS(chl_combined, file = "./Outputs/Pre-process/Chl/monthly_chl_Black_Sea_1992_2022.RData")

# If the object was previously saved, load the saved object and assign it to a variable.
#chl_combined <- readRDS("./Outputs/Pre-process/Chl/monthly_chl_Black_Sea_1992_2022.RData")
#summary(chl_combined)
#str(chl_combined)


# -------------------------------------- #
# 2. CHLOROPHYLL DATA CALCULATIONS  ----
# -------------------------------------- #

# -------------------------------------- #
## 2.1. SURFACE DATA  ----
# -------------------------------------- #

# ----------------------------- #
### 2.1.1. Monthly matrix ----
# ---------------------------- #
#wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Create the empty array
str(chl_combined)
chl_1<-array(data=numeric(0),dim= c(length(lon),length(lat),total_time))
str(chl_1)

#loop
for (i in 1:total_time){  
  chl_1[,,i]<-chl_combined[,,1,i]
}

# Check the new object
str(chl_1)
summary(chl_1)
image.plot(chl_1[,,1])

#SAVE matrix.This matrix can be used for the Cumulative Impact analysis based on monthly anomalies
#saveRDS(chl_1, file = "./Outputs/Pre-process/Chl/surface_chl_matrix_Black_Sea_1992_2022.RData")

# --------------------------- #
### 2.1.2. Yearly matrix ----
# --------------------------- #

# Summarize monthly data to obtain annual averages
years <- seq(start_year, end_year, 1)
start_date <- seq(1, length(dates_combined), 12)
end_date <- seq(12, length(dates_combined), 12)

# Create an empty 3-dimensional array to store yearly means
surface_chl_yearly_matrix <- array(NA, dim = c(dim(chl_1)[1], dim(chl_1)[2], length(years)))

#loop
for(i in 1:length(years)){
  chl_surface_year <- chl_1[,,start_date[i]:end_date[i]]
  chl_surface_year <-apply(chl_surface_year,c(1,2),"mean")
  surface_chl_yearly_matrix[,,i] <- chl_surface_year
  print(i)
}

# Check the new object
surface_chl_yearly_matrix
str(surface_chl_yearly_matrix) # third dimension should correspond to umber of years
image.plot(surface_chl_yearly_matrix[,,4])

# calculate the yearly mean
surface_chl_mean<- apply(surface_chl_yearly_matrix,c(1,2),"mean")
image.plot(surface_chl_mean, main="Mean Surface Chl 1992-2022")

#SAVE surface data matrix. This matrix can be used for the Cumulative Impact analysis based on mean annual values.
#saveRDS(surface_chl_yearly_matrix, file = "./Outputs/Pre-process/Chl/surface_chl_matrix_YEARLY_Black_Sea_1992_2022.RData")

####  END ----