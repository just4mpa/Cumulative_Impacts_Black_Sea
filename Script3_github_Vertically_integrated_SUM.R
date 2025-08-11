# =============================================================== #
#   Unpacking cumulative impacts: exploring methodological choices for marine ecosystem analysis #
#   Supplementary material
#   Script 3- 
#    Pre-process: Prepare input data
#   Vertically integrated data weighted SUM
#   31/01/2025
# =============================================================== # 

#AIM OF THE SCRIPT

# Pre-process netcfd files to obtain Vertically integrated data: Integrated over the first 150 m of the water column
# Since the vertical resolution of the environmental variables is not homogeneous, the width of each depth layer was used as a weighting factor for vertical integration

# With this scripts we generate vertically integrated matrix per year and month (weighted SUM).
# This files will be used as input data in the different version of the Cumulative Impacts analysis (based on anual means and on monthly anomalies, respectively). 
# Here, this script is applied for Chlorophyll-a but it can be used for all the variables we want to include in the cumulative impacts analysis.


# Accomodate all the necessary parameters to your area in the first section of the script
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
# Set the working directory directly where your script is:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# ----------------------------------------------------- #
## 0.3. Define parameters that need to be changed -----
# ---------------------------------------------------- #
# Define max depth to extract data
z<-150

# Define the desired order of months
desired_order <- month.name
months <- month.name

# Determine the start and end year at the beginning of each script
# Your data needs to have complete years (12 months)!!
start_year <- 1993
end_year <- 2021

# Define longitude and latitude ranges for the analysis
# If you don't need a subset, you still need to define the min and max of your entire dataset
# These values must be EXACTLY the SAME as in the nc. (otherwise it won't work)
lon_range <- c(27.37000, 41.96259)  #all black sea
lat_range <- c(40.86000, 46.80444)  #all black_sea

# --------------------------------------------- #
# 1. READ DATASETS and EXTRACT PARAMETERS ----
# --------------------------------------------- #

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
#identify the depth layer for which we want to do a subset of your data (this value can be changed, normally from 200 m there is no light and no chl nor npp)
ind<-which.min(abs(depth-z)) #in this case z=150 as defined at the beginning of the script
#create a subset of depth from layer 1 up to the selected depth, in this case 150m
depth2<-depth[1:ind]
depth2 #chek it
#extract the "width"  of each depth layer
espesor_depth<-diff(depth)
str(espesor_depth)
#cut it to our depth of intrest
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
total_time <- num_files * length(dates)  
total_time

# Create an empty array with he dimensions fo all the files combined so that we can store it
chl_combined<-array(NA,dim=c(length(lon),length(lat), length(depth),total_time))
str(chl_combined)

#Create an empty vector to store the dates name if needed
dates_combined <- vector("list", total_time)

#Loop
for (i in 1:num_files) {
  name <- files[i]
  # Print the current time being processed
  print(paste("File", i))
  aa <- nc_open(name) # Open the .nc file
  chl <- ncvar_get(aa, varid = "chl") # Extract chl data
  str(chl)
  # Extract dates and append to dates_combined
  time <- ncvar_get(aa, varid = "time")
  dates <- as.POSIXct(1 * time, origin = "1970-01-01", tz = "GMT")
  dates_combined[i] <- list(dates)
  print(dates_combined[i]) #you can keep track if the files are chronologically in order
  # Subset the matrix using these indices
  chl2<- chl[lon_indices, lat_indices, , ]
  chl3<-chl2[,,1:ind,]
  # Assign chl data for the current file to the corresponding slice in chl_combined
  if (i == 1) {
    chl_combined <- chl3
  } else {
    chl_combined <- abind::abind(chl_combined, chl3, along = 4)
  }
  
  nc_close(aa) # Close the NetCDF file
}

# Check the new combined array
str(chl_combined)
class(chl_combined)
dim(chl_combined)
summary(chl_combined) # if the matrix is very big, this might take some time
image.plot(chl_combined[,,1,12]) #plot a random slice of the matrix to check it


# Accessing a specific slice of the array if needed
#slice <- chl_combined[,, ,1]

#SAVE (we recomend saving this so that you don't need to run it everytime you use hte variable)
#go back to the working directory where the script is so that the results are correclty saved in outputs
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#saveRDS(chl_combined, file = "./Outputs/Pre-process/Chl/monthly_chl_Black_Sea_1992_2022.RData")

 # If the object was previously saved, Load the saved object and assign it to a variable
#chl_combined <- readRDS("./Outputs/Pre-process/Chl/monthly_chl_Black_Sea_1992_2022.RData")
#summary(chl_combined)
#str(chl_combined)


# -------------------------------------- #
# 2. CHLOROPHYLL DATA CALCULATIONS  ----
# -------------------------------------- #

# -------------------------------------- #
## 2.1. VERTICALLY INTEGRATING ----
# -------------------------------------- #

# Create a 2D matrix (2D vertically integrated data)

# ----------------------------- #
### 2.1.1. Monthly matrix ----
# ---------------------------- #
#wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# With this first loop we obtain matrix of the spatially vertically integrated chl per month

# Create an empty matrix to store the outputs of the first loop
chl_matrix<-array(NA,dim=c(length(lon),length(lat),total_time))
str(chl_matrix)

#loop
# weighted SUM vertical integration
for (i in 1:total_time){  
  print(paste("Month",i))
  for (j in 1:length(lon)){
    for (k in 1:length(lat)){
      chl_ <- chl_combined[j,k,,i] #access the chla value of lon (j), lat(k) and month(i)
      chl_m<-sum(chl_*espesor_depth2, na.rm=TRUE)#multiply each cell (pixel) by the width of the layer, and sum all of them together
      chl_matrix[j,k,i]<-chl_m
    }
  }
}

# Check the new object
str(chl_matrix) 
summary(chl_matrix)
image.plot(chl_matrix[,,1])

#SAVE matrix. This matrix can be used for the Cumulative Impact analysis based on monthly anomalies.
#saveRDS(chl_matrix, file = "./Outputs/Pre-process/Chl/v_int_SUM_chl_Black_sea_1992_2022.RData")


# --------------------------- #
### 2.1.2. Yearly matrix ----
# --------------------------- #

# Summarize monthly data to obtain annual averages
years <- seq(start_year, end_year, 1)
start_date <- seq(1, length(dates_combined), 12)
end_date <- seq(12, length(dates_combined), 12)

# Create an empty 3-dimensional array to store yearly means
v_int_chl_yearly_matrix <- array(NA, dim = c(dim(chl_matrix)[1], dim(chl_matrix)[2], length(years)))

#loop 
for(i in 1:length(years)){
  v_int_chl_yearly <- chl_matrix[,,start_date[i]:end_date[i]]
  v_int_chl_yearly <-apply(v_int_chl_yearly,c(1,2),"mean")
  v_int_chl_yearly_matrix[,,i] <- v_int_chl_yearly
  print(i)
}

# Check the new object
v_int_chl_yearly_matrix 
str(v_int_chl_yearly_matrix)
summary(v_int_chl_yearly_matrix)
image.plot(v_int_chl_yearly_matrix[,,4]) 

# Calculate the yearly mean
chl_mean <- apply(v_int_chl_yearly_matrix,c(1,2),"mean")
image.plot(chl_mean, main="Mean v_int SUM Chl 1992-2022")

#SAVE matrix. his matrix can be used for the Cumulative Impact analysis based on mean annual values
#saveRDS(v_int_chl_yearly_matrix, file = "./Outputs/Pre-process/Chl/v_int_SUM_chl_matrix_YEARLY_Black_Sea_1992_2022.RData")

####  END ----