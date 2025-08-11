# Cumulative_Impacts_Black_Sea
Cumulative Impacts for the Black Sea methodological comparison

# Introduction
In this repository you can find data and code to replicate the comparison of cumulative impact approaches, using the Black Sea as a case study. We explore and compare two methodological approaches for highly resolved, spatially explicit environmental cumulative impact analysis. We combined two temporal resolutions (annual vs monthly) and two types of metrics (absolute means vs percent-based of anomalies) resulting in what we would refer to as two approaches: annual mean absolute values vs. percentages of monthly anomalies. Additionally, we employed three different considerations for the vertical component of the environmental input data: surface data, vertically integrated weighted mean, and vertically integrated weighted sum) and together, these factors produce six different scenarios for cumulative impact analysis.
Scripts 1 to 5 are used to prepare the input data. The files generated with these scripts will be used as input data in the different versions of the Cumulative Impacts analysis (based on annual mean absolute values and percentages of monthly anomalies, respectively). In these examples, the scripts are applied for Chlorophyll-a and Temperature, as an example, but they can be used for any variable of interest that will be included in the cumulative impacts analysis.
Scripts 6 to 8 provide the code to perform different cumulative impact approaches. In these scripts, we use SURFACE environmental variables as an example, but any of the previous pre-process outputs can be used instead (vertical integration weighted mean, vertical integration weighted sum, vertical integration 1m interpolation). All this processed data can be found in the repository under “Data” and “CI_annual_mean_absolute_values” or “CI_percentages_monthly_anomalies”, depending on the approach of interest.

# Must consider
-Keep in mind that data needs to have complete years (12 months) for the script to work properly.
-Latitude and longitude information of your study area are required.
- For scripts 1 to 5, data must be organized as a folder containing .nc files of the variable of interest per year with a monthly resolution. (If you are going to run the code for several variables, have one folder per variable).



# Packages required
install.packages(c ("raster", "terra", "fields", "plyr", "ggplot2", "sp", "sf", "ncdf4", "lubridate", "dplyr", "reshape", "tidyverse", "rstudioapi", "R.matlab", "mapview", "rnaturalearth", "rnaturalearthdata", "abind", "pracma", "circular", "CircStats", "arrayhelpers", "gstat", "reshape2", "nc", "plotly"))


# Scripts index
## 1.	Pre-process-scripts to prepare the input data
Script 1: Surface data (only for Chl-a as e.g.).
Script 2: Vertically integrated weighted MEAN (only for Temperature as e.g.).
Script 3: Vertically integrated weighted SUM (only for Chl-a as e.g.).
Script 4: Alternative vertically integrated with 1 m interpolations MEAN (Temp).
Script 5: Alternative vertically integrated with 1 m interpolations SUM (Chl-a).
*We always need to extract the data in two forms, as an annual mean matrix and as monthly mean matrix to have the outputs necessary for the two types of CI (annual mean absolute values and percentages of monthly anomalies, respectively).
*In the sample script we only process data for Chlrorphyll-a (script 1, 3 and 5) and for Temperature (script 2 and 4). But this needs to be done with all the data sets that need to be included in the subsequent CI analysis).

## 2.	Cumulative impacts analysis scripts
Script 6: Based on annual mean absolute values.
Script 7: Based on percentages of monthly anomalies.
Script 8: Alternative- Based on annual mean absolute values but with non-significant p- values set to 0.
*For all the CI analysis we will only use surface data of all the variables a simplification of the scripts. But the same can be applied to any of the other types of data generated with scripts 1 to 5).

# Additional information
We also provide a folder with the final figures of the six different scenarios for cumulative impact analysis compared in the paper, plus the two alternative trials (CI based on annual mean absolute values but with non-significant p- values set to 0 and using the 1 m interpolation for the scenario of vertically integrated SUM). 
This folder contains 8 other folders with the corresponding figures:
1.	CI_annual_mean_absolute_values_Surface
2.	CI_annual_mean_absolute_values_Vertically_int_MEAN
3.	CI_annual_mean_absolute_values_Vertically_int_SUM
4.	CI_percentages_monthly_anomalies_Surface
5.	-CI_percentages_monthly_anomalies_Vertically_int_MEAN
6.	CI_percentages_monthly_anomalies_Vertically_int_SUM
7.	Alternative_CI_annual_mean_absolute_values_Vertically_int_SUM_1m_interpolation
8.	Alternative_CI_annual_mean_absolute_values_Vertically_int_SUM_1m_interpolation

