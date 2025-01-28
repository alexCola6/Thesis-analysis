# Thesis analysis

The intensification of extreme bushfires around the world is a major issue, especially from an economic perspective. In order to estimate the evolution of bushfires due to climate change, an analysis of the effects of two fire weather components, surface temperature and rainfall, on burned areas is performed using a linear regression model. The results suggest that an increase in surface temperature has a significant effect on burned areas, but no conclusions can be reached regarding the impact of rainfall. To account for climate change, projections for the year 2080 are made following the IPCC intermediate scenario. The results indicate that bushfires are expected to be more frequent when accounting for surface temperature increases and constant rainfall, as well as for increases in both fire weather variables. On the other hand, an increase in rainfall alone is associated with a reduction in predicted burned areas. Overall, climate change is expected to exacerbate the occurrence of bushfires and they are set to become more intense.

# Data describtion
All data used expand from March 2013 to December 2020, and are retrieved from different sources.

Burned area data come from the Global Fire Emissions Database (GFED5) and can be downloaded from the following website : https://zenodo.org/records/7668424. 

Surface temperature data are collected from Digital Earth Africa, based on the United States Geological Survey’s (USGS) Landsat satellite program. A sandbox must be used to download the required data and the code is described in the "csv_export.ipynb" file.  

Rainfall data are also recovered from Digital Earth Africa, as it provides an open access copy of the monthly Climate Hazards Group InfraRed Precipitation with Station data (CHIRPS). Again, a sandbox must be used to download the required data and the code is described in the "test_rainfall.ipynb" file. 

Here is the link to access the Digital Earth Africa sandbox: https://docs.digitalearthafrica.org/en/latest/data_specs/CHIRPS_specs.html#Amazon-Web-Service-S3.

## Files 
- "R_code_thesis.R": this file contains the R script for the data analysis
- "csv_export.ipynb": this file contains the steps to retrieve and download the surface temperature data from Digital Earth Africa
- "monthly_mean_surface_temperature.csv": this file contains the data for monthly mean surface temperature following the download from Digital Earth Africa
- "rainfall_monthly_2013-2020.csv": this file contains the data for monthly mean rainfall following the download from Digital Earth Africa
- "test_rainfall.ipynb": this file contains the steps to retrieve and download the rainfall data from Digital Earth Africa

# Data treatment and analysis
Once all variables are recovered, the analysis is performed using RStudio. The description of the code used is detailed in the "R_code_thesis.R" file. 

In order to replicate the analysis, data have to be set up in R and organised uniformally between March 2013 and December 2020 with the same location. In this project, the selected coordinates were: latitude 16.625 and longitude -14.875. Once the files are available, the datasets are organised in columns, with: Time, "Variable of interest". Then, a merged_data dataset is created to have all variables together. It contains the following columns: Time, surface_temperature, rainfall, Burned_Area, Burned_Area_Lag. 

This allows to perform the linear regression model with this full dataset. When running the code for the model, results are provided in a table, allowing to analyse the value of the different coefficients and seeing whether they are significant or not. 

Since the results have turned out to be non-statistically significant, a cross-validation test is performed, more precisely a fixed window cross-validation test, to evaluate the validity and quality of the linear model. Nevertheless, the results are not significant. 

Projections of future burned areas due to climate change are performed using an increase in surface temperature of +1.7°C for 2080, compared to the period 2013-2020, and an increase in the amount of precipitations of +5% for 2080, compared to 2013-2020. The projections follow four scenarios: the first one with an increase in surface temperature while keeping rainfall constant; the second scenario with constant surface temperature and an increase in rainfall; the third scenario with an increase in both parameters; and the final scenario with both variables kept constant. 
The results indicate that burned areas are expected to increase by 2080, especially when considering an increase of 1.7°C in surface temperature. This result is in line with the estimated coefficient of the linear model. Nevertheless, other estimations contradict the coefficient obtained in the linear regression model. 

Overall, the results should be interpreted with caution since there a generally non-statistically significant. They can serve as an indicator, but the model would would benefit from more data and explanatory variables. 




