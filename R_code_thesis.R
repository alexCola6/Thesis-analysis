install.packages("plm")
install.packages("dplyr")
# library(plm)
library(dplyr)
library(ncdf4)
library(ggplot2)
library(tidyr)

############################### Merging datasets ############################### 

rainfall_data <- read.csv('/Users/alexandrecolassis/Documents/0_Documents MacBook/Uni/UNIL/UNIL Master 2022-2023/thesis /_DATA /rainfall_monthly_2013-2020.csv')
temperature_data <- read.csv('/Users/alexandrecolassis/Documents/0_Documents MacBook/Uni/UNIL/UNIL Master 2022-2023/thesis /_DATA /monthly_mean_surface_temperature.csv')

# Set 'time' columns as Date format
rainfall_data$time <- as.Date(rainfall_data$time)
temperature_data$time <- as.Date(temperature_data$time)

# Apply the conversion formula to surface_temperature to get values in celsius 
temperature_data$surface_temperature <- (temperature_data$surface_temperature * 0.00341802) - 124.15

# Calculate the mean surface_temperature for each timeframe
temperature_data <- temperature_data %>%
  group_by(time) %>%                              # Group by time
  summarize(surface_temperature = mean(surface_temperature, na.rm = TRUE),  # Calculate mean
            .groups = "drop")                     # Drop grouping after summarization

# Extract year and month for alignment
rainfall_data <- rainfall_data %>%
  mutate(year = format(time, "%Y"),
         month = format(time, "%m"))

temperature_data <- temperature_data %>%
  mutate(year = format(time, "%Y"),
         month = format(time, "%m"))

# Merge datasets by "time" and "spatial reference"
merged_data <- left_join(temperature_data, rainfall_data, 
                         by = c("year", "month"))

# Select only the columns needed
final_data <- merged_data %>% 
  select(time.x, `surface_temperature`, rainfall) %>%
  rename(time = time.x)

# Format the 'time' column to keep only year and month
final_data <- final_data %>%
  mutate(time = format(as.Date(time), "%Y-%m"))



##################### Trying to obtain BA from files (loop) #####################

# Define function to extract burned area and date from .nc file
extract_burned_area <- function(file_path) {
  # Open the NetCDF file
  nc <- nc_open(file_path)
  
  # Extract burned area data (replace "burned_area_variable" with the actual variable name)
  burned_area <- ncvar_get(nc, "burned_area_variable") 
  
  # Extract the date (replace "time_variable" with the actual time variable name)
  time <- ncvar_get(nc, "time_variable")
  
  # Convert the time to a proper Date format if needed
  # (adjust the date conversion depending on how the time is stored in the .nc files)
  date <- as.Date("1970-01-01") + time  # Replace as needed
  
  # Close the NetCDF file
  nc_close(nc)
  
  # Return data as a data frame
  data.frame(date = date, burned_area = burned_area)
}


# Define the folder where your NetCDF files are stored
folder_path <- '/Users/alexandrecolassis/Documents/0_Documents MacBook/Uni/UNIL/UNIL Master 2022-2023/thesis /_DATA /BA/BA_2013-2020'

# List all the NetCDF files (e.g., all files for 2020)
files <- list.files(folder_path, pattern = "^BA.*\\.nc$", full.names = TRUE)

# Define the target coordinates and time for extraction
target_lat <- 16.625
target_lon <- -14.875
target_time <- as.Date("2020-12-31")

# Initialize a list to store the burned area values for each file
monthly_burned_area <- list()

# Loop through each file
for (file in files) {
  # Open the NetCDF file
  nc_data <- nc_open(file)
  
  # Extract dimension variables (lat, lon, time)
  lat_vals <- ncvar_get(nc_data, "lat")
  lon_vals <- ncvar_get(nc_data, "lon")
  time_vals <- ncvar_get(nc_data, "time")
  
  # Find the closest latitude and longitude indices
  lat_idx <- which.min(abs(lat_vals - target_lat))
  lon_idx <- which.min(abs(lon_vals - target_lon))
  
  # Get the time reference and convert target_time to the number of days since
  time_units <- ncatt_get(nc_data, "time", "units")$value
  ref_date <- as.Date(sub("days since ", "", time_units))
  target_time_days <- as.numeric(difftime(target_time, ref_date, units = "days"))
  
  # Find the closest time index
  time_idx <- which.min(abs(time_vals - target_time_days))
  
  # Corrected start and count vectors for 3 dimensions (lon, lat, time)
  start <- c(lon_idx, lat_idx, time_idx)
  count <- c(1, 1, 1)
  
  # Extract the total burned area
  total_val <- ncvar_get(nc_data, "Total", start = start, count = count)
  
  # Store the extracted value for this file
  monthly_burned_area[[file]] <- total_val
  
  # Close the NetCDF file
  nc_close(nc_data)
  
  # Print the extracted value for this file
  print(paste("Burned area for file", file, ":", total_val))
}

# Print the list of total burned areas for each file
print(monthly_burned_area)


############# Create a dataset with the extracted BA from the loop #############

# Define the folder where your NetCDF files are stored
folder_path <- '/Users/alexandrecolassis/Documents/0_Documents MacBook/Uni/UNIL/UNIL Master 2022-2023/thesis /_DATA /BA/BA_2013-2020'

# List all NetCDF files that start with "BA" and end in ".nc"
files <- list.files(folder_path, pattern = "^BA.*\\.nc$", full.names = TRUE)

# Initialize an empty data frame to store results
burned_area_data <- data.frame(Time = character(), Burned_Area = numeric(), stringsAsFactors = FALSE)

# Define the target coordinates (lat, lon) for extraction
target_lat <- 16.625
target_lon <- -14.875

# Loop through each file
for (file in files) {
  # Extract year and month from the filename (assuming format like BA202012.nc)
  file_name <- basename(file)
  year_month <- substr(file_name, 3, 8)  # Extracts "202012" from "BA202012.nc"
  year <- substr(year_month, 1, 4)
  month <- substr(year_month, 5, 6)
  
  # Open the NetCDF file
  nc_data <- nc_open(file)
  
  # Extract dimension variables
  lat_vals <- ncvar_get(nc_data, "lat")
  lon_vals <- ncvar_get(nc_data, "lon")
  time_vals <- ncvar_get(nc_data, "time")
  
  # Find the closest indices for target coordinates
  lat_idx <- which.min(abs(lat_vals - target_lat))
  lon_idx <- which.min(abs(lon_vals - target_lon))
  
  # Extract the time reference and convert target time to the number of days since
  time_units <- ncatt_get(nc_data, "time", "units")$value
  ref_date <- as.Date(sub("days since ", "", time_units))
  target_time <- as.Date(paste0(year, "-", month, "-01"))
  target_time_days <- as.numeric(difftime(target_time, ref_date, units = "days"))
  time_idx <- which.min(abs(time_vals - target_time_days))
  
  # Extract the burned area value (assuming the variable name is "Total")
  total_val <- ncvar_get(nc_data, "Total", start = c(lon_idx, lat_idx, time_idx), count = c(1, 1, 1))
  
  # Close the NetCDF file
  nc_close(nc_data)
  
  # Store the result in the data frame
  burned_area_data <- rbind(burned_area_data, data.frame(Time = paste(year, month, sep = "-"), Burned_Area = total_val))
}

# Print the data frame
print(burned_area_data)


####################### Different plot options for the BA #######################

# Time Series Line Plot
ggplot(burned_area_data, aes(x = as.Date(paste0(Time, "-01")), y = Burned_Area)) +
  geom_line(color = "firebrick") +
  labs(title = "Monthly Burned Area Over Time", x = "Time", y = "Burned Area (km²)") +
  theme_minimal()

# Bar Plot of Monthly Burned Areas
ggplot(burned_area_data, aes(x = as.Date(paste0(Time, "-01")), y = Burned_Area)) +
  geom_bar(stat = "identity", fill = "orange") +
  labs(title = "Monthly Burned Area", x = "Time", y = "Burned Area (km²)") +
  theme_minimal()

# Seasonal Pattern Visualization (Box Plot by Month)
# Extract month as a separate variable
burned_area_data$Month <- format(as.Date(paste0(burned_area_data$Time, "-01")), "%B")

ggplot(burned_area_data, aes(x = Month, y = Burned_Area)) +
  geom_boxplot(fill = "darkorange", color = "black") +
  labs(title = "Monthly Burned Area Distribution by Month", x = "Month", y = "Burned Area (km²)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Heatmap for Monthly Burned Area Over Years 
# Extract Year and Month as separate variables
burned_area_data$Year <- format(as.Date(paste0(burned_area_data$Time, "-01")), "%Y")
burned_area_data$Month <- format(as.Date(paste0(burned_area_data$Time, "-01")), "%m")

ggplot(burned_area_data, aes(x = Month, y = Year, fill = Burned_Area)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightyellow", high = "darkred", name = "Burned Area (km²)") +
  labs(title = "Monthly Burned Area Heatmap", x = "Month", y = "Year") +
  theme_minimal() +
  theme(
    plot.background = element_rect(color = "black", fill = NA, size = 1), # Add a frame around the entire plot
    plot.margin = margin(10, 10, 10, 10) # Add spacing around the plot for neatness
  )

# Heatmap for Surface Temperature Over Years 
# Extract Year and Month as separate variables
temperature_data$year <- format(as.Date(paste0(temperature_data$time, "-01")), "%Y")
temperature_data$month <- format(as.Date(paste0(temperature_data$time, "-01")), "%m")

ggplot(temperature_data, aes(x = month, y = year, fill = surface_temperature)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightyellow", high = "darkred", name = "Surface Temperature (°C)") +
  labs(title = "Monthly Surface Temperature Heatmap", x = "Month", y = "Year") +
  theme_minimal() +
  theme(
    plot.background = element_rect(color = "black", fill = NA, size = 1), # Add a frame around the entire plot
    plot.margin = margin(10, 10, 10, 10) # Add spacing around the plot for neatness
  )

# Heatmap for Rainfall Over Years 
# Extract Year and Month as separate variables
rainfall_data$year <- format(as.Date(paste0(rainfall_data$time, "-01")), "%Y")
rainfall_data$month <- format(as.Date(paste0(rainfall_data$time, "-01")), "%m")

ggplot(rainfall_data, aes(x = month, y = year, fill = rainfall)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightyellow", high = "darkred", name = "Rainfall (mm)") +
  labs(title = "Monthly Rainfall Heatmap", x = "Month", y = "Year") +
  theme_minimal() +
  theme(
    plot.background = element_rect(color = "black", fill = NA, size = 1), # Add a frame around the entire plot
    plot.margin = margin(10, 10, 10, 10) # Add spacing around the plot for neatness
  )


########## make sure that my columns have the same name before merging ##########
names(final_data)[names(final_data) == "time"] <- "Time"

# Select only the Time and Burned_Area columns from burned_area_data
burned_area_data <- burned_area_data %>% select(Time, Burned_Area)

# Merge the two datasets using full_join to keep the largest time frame
merged_data <- full_join(final_data, burned_area_data, by = "Time")

# Arrange by date to ensure chronological order
merged_data <- merged_data %>% arrange(Time)

# Print the first few rows to verify
print(head(merged_data))
print(tail(merged_data))


############# Remove surface_temperature outliers from merged_data #############
library(dplyr)

# Define the lower and upper quantiles
lower_bound <- quantile(merged_data$surface_temperature, 0.10, na.rm = TRUE)  # 10th percentile
upper_bound <- quantile(merged_data$surface_temperature, 0.90, na.rm = TRUE)  # 90th percentile

# Filter the dataset to exclude outliers
merged_data <- merged_data %>%
  filter(surface_temperature >= lower_bound, surface_temperature <= upper_bound)

# Check the result
summary(merged_data$surface_temperature)



############################### Time series model #################################

model_time_series <- lm(
  Burned_Area ~ surface_temperature + rainfall + surface_temperature:rainfall,
  data = merged_data
)
summary(model_time_series)

##### Test with non-linear model, but does not improve the results
model_time_series_non_linear <- lm(
  Burned_Area ~ surface_temperature + rainfall + I(rainfall^2) + surface_temperature:rainfall,
  data = merged_data
)
summary(model_time_series_non_linear)


############################## Model including lag ############################## 

# Create lagged variable for Burned_Area
merged_data <- merged_data %>%
  mutate(Burned_Area_Lag = lag(Burned_Area))

# Fit model with lagged Burned_Area
model_time_series_lag <- lm(
  Burned_Area ~ surface_temperature + rainfall + surface_temperature:rainfall + Burned_Area_Lag,
  data = merged_data
)

summary(model_time_series_lag)


##################### Time-Series Cross-Validation setting ##################### 
######################## Rolling window Cross-Validation ########################

install.packages("caret")
library(caret)

merged_data <- merged_data %>% filter(!is.na(Burned_Area_Lag))

# Set up time series cross-validation
train_control_ts <- trainControl(
  method = "timeslice",  # Time series-specific cross-validation
  initialWindow = 48,   # Initial training window
  horizon = 12,          # Forecast horizon
  fixedWindow = TRUE     # Use fixed window (not rolling)
)

# Train the time series model with cross-validation
model_time_series_cv <- train(
  Burned_Area ~ surface_temperature + rainfall + surface_temperature:rainfall + Burned_Area_Lag,
  data = merged_data,
  method = "lm",
  trControl = train_control_ts
)

print(model_time_series_cv)


################################################################################
################################################################################
############################# Projection scenarios #############################

library(dplyr)
library(tidyr)
library(lubridate)

# Convert `Time` to Year-Month format
merged_data <- merged_data %>%
  mutate(Time = ym(Time))  # 'ym()' interprets "YYYY-MM"

# Step 1: Compute monthly mean surface_temperature and rainfall for 2013-2020
monthly_means <- merged_data %>%
  filter(year(Time) >= 2013 & year(Time) <= 2020) %>%
  mutate(Month = month(Time)) %>%
  group_by(Month) %>%
  summarise(
    mean_surface_temperature = mean(surface_temperature, na.rm = TRUE),
    mean_rainfall = mean(rainfall, na.rm = TRUE)
  )

# Step 2: Create future scenarios for 2080
future_scenarios <- monthly_means %>%
  mutate(
    Scenario1_surface_temperature = mean_surface_temperature + 1.7,  # Temp increase by 1.7°C
    Scenario1_rainfall = mean_rainfall,                              # Rainfall constant
    
    Scenario2_surface_temperature = mean_surface_temperature,       # Temp constant
    Scenario2_rainfall = mean_rainfall * 1.05,                      # Rainfall increases by 5%
    
    Scenario3_surface_temperature = mean_surface_temperature + 1.7, # Temp increase by 1.7°C
    Scenario3_rainfall = mean_rainfall * 1.05,                      # Rainfall increases by 5%
    
    Scenario4_surface_temperature = mean_surface_temperature,       # Temp constant
    Scenario4_rainfall = mean_rainfall                              # Rainfall constant
  )

# Step 3: Reshape data for projections
future_data <- future_scenarios %>%
  pivot_longer(cols = starts_with("Scenario"), 
               names_to = c(".value", "Scenario"),
               names_sep = "_") %>%
  mutate(
    Time = as.Date(sprintf("2080-%02d-01", Month)),
    Burned_Area_Lag = 0,  # Assuming no prior burned area for projections
  ) 


# Create separate datasets for each scenario
scenario1 <- future_data %>%
  select(Time, Month, surface_temperature = Scenario1, rainfall = mean_rainfall, Burned_Area_Lag) %>%
  distinct(Time, .keep_all = TRUE)

scenario2 <- future_data %>%
  select(Time, Month, surface_temperature = mean_surface_temperature, rainfall = Scenario2, Burned_Area_Lag) %>%
  group_by(Month) %>%  # Group by month
  slice(2) %>%         # Keep the second observation for each month
  ungroup()            # Ungroup to finish processing

scenario3 <- future_data %>%
  filter(Scenario == "surface") %>%
  select(Time, Month, surface_temperature = Scenario3, Burned_Area_Lag) %>%
  left_join(
    future_data %>%
      filter(Scenario == "rainfall") %>%
      select(Time, rainfall = Scenario3),
    by = "Time"
  )
scenario3 <- scenario3 %>%
  select(Time, Month, surface_temperature, rainfall, Burned_Area_Lag)

scenario4 <- future_data %>%
  select(Time, Month, surface_temperature = Scenario4, rainfall = mean_rainfall, Burned_Area_Lag) %>%
  distinct(Time, .keep_all = TRUE)



library(plm)

# Convert scenario1 to pdata.frame
scenario1_pdata <- pdata.frame(scenario1, index = c("Time"))

# Predict burned area for scenario1
scenario1 <- scenario1 %>%
  mutate(Predicted_Burned_Area = predict(model_time_series_lag, newdata = scenario1_pdata))

# Repeat for scenario2
scenario2_pdata <- pdata.frame(scenario2, index = c("Time"))
scenario2 <- scenario2 %>%
  mutate(Predicted_Burned_Area = predict(model_time_series_lag, newdata = scenario2_pdata))

# Repeat for scenario3
scenario3_pdata <- pdata.frame(scenario3, index = c("Time"))
scenario3 <- scenario3 %>%
  mutate(Predicted_Burned_Area = predict(model_time_series_lag, newdata = scenario3_pdata))

# Repeat for scenario4
scenario4_pdata <- pdata.frame(scenario4, index = c("Time"))
scenario4 <- scenario4 %>%
  mutate(Predicted_Burned_Area = predict(model_time_series_lag, newdata = scenario4_pdata))


projections <- bind_rows(
  scenario1 %>% mutate(Scenario = "Scenario 1"),
  scenario2 %>% mutate(Scenario = "Scenario 2"),
  scenario3 %>% mutate(Scenario = "Scenario 3"),
  scenario4 %>% mutate(Scenario = "Scenario 4")
)

