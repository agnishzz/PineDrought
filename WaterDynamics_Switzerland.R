##Analysing water dynamics for Switzerland. This script will be used for the future simulations


#Load the libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(corrplot)
library(RSQLite)
library(here)
library(cowplot)
library(patchwork)
library(xlsx)
library(readxl)
library(writexl)
library(grid)


# Define the directory containing the SQLite files
data_dir <- "D:/Outputs_Switzerland/"
# Get a list of all SQLite files in the directory
sqlite_files <- list.files(data_dir, pattern = "\\.sqlite$", full.names = TRUE)

# Initialize an empty list to store the data
data_list <- list()

# Loop through each file, read the data, and store it in the list
for (file_path in sqlite_files) {
  # Extract the base name of the file (without the directory path)
  file_name <- basename(file_path)
  # Create a connection to the SQLite database
  db_connection <- dbConnect(RSQLite::SQLite(), dbname = file_path)
  water_data <- dbReadTable(db_connection, "water")
  # Store the data in the list with appropriate names
  data_list[[paste0("water_", file_name)]] <- water_data
  dbDisconnect(db_connection)
}
# Print the names of the data frames stored in the list to verify
print(names(data_list))

####Add a simulation site column, similarly climate model column and rcp column will be added too.
# Function to determine Simulation_site based on table name
get_simsite <- function(table_name) {
  if (grepl("_Visp_", table_name)) {
    return("CH_Visp")
  } else {
    return(NA)  # Handle unexpected cases if any
  }
}

# Add Simulation_site column to each table in cumulative_volume_list
for (table_name in names(data_list)) {
  simsite <- get_simsite(table_name)
  data_list[[table_name]]$Simulation_site <- simsite
}

# Print tables with added Simulation_site column
for (table_name in names(data_list)) {
  cat("Table:", table_name, "\n")
  print(data_list[[table_name]])
  cat("\n")
}


##
get_cm <- function(table_name) {
  if (grepl("_CC1_", table_name)) {
    return("CC1")
  } else if (grepl("_CC2_", table_name)) {
    return("CC2")
  } else if (grepl("_CC3_", table_name)) {
    return("CC3")
  }else if (grepl("_CC4_", table_name)) {
    return("CC4")
  }else if (grepl("_CC5_", table_name)) {
    return("CC5")
  }else if (grepl("_Obs_", table_name)) {
    return("Reference")
  }
  else {
    return(NA)  # Handle unexpected cases if any
  }
}

# Add Simulation_site column to each table in cumulative_volume_list
for (table_name in names(data_list)) {
  cm <- get_cm(table_name)
  data_list[[table_name]]$Climate_model <- cm
}

# Print tables with added Simulation_site column
for (table_name in names(data_list)) {
  cat("Table:", table_name, "\n")
  print(data_list[[table_name]])
  cat("\n")
}

##
get_rcp <- function(table_name) {
  if (grepl("_rcp26.sqlite", table_name)) {
    return("RCP 2.6")
  } else if (grepl("_rcp45.sqlite", table_name)) {
    return("RCP 4.5")
  } else if (grepl("_rcp85.sqlite", table_name)) {
    return("RCP 8.5")
  } else if (grepl("_ref.sqlite", table_name)) {
    return("Refclim")
  }else {
    return(NA)  # Handle unexpected cases if any
  }
}

# Add Simulation_site column to each table in cumulative_volume_list
for (table_name in names(data_list)) {
  rcp <- get_rcp(table_name)
  data_list[[table_name]]$Scenario <- rcp
}

# Print tables with added Simulation_site column
for (table_name in names(data_list)) {
  cat("Table:", table_name, "\n")
  print(data_list[[table_name]])
  cat("\n")
}

print(names(data_list))


##Merge to get entire water table
water<- do.call(rbind, lapply(names(data_list), function(name) {
  if (grepl("^water_", name)) {
    return(data_list[[name]])
  } else {
    return(NULL)
  }
}))
print(water)
head(water)
str(water)
tail(water)
water$Year <- water$year + 1997
water$Country <- "Switzerland"
head(water)

writexl::write_xlsx(water, "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/WaterDynamics_Switzerland.xlsx")
