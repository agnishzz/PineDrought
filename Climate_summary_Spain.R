#Climate summaries, Spanish plots
#REVISED 7th October 2026:
#  - evapotranspiration is now FAO-56 Penman-Monteith (was temperature-based)
#  - CWB follows from the new PET
#  - Eto added to every period
#  - Growing Period (April-September, suffix _GP) added
#  - output directory changed

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
library(readxl)
library(writexl)
library(grid)
library(lubridate)
library(stringr)

# Define the directory containing the SQLite files
data_dir <- "D:/Climate/Spain"

# NEW: output directory
out_dir  <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Get a list of all SQLite files in the directory
sqlite_files <- list.files(data_dir, pattern = "\\.sqlite$", full.names = TRUE)

# Initialize an empty list to store the data
data_list <- list()

# Loop through each file, read the data, and store it in the list
for (file_path in sqlite_files) {
  file_name <- basename(file_path)
  db_connection <- dbConnect(RSQLite::SQLite(), dbname = file_path)
  S0_data <- dbReadTable(db_connection, "S0")
  S1_data <- dbReadTable(db_connection, "S1")
  S2_data <- dbReadTable(db_connection, "S2")
  S3_data <- dbReadTable(db_connection, "S3")
  S4_data <- dbReadTable(db_connection, "S4")
  data_list[[paste0("S0_", file_name)]] <- S0_data
  data_list[[paste0("S1_", file_name)]] <- S1_data
  data_list[[paste0("S2_", file_name)]] <- S2_data
  data_list[[paste0("S3_", file_name)]] <- S3_data
  data_list[[paste0("S4_", file_name)]] <- S4_data
  dbDisconnect(db_connection)
}
print(names(data_list))


#### Simulation site, climate model and rcp columns --------------------------
# Patterns anchored with ^ so they can only match the prefix we added above.
get_simsite <- function(table_name) {
  if (grepl("^S0_", table_name)) {
    return("SP_Prades")
  } else if (grepl("^S1_", table_name)) {
    return("SP_Arcalis1")
  } else if (grepl("^S2_", table_name)) {
    return("SP_Arcalis2")
  } else if (grepl("^S3_", table_name)) {
    return("SP_Arcalis3")
  } else if (grepl("^S4_", table_name)) {
    return("SP_Arcalis4")
  } else {
    return(NA)
  }
}

for (table_name in names(data_list)) {
  data_list[[table_name]]$Simulation_site <- get_simsite(table_name)
}

##
get_cm <- function(table_name) {
  if (grepl("_CC1_", table_name)) {
    return("CC1")
  } else if (grepl("_CC2_", table_name)) {
    return("CC2")
  } else if (grepl("_CC3_", table_name)) {
    return("CC3")
  } else if (grepl("_CC4_", table_name)) {
    return("CC4")
  } else if (grepl("_CC5_", table_name)) {
    return("CC5")
  } else if (grepl("_Obs_", table_name)) {
    return("Reference")
  } else {
    return(NA)
  }
}

for (table_name in names(data_list)) {
  data_list[[table_name]]$Climate_model <- get_cm(table_name)
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
  } else {
    return(NA)
  }
}

for (table_name in names(data_list)) {
  data_list[[table_name]]$Scenario <- get_rcp(table_name)
}

print(names(data_list))

##### Merge all the data into one single table for analysis ------------------
merged_data <- bind_rows(
  lapply(names(data_list), function(name) data_list[[name]])
)

head(merged_data)
tail(merged_data)

stopifnot(!anyNA(merged_data$Simulation_site),
          !anyNA(merged_data$Climate_model),
          !anyNA(merged_data$Scenario))

###############################################################################
### CHANGED BLOCK: FAO-56 Penman-Monteith reference evapotranspiration and CWB
###############################################################################

# Site coordinates (decimal degrees) and elevation (m a.s.l.)
site_meta <- tibble::tribble(
  ~Simulation_site, ~lat,        ~lon,         ~elev,
  "SP_Arcalis1",    42.36254117, 1.095519333,  660.3333333,
  "SP_Arcalis2",    42.34864743, 1.088869286,  987,
  "SP_Arcalis3",    42.3673052,  1.049716,     1216,
  "SP_Arcalis4",    42.34418886, 1.097959857,  1321.571429,
  "SP_Prades",      41.33691,    1.01,         928.3
)
stopifnot(all(unique(merged_data$Simulation_site) %in% site_meta$Simulation_site))

# Saturation vapour pressure (kPa) at temperature T (deg C)
e0 <- function(T) 0.6108 * exp(17.27 * T / (T + 237.3))

# Extraterrestrial radiation (MJ m-2 d-1) from latitude and day of year
ra_fao <- function(lat_deg, doy) {
  phi   <- lat_deg * pi / 180
  dr    <- 1 + 0.033 * cos(2 * pi * doy / 365)
  delta <- 0.409 * sin(2 * pi * doy / 365 - 1.39)
  ws    <- acos(pmin(pmax(-tan(phi) * tan(delta), -1), 1))
  (24 * 60 / pi) * 0.0820 * dr *
    (ws * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(ws))
}

u2 <- 2   # m/s, FAO-56 default where wind speed is not measured

merged_data <- merged_data %>%
  left_join(site_meta, by = "Simulation_site") %>%
  mutate(
    date  = make_date(year, month, day),
    doy   = yday(date),
    Tmean = (min_temp + max_temp) / 2,
    
    # psychrometrics
    P_atm  = 101.3 * ((293 - 0.0065 * elev) / 293)^5.26,   # kPa
    gamma  = 0.000665 * P_atm,
    Delta  = 4098 * e0(Tmean) / (Tmean + 237.3)^2,
    es     = (e0(max_temp) + e0(min_temp)) / 2,
    ea     = pmax(es - vpd, 0.01),        # actual vapour pressure, from VPD
    
    # radiation
    Ra     = ra_fao(lat, doy),
    Rso    = (0.75 + 2e-5 * elev) * Ra,
    Rs_Rso = pmin(pmax(rad / Rso, 0.25), 1),
    Rns    = (1 - 0.23) * rad,
    Rnl    = 4.903e-9 * (((max_temp + 273.16)^4 + (min_temp + 273.16)^4) / 2) *
      (0.34 - 0.14 * sqrt(ea)) * (1.35 * Rs_Rso - 0.35),
    Rn     = Rns - Rnl,
    
    # FAO-56 Penman-Monteith reference ET (mm/day)
    PET = (0.408 * Delta * Rn + gamma * (900 / (Tmean + 273)) * u2 * vpd) /
      (Delta + gamma * (1 + 0.34 * u2)),
    PET = pmax(PET, 0),
    
    CWB = prec - PET
  )

head(merged_data)

# Quick plausibility check against the observed data for the same sites
merged_data %>%
  dplyr::filter(Climate_model == "Reference") %>%
  group_by(Simulation_site, year) %>%
  summarise(prec = sum(prec), Eto = sum(PET), CWB = sum(CWB), .groups = "drop") %>%
  group_by(Simulation_site) %>%
  summarise(across(prec:CWB, mean))

###############################################################################
### Period summaries: whole year (_WY), growing period (_GP), summer (_MSEA)
###############################################################################

# One helper so all three periods are computed identically.
summ_period <- function(d, sfx, indices = TRUE) {
  d <- d %>% arrange(Simulation_site, Climate_model, Scenario, date)
  
  out <- d %>%
    group_by(Simulation_site, Climate_model, Scenario, year) %>%
    summarise(
      TMean = mean(Tmean),
      Prec  = sum(prec),
      VPD   = mean(vpd),
      Eto   = sum(PET),
      CWB   = sum(CWB),
      
      # CDD: Maximum length of dry spell (prec < 1 mm)
      CDD = {
        rle_dry <- rle(prec < 1)
        if (any(rle_dry$values)) max(rle_dry$lengths[rle_dry$values]) else 0
      },
      # SDII: Simple precipitation intensity index (sum on wet days / wet days)
      SDII = {
        wet_days_prec <- prec[prec >= 1]
        if (length(wet_days_prec) > 0) sum(wet_days_prec) / length(wet_days_prec) else 0
      },
      # Wdc: Count of days where precipitation is >= 1 mm
      Wdc = sum(prec >= 1),
      
      n_days = n(),
      .groups = "drop"
    )
  
  if (!indices) out <- out %>% dplyr::select(-CDD, -SDII, -Wdc)
  
  out %>%
    rename_with(~ paste0(.x, sfx),
                -c(Simulation_site, Climate_model, Scenario, year))
}

# Whole year: as before, no precipitation indices
wy_sum <- summ_period(merged_data, "_WY", indices = FALSE)
head(wy_sum)

# NEW: growing period, April to September
gp <- merged_data %>% dplyr::filter(month %in% 4:9)
gp_sum <- summ_period(gp, "_GP")
head(gp_sum)

# Midsummer to early autumn, July to September
msea <- merged_data %>% dplyr::filter(month %in% 7:9)
msea_sum <- summ_period(msea, "_MSEA")
head(msea_sum)

###############################################################################
### Merge the periods
###############################################################################
join_keys <- c("Simulation_site", "Climate_model", "Scenario", "year")

climate_merged <- wy_sum %>%
  left_join(gp_sum,   by = join_keys) %>%
  left_join(msea_sum, by = join_keys)

# Completeness check: 365/366 days per year, 183 for GP, 92 for MSEA
climate_merged %>%
  dplyr::filter(n_days_WY < 365 | n_days_GP < 183 | n_days_MSEA < 92) %>%
  dplyr::select(Simulation_site, Climate_model, Scenario, year,
                n_days_WY, n_days_GP, n_days_MSEA)

climate_merged <- climate_merged %>% dplyr::select(-starts_with("n_days"))

# Backward compatibility: the summer indices keep their original unsuffixed
# names, so anything downstream that expects CDD / SDII / Wdc still works.
# Comment these three lines out if you would rather have CDD_MSEA etc.
climate_merged <- climate_merged %>%
  rename(CDD = CDD_MSEA, SDII = SDII_MSEA, Wdc = Wdc_MSEA)

climate_merged <- climate_merged %>% mutate(Country = "Spain")

climate_merged <- climate_merged %>% dplyr::filter(between(year, 2007, 2100))
head(climate_merged)
tail(climate_merged)
min(climate_merged$year)
names(climate_merged)

writexl::write_xlsx(climate_merged,
                    file.path(out_dir, "ClimateSummaries_Spain.xlsx"))