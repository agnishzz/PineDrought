#Date 11th August 2025, revised 7th October 2026
#The code is being prepared to analyse the climate from Spanish plots
#REVISION: evapotranspiration is now FAO-56 Penman-Monteith (was temperature-based).
#          CWB follows from the new PET. Eto is now also summarised and anomalised.
#          Simulation sites are renamed to real site names. Output directory changed.

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


#### Add a simulation site column, climate model column and rcp column -------
# CHANGED: S0 = Prades, S1..S4 = Arcalis1..4. Anchored with ^ so the pattern
# can only match the prefix we added, never a substring of the file name.
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

# Sanity check: every row must have a site, model and scenario
stopifnot(!anyNA(merged_data$Simulation_site),
          !anyNA(merged_data$Climate_model),
          !anyNA(merged_data$Scenario))
dplyr::count(merged_data, Simulation_site, Climate_model, Scenario)

write.csv(merged_data, file.path(out_dir, "ClimateDataSpain_whole.csv"), row.names = FALSE)

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

# Quick plausibility check: mean annual PET and CWB per site, reference climate
merged_data %>%
  dplyr::filter(Climate_model == "Reference") %>%
  group_by(Simulation_site, year) %>%
  summarise(prec = sum(prec), PET = sum(PET), CWB = sum(CWB), .groups = "drop") %>%
  group_by(Simulation_site) %>%
  summarise(across(prec:CWB, mean))

# Isolate midsummer to early autumn (July, August, September)
msea <- merged_data %>% dplyr::filter(month %in% 7:9)
head(msea)

# Define the baseline period
baseline_start <- 2007
baseline_end   <- 2015

############################ CWB anomalies ####################################
sum_msea <- msea %>%
  dplyr::group_by(Simulation_site, Climate_model, Scenario, year) %>%
  dplyr::summarise(CWB_sum = sum(CWB), .groups = "drop")
head(sum_msea)

sum_msea_anomalies_baseline <- sum_msea %>%
  group_by(Simulation_site, Climate_model, Scenario) %>%
  mutate(
    baseline_mean = mean(CWB_sum[year >= baseline_start & year <= baseline_end]),
    CWB_anomaly   = CWB_sum - baseline_mean
  ) %>%
  ungroup()

print(sum_msea_anomalies_baseline)

ggplot(sum_msea_anomalies_baseline, aes(x = year, y = CWB_anomaly)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(aes(color = Climate_model,
                group = interaction(Simulation_site, Scenario, Climate_model)),
            alpha = 0.6) +
  stat_summary(fun = "mean", geom = "line", color = "black",
               linetype = "dotted", linewidth = 1) +
  facet_grid(Simulation_site ~ Scenario) +
  labs(
    title    = "Climatic Water Balance (CWB) Anomaly by Scenario",
    subtitle = "Penman-Monteith PET; anomalies relative to the 2007-2015 baseline",
    x = "Year", y = "CWB Anomaly (mm)", color = "Climate model"
  ) +
  theme_bw()

writexl::write_xlsx(sum_msea_anomalies_baseline,
                    file.path(out_dir, "CWB_Anomalies_Spain_fut_2.xlsx"))

################### NEW: Eto (reference ET) anomalies ##########################
sum_msea_e <- msea %>%
  dplyr::group_by(Simulation_site, Climate_model, Scenario, year) %>%
  dplyr::summarise(Eto_sum = sum(PET), .groups = "drop")
head(sum_msea_e)

sum_msea_anomalies_baseline_e <- sum_msea_e %>%
  group_by(Simulation_site, Climate_model, Scenario) %>%
  mutate(
    baseline_mean = mean(Eto_sum[year >= baseline_start & year <= baseline_end]),
    Eto_anomaly   = Eto_sum - baseline_mean
  ) %>%
  ungroup()

print(sum_msea_anomalies_baseline_e)

ggplot(sum_msea_anomalies_baseline_e, aes(x = year, y = Eto_anomaly)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(aes(color = Climate_model,
                group = interaction(Simulation_site, Scenario, Climate_model)),
            alpha = 0.6) +
  stat_summary(fun = "mean", geom = "line", color = "black",
               linetype = "dotted", linewidth = 1) +
  facet_grid(Simulation_site ~ Scenario) +
  labs(
    title    = "Reference Evapotranspiration (Eto) Anomaly by Scenario",
    subtitle = "FAO-56 Penman-Monteith; anomalies relative to the 2007-2015 baseline",
    x = "Year", y = "Eto Anomaly (mm)", color = "Climate model"
  ) +
  theme_bw()

writexl::write_xlsx(sum_msea_anomalies_baseline_e,
                    file.path(out_dir, "Eto_Anomalies_Spain_fut.xlsx"))

################## Precipitation indices (unchanged) ###########################
climate_indices <- msea %>%
  group_by(year, Simulation_site, Climate_model, Scenario) %>%
  summarise(
    # CDD: Maximum length of dry spell (prec < 1 mm)
    CDD = {
      rle_dry <- rle(prec < 1)
      if (any(rle_dry$values)) max(rle_dry$lengths[rle_dry$values]) else 0
    },
    # SDII: Simple precipitation intensity index (sum on wet days / count of wet days)
    SDII = {
      wet_days_prec <- prec[prec >= 1]
      if (length(wet_days_prec) > 0) sum(wet_days_prec) / length(wet_days_prec) else 0
    },
    # Wdc: Count of days where precipitation is >= 1 mm
    Wdc = sum(prec >= 1),
    .groups = 'drop'
  )

print(climate_indices)

climate_indices_anomalies <- climate_indices %>%
  group_by(Simulation_site, Climate_model, Scenario) %>%
  mutate(
    CDD_baseline_mean  = mean(CDD[year >= baseline_start & year <= baseline_end]),
    CDD_anomaly        = CDD - CDD_baseline_mean,
    SDII_baseline_mean = mean(SDII[year >= baseline_start & year <= baseline_end]),
    SDII_anomaly       = SDII - SDII_baseline_mean,
    Wdc_baseline_mean  = mean(Wdc[year >= baseline_start & year <= baseline_end]),
    Wdc_anomaly        = Wdc - Wdc_baseline_mean
  ) %>%
  ungroup()

print(climate_indices_anomalies)

ggplot(climate_indices_anomalies, aes(x = year, y = CDD_anomaly)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(aes(color = Climate_model,
                group = interaction(Simulation_site, Scenario, Climate_model)),
            alpha = 0.6) +
  stat_summary(fun = "mean", geom = "line", color = "black",
               linetype = "dotted", linewidth = 1) +
  facet_grid(Simulation_site ~ Scenario) +
  labs(
    title    = "Consecutive Dry Days (CDD) Anomaly by Scenario",
    subtitle = "Anomalies are relative to the 2007-2015 baseline period",
    x = "Year", y = "CDD Anomaly", color = "Climate model"
  ) +
  theme_bw()

writexl::write_xlsx(climate_indices_anomalies,
                    file.path(out_dir, "CI_Anomalies_Spain_fut.xlsx"))

############################# Precipitation anomaly ############################
sum_msea_p <- msea %>%
  dplyr::group_by(Simulation_site, Climate_model, Scenario, year) %>%
  dplyr::summarise(prec_sum = sum(prec), .groups = "drop")

sum_msea_anomalies_baseline_p <- sum_msea_p %>%
  group_by(Simulation_site, Climate_model, Scenario) %>%
  mutate(
    baseline_mean = mean(prec_sum[year >= baseline_start & year <= baseline_end]),
    prec_anomaly  = prec_sum - baseline_mean
  ) %>%
  ungroup()

print(sum_msea_anomalies_baseline_p)

writexl::write_xlsx(sum_msea_anomalies_baseline_p,
                    file.path(out_dir, "prec_Anomalies_Spain_fut.xlsx"))

############################### MaxTemp anomaly ################################
sum_msea_t <- msea %>%
  dplyr::group_by(Simulation_site, Climate_model, Scenario, year) %>%
  dplyr::summarise(temp_av = mean(max_temp), .groups = "drop")

sum_msea_anomalies_baseline_t <- sum_msea_t %>%
  group_by(Simulation_site, Climate_model, Scenario) %>%
  mutate(
    baseline_mean = mean(temp_av[year >= baseline_start & year <= baseline_end]),
    temp_anomaly  = temp_av - baseline_mean
  ) %>%
  ungroup()

print(sum_msea_anomalies_baseline_t)

writexl::write_xlsx(sum_msea_anomalies_baseline_t,
                    file.path(out_dir, "temp_Anomalies_Spain_fut.xlsx"))

####################### Making a single dataset ################################
cwb  <- readxl::read_excel(file.path(out_dir, "CWB_Anomalies_Spain_fut_2.xlsx"))
eto  <- readxl::read_excel(file.path(out_dir, "Eto_Anomalies_Spain_fut.xlsx"))
prec <- readxl::read_excel(file.path(out_dir, "prec_Anomalies_Spain_fut.xlsx"))
temp <- readxl::read_excel(file.path(out_dir, "temp_Anomalies_Spain_fut.xlsx"))
CI   <- readxl::read_excel(file.path(out_dir, "CI_Anomalies_Spain_fut.xlsx"))

head(cwb); head(eto); head(prec); head(temp); head(CI)

join_keys <- c("year", "Simulation_site", "Climate_model", "Scenario")

all_anomalies <- cwb %>%
  dplyr::select(-ends_with("baseline_mean")) %>%
  left_join(eto  %>% dplyr::select(-ends_with("baseline_mean")), by = join_keys) %>%
  left_join(prec %>% dplyr::select(-ends_with("baseline_mean")), by = join_keys) %>%
  left_join(temp %>% dplyr::select(-ends_with("baseline_mean")), by = join_keys) %>%
  left_join(CI   %>% dplyr::select(-ends_with("baseline_mean")), by = join_keys)

print(all_anomalies)

all_anomalies <- all_anomalies %>% mutate(Country = "Spain")

# Check: no rows lost or duplicated by the joins
stopifnot(nrow(all_anomalies) == nrow(cwb))

print(all_anomalies)
writexl::write_xlsx(all_anomalies, file.path(out_dir, "Spain_clim_fut.xlsx"))