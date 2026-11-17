library(DBI)
library(RSQLite)
library(dplyr)
library(purrr)

clim_dir <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate"
files    <- list.files(clim_dir, pattern = "\\.sqlite$", full.names = TRUE)

# Quick look at which tables each file contains
walk(files, function(f) {
  con <- dbConnect(SQLite(), f)
  cat(basename(f), ":", paste(dbListTables(con), collapse = ", "), "\n")
  dbDisconnect(con)
})

read_site <- function(f) {
  site <- sub("_[0-9_]+\\.sqlite$", "", basename(f))   # strip _07_15 / _1998_2020
  con  <- dbConnect(SQLite(), f)
  on.exit(dbDisconnect(con))
  map_dfr(dbListTables(con), function(tb) {
    dbReadTable(con, tb) |> mutate(Table = tb)
  }) |>
    mutate(Site = site, .before = 1)
}

clim_all <- map_dfr(files, read_site)

# Drop the Table column if every file had only one table
if (n_distinct(clim_all$Table) == 1) clim_all <- select(clim_all, -Table)

count(clim_all, Site)
glimpse(clim_all)
####

library(dplyr)
library(lubridate)

# ---- 1. Site metadata (decimal degrees, metres a.s.l.) -------------------
site_meta <- tibble::tribble(
  ~Site,          ~lat,        ~lon,         ~elev,
  "CH_Visp",      46 + 18/60,  7 + 52/60,    695,
  "SP_Arcalis1",  42.36254117, 1.095519333,  660.3333333,
  "SP_Arcalis2",  42.34864743, 1.088869286,  987,
  "SP_Arcalis3",  42.3673052,  1.049716,     1216,
  "SP_Arcalis4",  42.34418886, 1.097959857,  1321.571429,
  "SP_Prades",    41.33691,    1.01,         928.3
)
stopifnot("Fill in lat/elev for every site" = !anyNA(site_meta))

# ---- 2. Helper functions (FAO-56, Allen et al. 1998) ---------------------
e0 <- function(T) 0.6108 * exp(17.27 * T / (T + 237.3))        # kPa

ra_fao <- function(lat_deg, doy) {                               # MJ m-2 d-1
  phi   <- lat_deg * pi / 180
  dr    <- 1 + 0.033 * cos(2 * pi * doy / 365)
  delta <- 0.409 * sin(2 * pi * doy / 365 - 1.39)
  ws    <- acos(pmin(pmax(-tan(phi) * tan(delta), -1), 1))
  (24 * 60 / pi) * 0.0820 * dr *
    (ws * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(ws))
}

u2 <- 2  # m/s, FAO-56 default when wind is unavailable

# ---- 3. Daily PET and CWB -----------------------------------------------
clim_pet <- clim_all |>
  dplyr::select(-any_of("Table")) |>
  left_join(site_meta, by = "Site") |>
  mutate(
    date  = make_date(year, month, day),
    doy   = yday(date),
    tmean = (max_temp + min_temp) / 2,
    
    # psychrometrics
    P_atm = 101.3 * ((293 - 0.0065 * elev) / 293)^5.26,          # kPa
    gamma = 0.000665 * P_atm,
    Delta = 4098 * e0(tmean) / (tmean + 237.3)^2,
    es    = (e0(max_temp) + e0(min_temp)) / 2,
    ea    = pmax(es - vpd, 0.01),        # actual vapour pressure from VPD
    
    # radiation
    Ra    = ra_fao(lat, doy),
    Rso   = (0.75 + 2e-5 * elev) * Ra,
    Rs_Rso = pmin(pmax(rad / Rso, 0.25), 1),
    Rns   = (1 - 0.23) * rad,
    Rnl   = 4.903e-9 * (((max_temp + 273.16)^4 + (min_temp + 273.16)^4) / 2) *
      (0.34 - 0.14 * sqrt(ea)) * (1.35 * Rs_Rso - 0.35),
    Rn    = Rns - Rnl,
    
    # FAO-56 Penman-Monteith reference ET (mm/day); VPD used directly
    PET_PM = (0.408 * Delta * Rn + gamma * (900 / (tmean + 273)) * u2 * vpd) /
      (Delta + gamma * (1 + 0.34 * u2)),
    PET_PM = pmax(PET_PM, 0),
    
    # Hargreaves-Samani check (mm/day)
    PET_HS = 0.0023 * 0.408 * Ra * (tmean + 17.8) *
      sqrt(pmax(max_temp - min_temp, 0)),
    
    CWB    = prec - PET_PM,
    CWB_HS = prec - PET_HS
  ) |>
   dplyr:: select(Site, date, year, month, day, min_temp, max_temp, prec, rad, vpd,
         PET_PM, PET_HS, CWB, CWB_HS)

# ---- 4. Sanity checks ---------------------------------------------------
clim_pet |>
  group_by(Site, year) |>
  summarise(prec = sum(prec), PET_PM = sum(PET_PM), PET_HS = sum(PET_HS),
            CWB = sum(CWB), .groups = "drop") |>
  group_by(Site) |>
  summarise(across(prec:CWB, mean))
###################################################################################
#####################################################################################

clim_dir <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate"
out_file <- file.path(clim_dir, "AllSites_PET_CWB.sqlite")

out <- clim_pet |>
  dplyr::mutate(date = format(date, "%Y-%m-%d"))

con <- dbConnect(SQLite(), out_file)
dbWriteTable(con, "climate", out, overwrite = TRUE)

# Verify
dbListTables(con)
dbGetQuery(con, "SELECT Site, COUNT(*) AS n, MIN(date) AS first, MAX(date) AS last
                 FROM climate GROUP BY Site")
dbDisconnect(con)

######################################################################################
######################################################################################
library(dplyr)

wet_thr <- 1   # mm

# Longest run of TRUE values in a logical vector
longest_run <- function(x) {
  r <- rle(x)
  if (!any(r$values)) 0L else max(r$lengths[r$values])
}

summ_period <- function(d) {
  d |>
    arrange(Site, date) |>                      # CDD needs days in order
    group_by(Site, year) |>
    summarise(
      MaxTemp = mean(max_temp),
      MinTemp = mean(min_temp),
      VPD     = mean(vpd),
      Prec    = sum(prec),
      Eto     = sum(PET_PM),
      CWB     = sum(CWB),
      WDcount = sum(prec >= wet_thr),
      SDII    = if_else(WDcount > 0, Prec / WDcount, NA_real_),
      CDD     = longest_run(prec < wet_thr),
      n_days  = n(),
      .groups = "drop"
    )
}

yr <- summ_period(clim_pet)

gp <- clim_pet |>
  dplyr::filter(month %in% 4:9) |>
  summ_period() |>
  rename_with(~ paste0(.x, "_GP"), -c(Site, year))

ss <- clim_pet |>
  dplyr::filter(month %in% 7:9) |>
  summ_period() |>
  rename_with(~ paste0(.x, "_SS"), -c(Site, year))

clim_annual <- yr |>
  left_join(gp, by = c("Site", "year")) |>
  left_join(ss, by = c("Site", "year"))

# Completeness check: full year 365/366, GP 183, SS 92 days
clim_annual |>
  dplyr::filter(n_days < 365 | n_days_GP < 183 | n_days_SS < 92) |>
  dplyr::select(Site, year, n_days, n_days_GP, n_days_SS)

clim_annual <- clim_annual |> dplyr::select(-starts_with("n_days"))
glimpse(clim_annual)
head(clim_annual)
writexl::write_xlsx(clim_annual, "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/AllCLimateSummary.xlsx")

####################################################################################################################################
#Calculate anomalies
vars <- setdiff(names(clim_annual), c("Site", "year"))

# Site means (long-term reference per site)
site_means <- clim_annual |>
  group_by(Site) |>
  summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Anomalies: yearly value minus site mean
clim_anom <- clim_annual |>
  group_by(Site) |>
  mutate(across(all_of(vars),
                ~ .x - mean(.x, na.rm = TRUE),
                .names = "AN_{.col}")) |>
  ungroup() |>
  dplyr::select(Site, year, starts_with("AN_"))

glimpse(clim_anom)

# Check: anomalies should average to ~0 within each site
clim_anom |>
  group_by(Site) |>
  summarise(across(starts_with("AN_"), ~ round(mean(.x, na.rm = TRUE), 10))) |>
  dplyr::select(Site, AN_Prec, AN_CWB_GP, AN_CDD_SS)

writexl::write_xlsx(clim_anom, "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/AllCLimateAnomalies.xlsx")


library(dplyr)

clim_lag <- clim_anom |>
  dplyr::select(Site, year) |>
  left_join(
    clim_anom |>
      mutate(year = year + 1) |>
      rename_with(~ sub("^AN_", "YL_", .x), starts_with("AN_")),
    by = c("Site", "year")
  )

glimpse(clim_lag)

# Check: 1999 YL_ values should equal 1998 AN_ values for Visp
bind_rows(
  clim_anom |> dplyr::filter(Site == "CH_Visp", year == 1998) |>
    dplyr::select(Site, year, AN_Prec, AN_CWB_GP),
  clim_lag  |> dplyr::filter(Site == "CH_Visp", year == 1999) |>
    dplyr::select(Site, year, YL_Prec, YL_CWB_GP)
)

# First year of each site has no previous year, so it should be all NA
clim_lag |>
  group_by(Site) |>
  slice_min(year) |>
  dplyr::select(Site, year, YL_Prec, YL_CWB)

writexl::write_xlsx(clim_lag, "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/ObsMor_anomaly.xlsx")
