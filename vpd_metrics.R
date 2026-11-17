library(here)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(fields)

path0 <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/AllSites_PET_CWB.sqlite"

file0 <- here::here(path0)




# Read in data using the RSQLite package
sqlite.driver <- RSQLite::dbDriver("SQLite")
db0 <-
  RSQLite::dbConnect(sqlite.driver, dbname = file0) # connect to the file
tables.in.the.file0 <-
  RSQLite::dbListTables(db0) # explore the tables in the file1
print(tables.in.the.file0)
climate0 <- RSQLite::dbReadTable(db0, "climate")
RSQLite::dbDisconnect(db0)

climate <- climate0

head(climate)
# ---------------------------------------------------------------------------
# Extreme VPD metric (reviewer response): mean VPD of the days at or above
# the 95th percentile, per Site x year x period.
#   GP = growing period (April-September)
#   SS = summer season  (July-September)
# The 95th-percentile threshold is computed within each Site x year x period,
# ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)

periods <- list(GP = 4:9, SS = 7:9)

vpd_p95 <- bind_rows(lapply(names(periods), function(p) {
  climate %>%
    filter(month %in% periods[[p]], !is.na(vpd)) %>%
    mutate(period = p)
})) %>%
  group_by(Site, year, period) %>%
  mutate(thr_p95 = quantile(vpd, 0.95, type = 7)) %>%
  summarise(
    vpd_thr_p95  = first(thr_p95),              # threshold (kPa)
    vpd_p95_mean = mean(vpd[vpd >= thr_p95]),   # metric (kPa)
    n_days_p95   = sum(vpd >= thr_p95),         # days used
    n_days_total = n(),                         # completeness check
    .groups = "drop"
  ) %>%
  arrange(Site, period, year)

# Wide format for joining to the mortality series (Site, year)
vpd_p95_wide <- vpd_p95 %>%
  select(Site, year, period, vpd_p95_mean) %>%
  pivot_wider(names_from = period, values_from = vpd_p95_mean,
              names_prefix = "VPD95_")          # -> VPD95_GP, VPD95_SS

# Flag incomplete years (GP should have 183 days, SS 92)
incomplete <- vpd_p95 %>%
  filter((period == "GP" & n_days_total < 183) |
           (period == "SS" & n_days_total < 92))
if (nrow(incomplete) > 0) {
  message("Incomplete site-years:"); print(incomplete)
}

vpd_p95_wide
##################################
# ---------------------------------------------------------------------------
# Anomalies: deviation from each Site's long-term mean (all available years)
# ---------------------------------------------------------------------------
vpd_p95_anom <- vpd_p95_wide %>%
  group_by(Site) %>%
  mutate(across(c(VPD95_GP, VPD95_SS),
                ~ .x - mean(.x, na.rm = TRUE),
                .names = "{.col}_anom")) %>%
  ungroup() %>%
  select(Site, year, VPD95_GP_anom, VPD95_SS_anom)

# Year-lagged anomalies: previous year's anomaly assigned to current year
vpd_p95_anom_lag <- vpd_p95_anom %>%
  arrange(Site, year) %>%
  group_by(Site) %>%
  mutate(across(c(VPD95_GP_anom, VPD95_SS_anom),
                ~ if_else(year - lag(year) == 1, lag(.x), NA_real_),
                .names = "{.col}_lag1")) %>%
  ungroup() %>%
  select(Site, year, VPD95_GP_anom_lag1, VPD95_SS_anom_lag1)

writexl::write_xlsx(vpd_p95_wide, "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/VPD_p95.xlsx")
writexl::write_xlsx(vpd_p95_anom, "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/VPD_anom.xlsx")
writexl::write_xlsx(vpd_p95_anom_lag, "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/VPD_anom_lag.xlsx")
