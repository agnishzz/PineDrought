library(tidyverse)

data <- readxl::read_xlsx("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/ObsMor_anomaly.xlsx")

# ---------------------------------------------------------------
# 1. Switzerland data
# ---------------------------------------------------------------
ch <- data %>%
  dplyr::filter(Country == "Switzerland") %>%
  dplyr::select(Site, Year, ObsMor, SimMor) %>%
  tidyr::drop_na(ObsMor, SimMor) %>%
  dplyr::arrange(Year)

ch

# ---------------------------------------------------------------
# 2. Nash-Sutcliffe Efficiency
# ---------------------------------------------------------------
NSE_CH <- 1 - sum((ch$ObsMor - ch$SimMor)^2) /
  sum((ch$ObsMor - mean(ch$ObsMor))^2)

R2_CH   <- cor(ch$ObsMor, ch$SimMor)^2
bias_CH <- mean(ch$SimMor - ch$ObsMor)

NSE_CH
R2_CH
bias_CH

# ---------------------------------------------------------------
# 3. Contingency table, ABSOLUTE threshold (3% for both series)
# ---------------------------------------------------------------
obs_thr <- 3

ch$obs_event     <- ch$ObsMor > obs_thr
ch$sim_event_abs <- ch$SimMor > obs_thr

hits_abs    <- sum( ch$obs_event &  ch$sim_event_abs)
fa_abs      <- sum(!ch$obs_event &  ch$sim_event_abs)
miss_abs    <- sum( ch$obs_event & !ch$sim_event_abs)
cn_abs      <- sum(!ch$obs_event & !ch$sim_event_abs)

table(Observed = ch$obs_event, Simulated = ch$sim_event_abs)

POD_abs  <- hits_abs / (hits_abs + miss_abs)
FAR_abs  <- fa_abs / (hits_abs + fa_abs)
POFD_abs <- fa_abs / (fa_abs + cn_abs)
PSS_abs  <- POD_abs - POFD_abs
HSS_abs  <- 2 * (hits_abs * cn_abs - fa_abs * miss_abs) /
  ((hits_abs + miss_abs) * (miss_abs + cn_abs) +
     (hits_abs + fa_abs)   * (fa_abs   + cn_abs))

# ---------------------------------------------------------------
# 4. Contingency table, FREQUENCY-MATCHED threshold
#    (simulated threshold gives same number of events as observed)
# ---------------------------------------------------------------
freq_obs  <- mean(ch$obs_event)
sim_thr_q <- unname(quantile(ch$SimMor, probs = 1 - freq_obs))

ch$sim_event_q <- ch$SimMor > sim_thr_q

hits_q <- sum( ch$obs_event &  ch$sim_event_q)
fa_q   <- sum(!ch$obs_event &  ch$sim_event_q)
miss_q <- sum( ch$obs_event & !ch$sim_event_q)
cn_q   <- sum(!ch$obs_event & !ch$sim_event_q)

table(Observed = ch$obs_event, Simulated = ch$sim_event_q)

POD_q  <- hits_q / (hits_q + miss_q)
FAR_q  <- fa_q / (hits_q + fa_q)
POFD_q <- fa_q / (fa_q + cn_q)
PSS_q  <- POD_q - POFD_q
HSS_q  <- 2 * (hits_q * cn_q - fa_q * miss_q) /
  ((hits_q + miss_q) * (miss_q + cn_q) +
     (hits_q + fa_q)   * (fa_q   + cn_q))

# ---------------------------------------------------------------
# 5. Summary table
# ---------------------------------------------------------------
skill_CH <- tibble::tibble(
  Threshold     = c("Absolute (3%)", "Frequency-matched"),
  Sim_threshold = c(obs_thr, sim_thr_q),
  Hits          = c(hits_abs, hits_q),
  False_alarms  = c(fa_abs, fa_q),
  Misses        = c(miss_abs, miss_q),
  Correct_neg   = c(cn_abs, cn_q),
  POD           = c(POD_abs, POD_q),
  FAR           = c(FAR_abs, FAR_q),
  POFD          = c(POFD_abs, POFD_q),
  PSS           = c(PSS_abs, PSS_q),
  HSS           = c(HSS_abs, HSS_q),
  NSE           = NSE_CH,
  R2            = R2_CH,
  Bias          = bias_CH
)

print(skill_CH, width = Inf)

# ---------------------------------------------------------------
# 6. Year-by-year classification (appendix table)
# ---------------------------------------------------------------
ch$class_abs <- ifelse( ch$obs_event &  ch$sim_event_abs, "Hit",
                        ifelse(!ch$obs_event &  ch$sim_event_abs, "False alarm",
                               ifelse( ch$obs_event & !ch$sim_event_abs, "Miss", "Correct negative")))

ch$class_q   <- ifelse( ch$obs_event &  ch$sim_event_q, "Hit",
                        ifelse(!ch$obs_event &  ch$sim_event_q, "False alarm",
                               ifelse( ch$obs_event & !ch$sim_event_q, "Miss", "Correct negative")))

print(ch[, c("Site", "Year", "ObsMor", "SimMor", "class_abs", "class_q")], n = Inf)

#######################################################################################################
#######################################################################################################
#######################################################################################################

