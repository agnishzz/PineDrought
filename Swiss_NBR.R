# =============================================================================
# Switzerland (Visp): annual mortality allocation using NBR + NDVI AGREEMENT
# -----------------------------------------------------------------------------
# Same hazard model as before:
#   m_t = 1 - exp(-beta * s_t),  beta = ln(N0/NT) / sum(s_t)   (74 -> 20 trees exact)
# What changes is the signal s_t. Canopy-loss drops are computed per index as
#   d_t = max(0, max(X_{t-1}, X_{t-2}) - X_t)                   (2-yr-max drop)
# and combined with a fixed set of rules (defined before looking at results):
#
#   NBR         NBR drop only (reference; previously validated)
#   NDVI        NDVI drop only (reference)
#   AND_NBR     year counts only if NBR AND NDVI both drop; weight = NBR drop
#   AND_meanz   both drop; weight = mean of standardised drops
#   AND_geomz   both drop; weight = geometric mean of standardised drops
#   ANDlag_NBR  NBR drops in t AND NDVI drops in t or t+1 (needle browning lag);
#               weight = NBR drop
#   OR_meanz    either drops; weight = mean of standardised drops
#   OR_maxz     either drops; weight = larger standardised drop
#
# Standardised drop z = d_t / SD(year-to-year differences of that index), so the
# two indices are on a comparable noise scale before combining.
#
# Validation: annual mortality % from each rule vs observed ObsMor (%).
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
})

# ---- Settings ---------------------------------------------------------------
infile     <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/RS/NBR_swiss.xlsx"          # Year, NBR, ObsMor, NDVI, NDWI
out_dir    <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/RS/swiss_consensus_output"
N0         <- 74
NT         <- 20
year_start <- 1999
year_end   <- 2020
window     <- 2
dir.create(out_dir, showWarnings = FALSE)

# ---- Functions --------------------------------------------------------------
drop2 <- function(x, w = 2) {
  s <- rep(NA_real_, length(x))
  for (t in seq_along(x)) if (t > w) s[t] <- max(0, max(x[(t - w):(t - 1)]) - x[t])
  s
}

# All consensus signals for one pair of NBR / NDVI series (ordered by year)
consensus_signals <- function(nbr, ndvi, w = 2) {
  dn <- drop2(nbr, w);  dv <- drop2(ndvi, w)
  zn <- dn / sd(diff(nbr));  zv <- dv / sd(diff(ndvi))
  dn0 <- ifelse(is.na(dn), 0, dn); dv0 <- ifelse(is.na(dv), 0, dv)
  zn0 <- ifelse(is.na(zn), 0, zn); zv0 <- ifelse(is.na(zv), 0, zv)
  both   <- dn0 > 0 & dv0 > 0
  either <- dn0 > 0 | dv0 > 0
  dv_next <- c(dv0[-1], 0)                       # NDVI drop in year t+1
  lagconf <- dn0 > 0 & (dv0 > 0 | dv_next > 0)
  list(
    NBR        = dn0,
    NDVI       = dv0,
    AND_NBR    = ifelse(both, dn0, 0),
    AND_meanz  = ifelse(both, (zn0 + zv0) / 2, 0),
    AND_geomz  = ifelse(both, sqrt(zn0 * zv0), 0),
    ANDlag_NBR = ifelse(lagconf, dn0, 0),
    OR_meanz   = ifelse(either, (zn0 + zv0) / 2, 0),
    OR_maxz    = ifelse(either, pmax(zn0, zv0), 0))
}

cum_round <- function(d) diff(c(0, round(cumsum(d))))

allocate_hazard <- function(s, N0, NT) {
  s[is.na(s)] <- 0
  if (sum(s) <= 0) stop("Signal is zero in every year.")
  beta  <- log(N0 / NT) / sum(s)
  rate  <- 1 - exp(-beta * s)
  N_end <- N0 * cumprod(1 - rate)
  N_beg <- c(N0, head(N_end, -1))
  data.frame(signal = s, rate = rate, N_start = N_beg, deaths = N_beg - N_end,
             deaths_int = cum_round(N_beg - N_end), mort_pct = 100 * rate)
}

# ---- Data -------------------------------------------------------------------
dat <- as.data.frame(read_excel(infile))
dat <- dat[order(dat$Year), ]
sel <- dat$Year >= year_start & dat$Year <= year_end
yrs <- dat$Year[sel]

sig <- consensus_signals(dat$NBR, dat$NDVI, window)
alloc <- lapply(sig, function(s) allocate_hazard(s[sel], N0, NT))

# ---- Validation vs observed mortality (%) ------------------------------------
obs <- dat$ObsMor[sel]
val <- do.call(rbind, lapply(names(alloc), function(m) {
  a  <- alloc[[m]]
  ct <- cor.test(a$mort_pct, obs)
  data.frame(Method = m,
             signal_years = sum(a$signal > 0),
             R2           = unname(ct$estimate)^2,
             pearson_r    = unname(ct$estimate),
             p_value      = ct$p.value,
             spearman_rho = cor(a$mort_pct, obs, method = "spearman"),
             RMSE_pct     = sqrt(mean((a$mort_pct - obs)^2)))
}))
val <- val[order(-val$R2), ]
write.csv(val, file.path(out_dir, "swiss_consensus_validation.csv"), row.names = FALSE)
cat("Validation against observed annual mortality (ObsMor, %), 1999-2020, n = ", length(obs), "\n", sep = "")
vv <- val; vv[, 3:7] <- round(vv[, 3:7], 3); print(vv, row.names = FALSE)

# ---- Output: annual mortality % per method (long + wide) ----------------------
wide <- data.frame(Year = yrs, ObsMor = obs,
                   sapply(alloc, function(a) a$mort_pct))
write.csv(wide, file.path(out_dir, "swiss_annual_mortpct_consensus_wide.csv"), row.names = FALSE)

long <- do.call(rbind, lapply(names(alloc), function(m)
  data.frame(Year = yrs, Method = m, alloc[[m]], ObsMor = obs)))
write.csv(long, file.path(out_dir, "swiss_annual_mortpct_consensus_long.csv"), row.names = FALSE)

# Which years each rule keeps
kept <- data.frame(Year = yrs, NBR_drop = round(sig$NBR[sel], 3), NDVI_drop = round(sig$NDVI[sel], 3),
                   sapply(sig, function(s) as.integer(s[sel] > 0)), ObsMor = round(obs, 1))
write.csv(kept, file.path(out_dir, "swiss_consensus_years_kept.csv"), row.names = FALSE)

# ---- Figures ----------------------------------------------------------------
lab <- setNames(sprintf("%s (R² = %.2f)", val$Method, val$R2), val$Method)
long$Method <- factor(long$Method, levels = val$Method)
pf <- ggplot(long, aes(Year)) +
  geom_col(aes(y = ObsMor), fill = "grey80", width = 0.8) +
  geom_line(aes(y = mort_pct), colour = "#2C6E9B", linewidth = 0.7) +
  geom_point(aes(y = mort_pct), colour = "#2C6E9B", size = 1.3) +
  facet_wrap(~ Method, ncol = 2, labeller = as_labeller(lab)) +
  labs(x = NULL, y = "Annual mortality (%)",
       title = "Visp: allocated annual mortality (line) vs observed ObsMor (bars)") +
  theme_bw()
pf
ggsave(file.path(out_dir, "swiss_consensus_timeseries.png"), pf, width = 10, height = 9, dpi = 300)

cat("\nOutputs written to", normalizePath(out_dir), "\n")
