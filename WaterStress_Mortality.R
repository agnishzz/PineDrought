library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

in_file <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/WR_Summary_Mor.xlsx"   # adjust path if needed

dat <- read_excel(in_file, sheet = "YearlyAverages") |>
  arrange(Simulation_site, Year) |>
  group_by(Simulation_site) |>
  mutate(WS     = 1 - WR,           # water stress, same year (t)
         WS_lag = lag(WS)) |>       # water stress, previous year (t - 1)
  ungroup()

# Years with both observed and simulated mortality
# (CH_Visp 1999-2020, Spanish sites 2008-2015)
mor <- dat |> filter(!is.na(ObsMor), !is.na(SimMor))

# Site groups: each site, all Spanish sites pooled, all six sites pooled
groups <- bind_rows(
  mor,
  mor |> filter(Simulation_site != "CH_Visp") |>
    mutate(Simulation_site = "Pooled, Spanish sites"),
  mor |> mutate(Simulation_site = "Pooled, all sites")
) |>
  mutate(Simulation_site = factor(Simulation_site,
                                  levels = c("CH_Visp", "SP_Arcalis1", "SP_Arcalis2", "SP_Arcalis3",
                                             "SP_Arcalis4", "SP_Prades",
                                             "Pooled, Spanish sites", "Pooled, all sites")))

# Correlation with p value. Spearman uses exact = FALSE because of tied values
# (e.g. several years with zero observed mortality).
cor_p <- function(x, y, method) {
  ok <- complete.cases(x, y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method,
                                  exact = if (method == "spearman") FALSE else NULL))
  c(est = unname(ct$estimate), p = ct$p.value)
}

# -----------------------------------------------------------------------------
# 1. Table L1 — observed vs. simulated mortality, now with RMSE for all sites
#    (p values are printed for reference only)
# -----------------------------------------------------------------------------
tab_L1 <- groups |>
  group_by(Simulation_site) |>
  summarise(
    Period  = if (first(as.character(Simulation_site)) == "CH_Visp") "1999-2020"
    else if (grepl("Pooled", first(as.character(Simulation_site)))) "-"
    else "2008-2015",
    n       = n(),
    r       = cor_p(ObsMor, SimMor, "pearson")["est"],
    p_r     = cor_p(ObsMor, SimMor, "pearson")["p"],
    rho     = cor_p(ObsMor, SimMor, "spearman")["est"],
    p_rho   = cor_p(ObsMor, SimMor, "spearman")["p"],
    RMSE    = sqrt(mean((ObsMor - SimMor)^2)),
    .groups = "drop")

print(tab_L1 |> mutate(across(c(r, rho, RMSE), \(x) round(x, 2)),
                       across(c(p_r, p_rho), \(x) signif(x, 2))), n = Inf)
write.csv(tab_L1, "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Table_L1_with_RMSE.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. Table L2 — correlations of observed and simulated mortality with
#    water stress (1 - WR), same year (t) and previous year (t - 1);
#    Pearson r and Spearman rho with p values; per site and pooled
# -----------------------------------------------------------------------------
tab_L2 <- groups |>
  group_by(Simulation_site) |>
  group_modify(\(d, key) {
    combos <- expand_grid(Mortality = c("Observed", "Simulated"),
                          Lag       = c("t", "t-1"))
    bind_rows(lapply(seq_len(nrow(combos)), \(i) {
      m  <- if (combos$Mortality[i] == "Observed") d$ObsMor else d$SimMor
      ws <- if (combos$Lag[i] == "t") d$WS else d$WS_lag
      pr <- cor_p(m, ws, "pearson"); sp <- cor_p(m, ws, "spearman")
      tibble(Mortality = combos$Mortality[i], WS_year = combos$Lag[i],
             n = sum(complete.cases(m, ws)),
             r = pr["est"], p_r = pr["p"], rho = sp["est"], p_rho = sp["p"])
    }))
  }) |>
  ungroup()

print(tab_L2 |> mutate(across(c(r, rho), \(x) round(x, 2)),
                       across(c(p_r, p_rho), \(x) signif(x, 2))), n = Inf)
write.csv(tab_L2, "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Table_L2_WS_correlations.csv", row.names = FALSE)

# Wide version laid out like Table L1 (one row per site)
tab_L2_wide <- tab_L2 |>
  mutate(key = paste(Mortality, WS_year, sep = "_")) |>
  select(Simulation_site, key, r, p_r, rho, p_rho) |>
  pivot_wider(names_from = key, values_from = c(r, p_r, rho, p_rho),
              names_glue = "{key}_{.value}")
write.csv(tab_L2_wide, "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Table_L2_WS_correlations_wide.csv", row.names = FALSE)

tab_L2_wide

# -----------------------------------------------------------------------------
# 3. Figure — mortality vs. water stress (1 - WR)
#    Panel (a) observed, panel (b) simulated mortality.
#    Points: each site-year once (Spanish sites vs. CH_Visp); filled = WS in
#            year t, open = WS in year t - 1.
#    Lines : Theil-Sen fits for the Spanish sites pooled and all sites pooled;
#            solid = t, dashed = t - 1.
#    Labels: Spearman rho for each line, taken from Table L2.
# -----------------------------------------------------------------------------

# Theil-Sen line: slope = median of all pairwise slopes,
# intercept = median(y - slope * x)
theil_sen <- function(x, y) {
  ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]
  ij <- combn(length(x), 2)
  dx <- x[ij[2, ]] - x[ij[1, ]]
  dy <- y[ij[2, ]] - y[ij[1, ]]
  b  <- median(dy[dx != 0] / dx[dx != 0])
  a  <- median(y - b * x)
  c(intercept = a, slope = b)
}

lag_lab  <- c(t = "t", t1 = "t - 1")
pool_lab <- c("Spanish sites pooled", "All sites pooled")
mor_lab  <- c(ObsMor = "(a) Observed mortality", SimMor = "(b) Simulated mortality")

# Points: one row per site-year, mortality type and water-stress year
pts <- mor |>
  mutate(Region = ifelse(Simulation_site == "CH_Visp", "CH_Visp", "Spanish sites")) |>
  select(Simulation_site, Year, Region, ObsMor, SimMor, WS, WS_lag) |>
  pivot_longer(c(ObsMor, SimMor), names_to = "Mortality", values_to = "Mor") |>
  pivot_longer(c(WS, WS_lag), names_to = "WS_year", values_to = "WS_value") |>
  filter(!is.na(WS_value)) |>
  mutate(Lag       = factor(ifelse(WS_year == "WS", lag_lab["t"], lag_lab["t1"]),
                            levels = lag_lab),
         Mortality = factor(mor_lab[Mortality], levels = mor_lab))

# Theil-Sen lines for the two pooled groups
ts_lines <- bind_rows(
  pts |> filter(Region == "Spanish sites") |> mutate(Pool = pool_lab[1]),
  pts |> mutate(Pool = pool_lab[2])
) |>
  group_by(Pool, Mortality, Lag) |>
  summarise(intercept = theil_sen(WS_value, Mor)["intercept"],
            slope     = theil_sen(WS_value, Mor)["slope"],
            x_min     = min(WS_value), x_max = max(WS_value),
            .groups   = "drop") |>
  mutate(Pool  = factor(Pool, levels = pool_lab),
         y_min = intercept + slope * x_min,
         y_max = intercept + slope * x_max)
print(ts_lines)

# Spearman rho labels from Table L2 (pooled rows only)
rho_lab <- tab_L2 |>
  filter(grepl("Pooled", Simulation_site)) |>
  mutate(Pool      = factor(ifelse(Simulation_site == "Pooled, Spanish sites",
                                   pool_lab[1], pool_lab[2]), levels = pool_lab),
         Lag       = factor(ifelse(WS_year == "t", lag_lab["t"], lag_lab["t1"]),
                            levels = lag_lab),
         Mortality = factor(ifelse(Mortality == "Observed", mor_lab["ObsMor"],
                                   mor_lab["SimMor"]), levels = mor_lab),
         short     = ifelse(Pool == pool_lab[1], "Spanish", "All sites"),
         label     = sprintf("'%s, %s:'~rho == '%.2f'", short, Lag, rho),
         order     = (as.integer(Pool) - 1) * 2 + as.integer(Lag),
         vj        = 1.4 + 1.3 * (order - 1))              # stack labels top-left

# Colours: points by region, lines by pooled group
col_vals <- c("Spanish sites"        = "#7FA7D6",   # points, Spanish site-years
              "CH_Visp"              = "#6BAA75",   # points, Visp site-years
              "Spanish sites pooled" = "#1F4E9C",   # line + label
              "All sites pooled"     = "#000000")   # line + label

p_ws <- ggplot() +
  geom_point(data = pts,
             aes(x = WS_value, y = Mor, colour = Region, shape = Lag),
             size = 1.9, alpha = 0.8, stroke = 0.7) +
  geom_segment(data = ts_lines,
               aes(x = x_min, xend = x_max, y = y_min, yend = y_max,
                   colour = Pool, linetype = Lag),
               linewidth = 0.9) +
  geom_text(data = rho_lab,
            aes(x = -Inf, y = Inf, label = label, vjust = vj, colour = Pool),
            hjust = -0.06, size = 3, parse = TRUE, show.legend = FALSE) +
  facet_wrap(~ Mortality, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = col_vals,
                      breaks = c("Spanish sites", "CH_Visp",
                                 "Spanish sites pooled", "All sites pooled"),
                      labels = c("Spanish site-years", "CH_Visp site-years",
                                 "Theil-Sen: Spanish sites pooled",
                                 "Theil-Sen: all sites pooled"),
                      name = NULL) +
  scale_shape_manual(values = c(16, 1), name = "Water stress year") +
  scale_linetype_manual(values = c("solid", "22"), name = "Water stress year") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.5))) +   # headroom for labels
  guides(colour = guide_legend(order = 1, ncol = 2,
                               override.aes = list(shape    = c(16, 16, NA, NA),
                                                   linetype = c(0, 0, 1, 1))),
         shape  = guide_legend(order = 2), linetype = guide_legend(order = 2)) +
  labs(x = "Water stress (1 - WR)",
       y = expression(Mortality~("%"~yr^{-1}))) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        strip.text       = element_text(hjust = 0),
        legend.position  = "bottom",
        legend.title     = element_text(size = 9),
        legend.key.width = unit(1.1, "cm"))

print(p_ws)
ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Fig_L_mortality_vs_water_stress_pooled1.png",
       p_ws, width = 9, height = 4.8, dpi = 300)