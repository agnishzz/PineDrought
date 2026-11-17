# =============================================================================
# Soil water check for the timing of simulated mortality drivers
#   Spain : growing-season SWC vs CURRENT-year summer precipitation anomaly
#   Visp  : growing-season SWC vs PREVIOUS-year summer CWB anomaly,
#           adjusted for the current summer's CWB
# Inputs : AllCLimateAnomalies.xlsx, SoilWaterContent.xlsx (iLand output)
# Method : site fixed effects (Spain), all variables centred within site and
#          z-scored -> standardised coefficients (b per SD)
# =============================================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

f_clim <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/AllCLimateAnomalies.xlsx"
f_swc  <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/SoilWaterContent.xlsx"
swc_var <- "mean_swc_gs_mm"          # growing season; "mean_swc_mm" = whole year

# ---- 1. Data ---------------------------------------------------------------
clim <- read_excel(f_clim) |>
  filter(!is.na(Site)) |>
  mutate(year = as.integer(year)) |>
  select(Site, year, AN_Prec_SS, AN_CWB_SS)

swc <- read_excel(f_swc) |>
  rename(Site = Simulation_site, year = Year) |>
  mutate(year = as.integer(year)) |>
  select(Site, year, SWC = all_of(swc_var))

# SWC joined onto the full climate grid, so lag() always uses the previous
# calendar year (missing SWC years stay NA instead of shifting the series)
d <- clim |>
  left_join(swc, by = c("Site", "year")) |>
  arrange(Site, year) |>
  group_by(Site) |>
  mutate(AN_Prec_SS_L1 = lag(AN_Prec_SS),
         AN_CWB_SS_L1  = lag(AN_CWB_SS)) |>
  ungroup() |>
  mutate(Region = if_else(grepl("^SP", Site), "Spain", "Visp"))

# ---- 2. Model helper ------------------------------------------------------
# centre within site, z-score, fit with Site factor (keeps fixed-effect df)
fit_std <- function(df, y, xs) {
  df <- df |>
    filter(if_all(all_of(c(y, xs)), ~ !is.na(.x))) |>
    group_by(Site) |>
    mutate(across(all_of(c(y, xs)), ~ .x - mean(.x))) |>
    ungroup() |>
    mutate(across(all_of(c(y, xs)), ~ as.numeric(scale(.x))))
  f <- reformulate(c(xs, if (n_distinct(df$Site) > 1) "Site"), y)
  m <- lm(f, data = df)
  ct <- summary(m)$coefficients[xs, , drop = FALSE]
  data.frame(term = xs, b = ct[, 1], se = ct[, 2], p = ct[, 4],
             n = nrow(df), row.names = NULL)
}

spain <- filter(d, Region == "Spain")
visp  <- filter(d, Region == "Visp")

res <- bind_rows(
  # Spain: current-year precipitation (all available years)
  cbind(Region = "Spain", model = "SWC ~ Prec_SS",
        fit_std(spain, "SWC", "AN_Prec_SS")),
  # Spain check: does the previous year add anything?
  cbind(Region = "Spain", model = "SWC ~ Prec_SS + Prec_SS_t-1",
        fit_std(spain, "SWC", c("AN_Prec_SS", "AN_Prec_SS_L1"))),
  # Visp: previous-year CWB, alone and adjusted for current CWB
  cbind(Region = "Visp", model = "SWC ~ CWB_SS_t-1",
        fit_std(visp, "SWC", "AN_CWB_SS_L1")),
  cbind(Region = "Visp", model = "SWC ~ CWB_SS + CWB_SS_t-1",
        fit_std(visp, "SWC", c("AN_CWB_SS", "AN_CWB_SS_L1")))
)
print(res, digits = 3)

# Weather persistence at Visp (is summer CWB autocorrelated?)
pers <- with(filter(visp, !is.na(AN_CWB_SS_L1)), cor.test(AN_CWB_SS, AN_CWB_SS_L1))
cat(sprintf("\nVisp CWB_SS lag-1 autocorrelation: r = %.2f, p = %.2f\n",
            pers$estimate, pers$p.value))

write.csv(res, "SWC_timing_results.csv", row.names = FALSE)

# ---- 3. Figure -------------------------------------------------------------
ink   <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e1e0d9"
col_sp <- "#eb6834"; col_ch <- "#2a78d6"

theme_fig <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
        axis.text  = element_text(colour = ink2),
        axis.title = element_text(colour = ink),
        plot.title = element_text(colour = ink, face = "bold", size = 10),
        plot.subtitle = element_text(colour = ink2, size = 8.5),
        plot.background = element_rect(fill = "white", colour = NA))

lab_stat <- function(r) sprintf("b = %.2f, p = %.3f, n = %d", r$b, r$p, r$n)

# (a) Spain: within-site anomalies, current-year summer precipitation
pa <- spain |>
  filter(!is.na(SWC), !is.na(AN_Prec_SS)) |>
  group_by(Site) |>
  mutate(x = AN_Prec_SS - mean(AN_Prec_SS), y = SWC - mean(SWC)) |>
  ungroup()
ra <- filter(res, model == "SWC ~ Prec_SS")

p1 <- ggplot(pa, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "#c3c2b7", linewidth = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, colour = col_sp,
              fill = col_sp, alpha = 0.15, linewidth = 0.8) +
  geom_point(colour = col_sp, size = 2, alpha = 0.85) +
  labs(title = "(a) Spain: current-year precipitation",
       subtitle = lab_stat(ra),
       x = "MSEA precipitation anomaly, year t (mm)",
       y = "Growing-season SWC anomaly, year t (mm)") +
  theme_fig

# (b) Visp: added-variable plot of previous-year CWB, given current-year CWB
pb <- filter(visp, !is.na(SWC), !is.na(AN_CWB_SS), !is.na(AN_CWB_SS_L1))
pb$x <- resid(lm(AN_CWB_SS_L1 ~ AN_CWB_SS, data = pb))
pb$y <- resid(lm(SWC ~ AN_CWB_SS, data = pb))
rb <- filter(res, model == "SWC ~ CWB_SS + CWB_SS_t-1", term == "AN_CWB_SS_L1")

p2 <- ggplot(pb, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "#c3c2b7", linewidth = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, colour = col_ch,
              fill = col_ch, alpha = 0.15, linewidth = 0.8) +
  geom_point(colour = col_ch, size = 2, alpha = 0.85) +
  labs(title = "(b) Visp: previous-year CWB",
       subtitle = paste0(lab_stat(rb), ", adjusted for year t CWB"),
       x = "MSEA CWB anomaly, year t-1 (mm, adjusted)",
       y = "Growing-season SWC, year t (mm, adjusted)") +
  theme_fig

fig <- p1 + p2
fig
ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Figures/Fig_SWC_timing.png", fig, width = 180, height = 85, units = "mm", dpi = 300)
ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Figures/Fig_SWC_timing.pdf", fig, width = 180, height = 85, units = "mm")
