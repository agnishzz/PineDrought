library(readxl)
library(dplyr)


in_file <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/TreeDensitySimMor.xlsx"   # adjust path if needed

td <- read_excel(in_file) |>
  arrange(Site, Year) |>
  group_by(Site) |>
  mutate(first3 = Year <= min(Year) + 2) |>   # first three simulated years per site
  ungroup()

site_order <- c("CH_Visp", "SP_Arcalis1", "SP_Arcalis2", "SP_Arcalis3",
                "SP_Arcalis4", "SP_Prades")

# Regression + correlations for one data subset.
# With n = 3, the regression has 1 residual df and Spearman rho can only take
# the values -1, -0.5, 0.5 or 1; p values are reported but are of limited use.
fit_one <- function(d) {
  n   <- nrow(d)
  m   <- lm(SimMor ~ SimTD, data = d)
  cf  <- summary(m)$coefficients
  pr  <- cor.test(d$SimTD, d$SimMor, method = "pearson")
  # exact Spearman p value when there are no ties, normal approximation otherwise
  ties <- anyDuplicated(d$SimTD) > 0 || anyDuplicated(d$SimMor) > 0
  sp  <- suppressWarnings(cor.test(d$SimTD, d$SimMor, method = "spearman",
                                   exact = !ties))
  tibble(
    n         = n,
    intercept = cf["(Intercept)", "Estimate"],
    slope     = cf["SimTD", "Estimate"],
    slope_SE  = cf["SimTD", "Std. Error"],
    slope_p   = cf["SimTD", "Pr(>|t|)"],
    R2        = summary(m)$r.squared,
    r         = unname(pr$estimate),
    p_r       = pr$p.value,
    rho       = unname(sp$estimate),
    p_rho     = sp$p.value
  )
}

run_set <- function(d, label) {
  bind_rows(
    fit_one(d) |> mutate(Site = "Pooled, all sites"),
    #fit_one(filter(d, Site != "CH_Visp")) |> mutate(Site = "Pooled, Spanish sites"),
    d |> group_by(Site) |> group_modify(\(x, k) fit_one(x)) |> ungroup()
  ) |>
    mutate(Period = label,
           Site   = factor(Site, levels = c("Pooled, all sites", "Pooled, Spanish sites", site_order))) |>
    arrange(Site) |>
    relocate(Period, Site)
}

tab_L3 <- bind_rows(
  run_set(td,                  "All simulated years"),
  run_set(filter(td, first3),  "First two simulated years")
)

print(tab_L3 |>
        mutate(across(c(intercept, R2, r, rho), \(x) round(x, 2)),
               across(c(slope, slope_SE), \(x) signif(x, 3)),
               across(c(slope_p, p_r, p_rho), \(x) signif(x, 2))),
      n = Inf, width = Inf)

write.csv(tab_L3, "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Tabledensity_mortality2.csv", row.names = FALSE)


d3  <- filter(td, first3) |> mutate(Site = factor(Site, levels = site_order))
fit <- lm(SimMor ~ SimTD, data = d3)
r2  <- summary(fit)$r.squared

p_td <- ggplot(d3, aes(x = SimTD, y = SimMor)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey80", linewidth = 0.7) +
  geom_point(aes(colour = Site, shape = Site), size = 2.8) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, size = 4,
           label = sprintf("R^2 == '%.2f'", r2), parse = TRUE) +
  scale_colour_manual(values = c(CH_Visp = "#24402F", SP_Arcalis1 = "#E08E2D",
                                 SP_Arcalis2 = "#3B7DD8", SP_Arcalis3 = "#AECF7E",
                                 SP_Arcalis4 = "#B5446E", SP_Prades = "#7A5C3E")) +
  scale_shape_manual(values = c(17, 16, 15, 18, 25, 8)) +
  labs(x = "Simulated tree density",            # add unit, e.g. (trees ha-1)
       y = expression(Simulated~mortality~("%"~yr^{-1})),
       colour = NULL, shape = NULL) +
  theme_classic(base_size = 11) +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        legend.position = "right")

print(p_td)
ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/TD_SimMor2.png",
       p_td, width = 6, height = 4, dpi = 300)
###############################################################################################################

td <- td |>
  group_by(Site) |>
  mutate(first2 = Year <= min(Year) + 1) |>   # first two simulated years per site
  ungroup()

d2 <- filter(td, first2) |> mutate(Site = factor(Site, levels = site_order))

tab_first2 <- fit_one(d2) |>
  mutate(Period = "First two simulated years", Site = "Pooled, all sites") |>
  relocate(Period, Site)

print(tab_first2 |>
        mutate(across(c(intercept, R2, r, rho), \(x) round(x, 2)),
               across(c(slope, slope_SE), \(x) signif(x, 3)),
               across(c(slope_p, p_r, p_rho), \(x) signif(x, 2))),
      width = Inf)

# write.csv(tab_first2, "Table_density_mortality_first2yrs.csv", row.names = FALSE)

fit2 <- lm(SimMor ~ SimTD, data = d2)
r2_2 <- summary(fit2)$r.squared

p_td2 <- ggplot(d2, aes(x = SimTD, y = SimMor)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey80", linewidth = 0.7) +
  geom_point(aes(colour = Site, shape = Site), size = 2.8) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, size = 4,
           label = sprintf("R^2 == '%.2f'", r2_2), parse = TRUE) +
  scale_colour_manual(values = c(CH_Visp = "#24402F", SP_Arcalis1 = "#E08E2D",
                                 SP_Arcalis2 = "#3B7DD8", SP_Arcalis3 = "#AECF7E",
                                 SP_Arcalis4 = "#B5446E", SP_Prades = "#7A5C3E")) +
  scale_shape_manual(values = c(17, 16, 15, 18, 25, 8)) +
  labs(x = "Simulated tree density",            # add unit, e.g. (trees ha-1)
       y = expression(Simulated~mortality~("%"~yr^{-1})),
       colour = NULL, shape = NULL) +
  theme_classic(base_size = 11) +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        legend.position = "right")

print(p_td2)

#####################################################################################
#####################################################################################

library(readxl)
library(dplyr)
library(ggplot2)

in_file <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/TreeDensitySimMor.xlsx"   # adjust path if needed

td <- read_excel(in_file) |>
  arrange(Site, Year) |>
  group_by(Site) |>
  mutate(first2 = Year <= min(Year) + 1) |>   # first two simulated years per site
  ungroup()

site_order <- c("CH_Visp", "SP_Arcalis1", "SP_Arcalis2", "SP_Arcalis3",
                "SP_Arcalis4", "SP_Prades")

# Regression + correlations for one data subset.
# With n = 2 (one site, first two years) the line passes exactly through both
# points: there are no residual df, so SE, p, R2 and the correlations are
# undefined (cor.test() errors for n < 3). In that case only the intercept and
# slope of the two-point line are returned; everything else is NA.
fit_one <- function(d) {
  n <- nrow(d)
  m <- lm(SimMor ~ SimTD, data = d)
  
  out <- tibble(
    n         = n,
    intercept = unname(coef(m)["(Intercept)"]),
    slope     = unname(coef(m)["SimTD"]),     # NA if SimTD identical in both years
    slope_SE  = NA_real_,
    slope_p   = NA_real_,
    R2        = NA_real_,
    r         = NA_real_,
    p_r       = NA_real_,
    rho       = NA_real_,
    p_rho     = NA_real_
  )
  if (n < 3) return(out)
  
  cf <- summary(m)$coefficients
  pr <- cor.test(d$SimTD, d$SimMor, method = "pearson")
  # exact Spearman p value when there are no ties, normal approximation otherwise
  ties <- anyDuplicated(d$SimTD) > 0 || anyDuplicated(d$SimMor) > 0
  sp <- suppressWarnings(cor.test(d$SimTD, d$SimMor, method = "spearman",
                                  exact = !ties))
  out |>
    mutate(
      slope_SE = cf["SimTD", "Std. Error"],
      slope_p  = cf["SimTD", "Pr(>|t|)"],
      R2       = summary(m)$r.squared,
      r        = unname(pr$estimate),
      p_r      = pr$p.value,
      rho      = unname(sp$estimate),
      p_rho    = sp$p.value
    )
}

run_set <- function(d, label) {
  bind_rows(
    fit_one(d) |> mutate(Site = "Pooled, all sites"),
    #fit_one(filter(d, Site != "CH_Visp")) |> mutate(Site = "Pooled, Spanish sites"),
    d |> group_by(Site) |> group_modify(\(x, k) fit_one(x)) |> ungroup()
  ) |>
    mutate(Period = label,
           Site   = factor(Site, levels = c("Pooled, all sites", "Pooled, Spanish sites", site_order))) |>
    arrange(Site) |>
    relocate(Period, Site)
}

tab_L2 <- bind_rows(
  run_set(td,                  "All simulated years"),
  run_set(filter(td, first2),  "First two simulated years")
)

print(tab_L2 |>
        mutate(across(c(intercept, R2, r, rho), \(x) round(x, 2)),
               across(c(slope, slope_SE), \(x) signif(x, 3)),
               across(c(slope_p, p_r, p_rho), \(x) signif(x, 2))),
      n = Inf, width = Inf)

write.csv(tab_L2, "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Tabledensity_mortality2.csv", row.names = FALSE)


d2  <- filter(td, first2) |> mutate(Site = factor(Site, levels = site_order))
fit <- lm(SimMor ~ SimTD, data = d2)
r2  <- summary(fit)$r.squared

p_td <- ggplot(d2, aes(x = SimTD, y = SimMor)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey80", linewidth = 0.7) +
  geom_point(aes(colour = Site, shape = Site), size = 2.8) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, size = 4,
           label = sprintf("R^2 == '%.2f'", r2), parse = TRUE) +
  scale_colour_manual(values = c(CH_Visp = "#24402F", SP_Arcalis1 = "#E08E2D",
                                 SP_Arcalis2 = "#3B7DD8", SP_Arcalis3 = "#AECF7E",
                                 SP_Arcalis4 = "#B5446E", SP_Prades = "#7A5C3E")) +
  scale_shape_manual(values = c(17, 16, 15, 18, 25, 8)) +
  labs(x = "Simulated tree density",            # add unit, e.g. (trees ha-1)
       y = expression(Simulated~mortality~("%"~yr^{-1})),
       colour = NULL, shape = NULL) +
  theme_classic(base_size = 11) +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        legend.position = "right")

print(p_td)
ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/TD_SimMor2.png",
       p_td, width = 6, height = 4, dpi = 300)
