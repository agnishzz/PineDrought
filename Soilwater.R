library (readxl)

data <- readxl::read_xlsx(in_file   <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/ObsMor_anomaly.xlsx")
head(data)

library(tidyverse)
library(lme4)
library(lmerTest)
library(MuMIn)
library(broom.mixed)

# ---------------------------------------------------------------
# 1. Align SWC with the lagged climate predictors
#    YL_* at row Year = t describe climate of year t-1, so the
#    matched mediator is SWC of year t-1, not SWC of year t.
#    Join on Year + 1 rather than dplyr::lag() so year gaps are safe.
# ---------------------------------------------------------------
swc_lagged <- data %>%
  transmute(Site, Year = Year + 1, SWC_GP_YL = SWC_GP_mm)

d <- data %>%
  left_join(swc_lagged, by = c("Site", "Year")) %>%
  arrange(Country, Site, Year)

# quick check on what the lag cost you
d %>%
  group_by(Country, Site) %>%
  summarise(n_years = n(), n_matched = sum(!is.na(SWC_GP_YL)), .groups = "drop") %>%
  print(n = Inf)

# ---------------------------------------------------------------
# 2. Standardise within country so betas are comparable across
#    countries despite different units, variances and soil depths
# ---------------------------------------------------------------
z_vars <- c("SWC_GP_YL", "SWC_GP_mm",
            "YL_WDcount_GP", "YL_SDII_GP", "YL_CWB_GP", "YL_Prec_GP")

d_std <- d %>%
  group_by(Country) %>%
  mutate(across(all_of(z_vars), ~ as.numeric(scale(.x)), .names = "z_{.col}")) %>%
  ungroup()

es <- d_std %>% filter(Country == "Spain",       !is.na(z_SWC_GP_YL))
ch <- d_std %>% filter(Country == "Switzerland", !is.na(z_SWC_GP_YL))
conflicted::conflicts_prefer(lmerTest::lmer)
# ---------------------------------------------------------------
# 3. Mediation models: own-country predictors
# ---------------------------------------------------------------
m_es_own <- lmer(z_SWC_GP_YL ~ z_YL_WDcount_GP + z_YL_SDII_GP + (1 | Site),
                 data = es, REML = FALSE)

m_ch_own <- lm(z_SWC_GP_YL ~ z_YL_CWB_GP, data = ch)

# ---------------------------------------------------------------
# 4. Crossover models: each country fitted with the OTHER country's
#    mortality predictors. This is the test that makes the soil
#    attribution non-speculative.
# ---------------------------------------------------------------
m_es_cross <- lmer(z_SWC_GP_YL ~ z_YL_CWB_GP + (1 | Site),
                   data = es, REML = FALSE)

m_ch_cross <- lm(z_SWC_GP_YL ~ z_YL_WDcount_GP + z_YL_SDII_GP, data = ch)

# ---------------------------------------------------------------
# 5. Comparison table
# ---------------------------------------------------------------
fit_row <- function(model, country, label) {
  if (inherits(model, "merMod")) {
    r2 <- MuMIn::r.squaredGLMM(model)
    tibble(Country = country, Predictors = label,
           n = nobs(model),
           R2_marginal = unname(r2[1, "R2m"]),
           R2_conditional = unname(r2[1, "R2c"]),
           AICc = MuMIn::AICc(model),
           singular = lme4::isSingular(model))
  } else {
    tibble(Country = country, Predictors = label,
           n = nobs(model),
           R2_marginal = summary(model)$r.squared,
           R2_conditional = NA_real_,
           AICc = MuMIn::AICc(model),
           singular = NA)
  }
}

comparison <- bind_rows(
  fit_row(m_es_own,   "Spain",       "WDcount_GP + SDII_GP (own)"),
  fit_row(m_es_cross, "Spain",       "CWB_GP (cross)"),
  fit_row(m_ch_own,   "Switzerland", "CWB_GP (own)"),
  fit_row(m_ch_cross, "Switzerland", "WDcount_GP + SDII_GP (cross)")
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(comparison, n = Inf)

# ---------------------------------------------------------------
# 6. Coefficients
# ---------------------------------------------------------------
coefs <- bind_rows(
  broom.mixed::tidy(m_es_own,   effects = "fixed") %>% mutate(Country = "Spain",       Model = "own"),
  broom.mixed::tidy(m_es_cross, effects = "fixed") %>% mutate(Country = "Spain",       Model = "cross"),
  broom::tidy(m_ch_own)                            %>% mutate(Country = "Switzerland", Model = "own"),
  broom::tidy(m_ch_cross)                          %>% mutate(Country = "Switzerland", Model = "cross")
) %>%
  filter(term != "(Intercept)") %>%
  select(Country, Model, term, estimate, std.error, statistic, p.value) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(coefs, n = Inf)

# ---------------------------------------------------------------
# 7. Secondary: legacy chain, climate(t-1) -> SWC(t)
#    Cheap to run, uses the unlagged SWC and the full n
# ---------------------------------------------------------------
m_es_legacy <- lmer(z_SWC_GP_mm ~ z_YL_WDcount_GP + z_YL_SDII_GP + (1 | Site),
                    data = d_std %>% filter(Country == "Spain"), REML = FALSE)
m_ch_legacy <- lm(z_SWC_GP_mm ~ z_YL_CWB_GP,
                  data = d_std %>% filter(Country == "Switzerland"))

summary(m_es_legacy)
summary(m_ch_legacy)
