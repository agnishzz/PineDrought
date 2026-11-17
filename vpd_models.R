library(dplyr)
library(ggstats)
library(glmmTMB)
library(broom.mixed)
library(performance)

vpd_mor <- readxl::read_xlsx("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/Mor_VPD.xlsx")
head(vpd_mor)

# =============================================================================
# VPD single-predictor ordered beta models for observed mortality (ObsMor)
#   Pooled      : both countries, random intercept for Site
#   Spain       : Spanish sites, random intercept for Site
#   Switzerland : CH_Visp only, no random effect
# Same framework as M1a / M2a / M3a: glmmTMB ordbeta (exact zeros kept),
# predictors standardised within each subset, common complete-case rows,
# AICc and LRT against the null, DHARMa diagnostics, R2.
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
library(performance)

# ---- Settings ---------------------------------------------------------------
in_dat   <- vpd_mor
response <- "ObsMor"
out_dir  <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Models/VPD_single_ObsMor"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

preds <- c("YL_VPD_GP_p95", "YL_VPD_SS_p95",
           "AN_VPD_GP_p95", "AN_VPD_SS_p95",
           "VPD_GP_p95",    "VPD_SS_p95")

subsets <- list(
  Pooled      = list(countries = c("Spain", "Switzerland"), re = TRUE),
  Spain       = list(countries = "Spain",                   re = TRUE),
  Switzerland = list(countries = "Switzerland",             re = FALSE)
)

alpha_diag <- 0.05
n_sim      <- 1000

# ---- Helpers ----------------------------------------------------------------
fit_ob <- function(d, vars, re) {
  rhs <- if (length(vars) > 0) paste(vars, collapse = " + ") else "1"
  if (re) rhs <- paste(rhs, "+ (1 | Site)")
  f <- as.formula(paste("y ~", rhs))
  # do.call embeds formula and data in the model call, so performance::r2_nakagawa
  # can refit the null model from inside functions (otherwise R2 is wrong)
  tryCatch(do.call(glmmTMB, list(formula = f, data = d, family = ordbeta())),
           error = function(e) NULL)
}

is_ok <- function(m) {
  !is.null(m) && isTRUE(m$sdr$pdHess)
}

safe_cor2 <- function(a, b) {
  if (length(a) < 3 || sd(a) == 0 || sd(b) == 0) return(NA_real_)
  cor(a, b)^2
}

# R2 measures
#   R2_marginal / R2_conditional : Nakagawa, random-effect models only
#   R2_cor    : squared correlation of observed and fitted (response scale, all years)
#   R2_efron  : 1 - SS_res / SS_tot (response scale, all years)
#   R2_FCN_nonzero : Ferrari & Cribari-Neto (2004), cor(eta, logit(y))^2 on 0 < y < 1
get_r2 <- function(m, d, re) {
  y   <- d$y
  mu  <- as.numeric(fitted(m))
  eta <- as.numeric(predict(m, type = "link"))
  nz  <- y > 0 & y < 1
  
  r2m <- NA_real_
  r2c <- NA_real_
  if (re) {
    r <- tryCatch(suppressWarnings(performance::r2_nakagawa(m)),
                  error = function(e) NULL)
    if (is.list(r)) {
      r2m <- unname(as.numeric(r$R2_marginal))
      r2c <- unname(as.numeric(r$R2_conditional))
    }
  }
  
  tibble(
    R2_marginal    = r2m,
    R2_conditional = r2c,
    R2_cor         = safe_cor2(mu, y),
    R2_efron       = 1 - sum((y - mu)^2) / sum((y - mean(y))^2),
    R2_FCN_nonzero = safe_cor2(eta[nz], qlogis(y[nz])),
    n_nonzero      = sum(nz)
  )
}

# Within-site temporal autocorrelation (smallest p across sites)
autocorr_p <- function(res, d) {
  ps <- vapply(unique(d$Site), function(s) {
    idx <- which(d$Site == s)
    if (length(idx) < 4) return(NA_real_)
    r <- recalculateResiduals(res, sel = d$Site == s)
    tryCatch(testTemporalAutocorrelation(r, time = d$Year[idx], plot = FALSE)$p.value,
             error = function(e) NA_real_)
  }, numeric(1))
  if (all(is.na(ps))) NA_real_ else min(ps, na.rm = TRUE)
}

diagnose <- function(m, d, re) {
  set.seed(123)
  res <- simulateResiduals(m, n = n_sim)
  
  p_unif <- testUniformity(res, plot = FALSE)$p.value
  p_disp <- testDispersion(res, plot = FALSE)$p.value
  p_out  <- testOutliers(res, type = "bootstrap", plot = FALSE)$p.value
  p_zi   <- testZeroInflation(res, plot = FALSE)$p.value
  p_ac   <- if (re) autocorr_p(res, d) else NA_real_
  sing   <- if (re) isTRUE(performance::check_singularity(m)) else NA
  
  pvals <- c(p_unif, p_disp, p_out, p_zi, p_ac)
  pass  <- all(pvals[!is.na(pvals)] > alpha_diag) && !isTRUE(sing)
  
  list(
    tab = tibble(p_uniformity = p_unif, p_dispersion = p_disp,
                 p_outliers = p_out, p_zeroinfl = p_zi,
                 p_autocorr_min_site = p_ac, singular = sing,
                 diag_pass = pass),
    res = res
  )
}

get_coef <- function(m, v) {
  cf <- summary(m)$coefficients$cond
  if (!v %in% rownames(cf)) {
    return(tibble(estimate = NA_real_, se = NA_real_, p_wald = NA_real_))
  }
  tibble(estimate = cf[v, "Estimate"],
         se       = cf[v, "Std. Error"],
         p_wald   = cf[v, "Pr(>|z|)"])
}

# ---- Run one subset ---------------------------------------------------------
run_subset <- function(name, cfg) {
  d <- in_dat %>%
    filter(Country %in% cfg$countries) %>%
    mutate(y = .data[[response]] / 100) %>%
    drop_na(all_of(c("y", preds))) %>%
    mutate(across(all_of(preds), ~ as.numeric(scale(.x))),
           Site = factor(Site))
  
  if (nrow(d) == 0) stop("No rows for subset ", name, ". Check Country labels.")
  
  cat("\n==========", name, "==========\n")
  cat("Rows:", nrow(d), "| Sites:", nlevels(d$Site),
      "| Zeros:", sum(d$y == 0), "| Random effect:", cfg$re, "\n")
  
  m_null <- fit_ob(d, character(0), cfg$re)
  if (!is_ok(m_null)) {
    warning("Null model did not converge for ", name, "; subset skipped.")
    return(tibble(subset = name, predictor = "NULL", converged = FALSE))
  }
  aicc_null <- AICc(m_null)
  
  pdf(file.path(out_dir, paste0("VPD_", name, "_DHARMa.pdf")), width = 10, height = 5)
  on.exit(dev.off(), add = TRUE)
  
  rows <- map(c("NULL", preds), function(v) {
    vars <- if (v == "NULL") character(0) else v
    m    <- if (v == "NULL") m_null else fit_ob(d, vars, cfg$re)
    
    if (!is_ok(m)) {
      return(tibble(subset = name, predictor = v, converged = FALSE))
    }
    
    lrt <- if (v == "NULL") NA_real_ else anova(m_null, m)[["Pr(>Chisq)"]][2]
    dg  <- diagnose(m, d, cfg$re)
    plot(dg$res, title = paste(name, "|", v))
    
    bind_cols(
      tibble(subset = name, predictor = v, converged = TRUE, n = nrow(d),
             AICc = AICc(m), dAICc_null = AICc(m) - aicc_null, LRT_p = lrt),
      if (v == "NULL") tibble(estimate = NA_real_, se = NA_real_, p_wald = NA_real_)
      else get_coef(m, v),
      get_r2(m, d, cfg$re),
      dg$tab
    )
  })
  
  tab <- bind_rows(rows) %>%
    arrange(AICc) %>%
    mutate(dAICc_best = AICc - min(AICc, na.rm = TRUE))
  
  if (!cfg$re) tab <- select(tab, -any_of(c("R2_marginal", "R2_conditional",
                                            "p_autocorr_min_site", "singular")))
  
  print(tab, width = Inf, n = Inf)
  write.csv(tab, file.path(out_dir, paste0("VPD_", name, "_single_predictor.csv")),
            row.names = FALSE)
  tab
}

# ---- Run all three ----------------------------------------------------------
results <- imap(subsets, function(cfg, name) run_subset(name, cfg))

all_tab <- bind_rows(results)
write.csv(all_tab, file.path(out_dir, "VPD_all_subsets_single_predictor.csv"),
          row.names = FALSE)

# Compact summary: one row per subset x predictor
summary_tab <- all_tab %>%
  filter(predictor != "NULL") %>%
  select(any_of(c("subset", "predictor", "converged", "estimate", "LRT_p", "AICc",
                  "dAICc_null", "R2_marginal", "R2_conditional", "R2_cor",
                  "R2_efron", "R2_FCN_nonzero", "diag_pass")))
print(summary_tab, width = Inf, n = Inf)