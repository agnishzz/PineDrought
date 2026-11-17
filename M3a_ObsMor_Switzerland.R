# =============================================================================
# Model 3a: Observed Scots pine mortality ~ climate anomalies
# Switzerland only (CH_Visp), single site: GLM without random effect
# Family: ordered beta (glmmTMB::ordbeta), handles exact 0 and 1 directly
#
# Steps
#   1. Data preparation
#   2. Window selection (GP vs SS) using single-predictor models
#   3. Predictor screening within the chosen window (LRT vs null)
#   4. Global model and all admissible subsets (|r| <= 0.7 constraint)
#   5. Residual diagnostics on candidate models (DHARMa, VIF, singularity)
#   6. Final model selection and export
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(glmmTMB)
  library(MuMIn)
  library(DHARMa)
  library(performance)
})

if (packageVersion("glmmTMB") < "1.1.8") stop("glmmTMB >= 1.1.8 required for ordbeta()")

# ---- 0. Settings -------------------------------------------------------------
in_file   <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/ObsMor_anomaly.xlsx"
out_dir   <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Models/M3a_Obs_Switzerland"
response  <- "ObsMor"
model_tag <- "M3a_Obs_Switzerland_(A)"
use_site_re <- FALSE         # single site: no random effect
site_prefix <- "^CH_"         # CH_Visp only

p_screen   <- 0.10           # LRT threshold for screening
r_max      <- 0.70           # max |r| between predictors in one model
delta_cand <- 4              # dAICc range of candidates sent to diagnostics
max_cand   <- 10             # max number of candidates diagnosed
delta_pref <- 2              # among passing models within this dAICc, prefer fewest terms
alpha_diag <- 0.05           # diagnostic test threshold
n_sim      <- 1000           # DHARMa simulations
obs_per_term <- 10           # cap: one predictor per 10 observations

set.seed(123)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base_vars <- c("VPD", "Prec", "Eto",
               "CWB", "WDcount", "SDII", "CDD") #"MaxTemp", "MinTemp", 
windows   <- c(GP = "_GP", SS = "_SS")   # YL_ = year-lagged prefix
var_name  <- function(base, w) paste0("YL_", base, windows[[w]])
base_of   <- function(v) sub("_(WY|GP|SS)$", "", sub("^YL_", "", v))
all_vars  <- unlist(lapply(names(windows), function(w) var_name(base_vars, w)))

# ---- 1. Data -----------------------------------------------------------------
raw <- read_excel(in_file)

dat <- raw %>%
  filter(grepl(site_prefix, Site)) %>%
  mutate(y = .data[[response]] / 100,
         Site = factor(Site)) %>%
  select(Site, Country, Year, y, all_of(all_vars)) %>%
  filter(if_all(c(y, all_of(all_vars)), ~ !is.na(.x))) %>%
  arrange(Site, Year) %>%
  mutate(across(all_of(all_vars), ~ as.numeric(scale(.x)))) %>%
  as.data.frame()

stopifnot(all(dat$y >= 0 & dat$y <= 1))
stopifnot(nrow(dat) > 0, !any(grepl("^SP_", dat$Site)))

cat("\n==== DATA ====\n")
cat("Observations:", nrow(dat), "\n")
print(dat %>% group_by(Site) %>%
        summarise(n = n(), years = paste(range(Year), collapse = "-"),
                  zeros = sum(y == 0), mean_pct = round(mean(y) * 100, 2),
                  max_pct = round(max(y) * 100, 2), .groups = "drop"))

# ---- Helpers -----------------------------------------------------------------
rand_part <- if (use_site_re) " + (1 | Site)" else ""

make_fml <- function(vars) {
  rhs <- if (length(vars)) paste(vars, collapse = " + ") else "1"
  as.formula(paste0("y ~ ", rhs, rand_part))
}

fit_ob <- function(vars) {
  tryCatch(
    suppressWarnings(glmmTMB(make_fml(vars), data = dat, family = ordbeta())),
    error = function(e) NULL
  )
}

is_ok <- function(m) {
  !is.null(m) && m$fit$convergence == 0 && isTRUE(m$sdr$pdHess)
}

get_aicc <- function(m) if (is_ok(m)) AICc(m) else NA_real_

# R2 for single-site models (no random effect, so Nakagawa R2 does not apply)
#   R2_FCN_nonzero: Ferrari & Cribari-Neto (2004) pseudo-R2, squared correlation
#                   between the linear predictor and logit(y). logit(y) is undefined
#                   at y = 0, so it is computed on non-zero years only (n_nonzero)
#   R2_cor:         squared correlation between observed and fitted values on the
#                   response scale, uses all years including zeros
#   R2_efron:       1 - SS_residual / SS_total on the response scale, all years
get_r2 <- function(m) {
  empty <- tibble(R2_FCN_nonzero = NA_real_, n_nonzero = NA_integer_,
                  R2_cor = NA_real_, R2_efron = NA_real_)
  if (!is_ok(m)) return(empty)
  y   <- dat$y
  mu  <- as.numeric(fitted(m))
  eta <- as.numeric(predict(m, type = "link"))
  nz  <- y > 0 & y < 1

  r2_fcn <- NA_real_
  if (sum(nz) > 2 && sd(eta[nz]) > 0) r2_fcn <- cor(eta[nz], qlogis(y[nz]))^2

  r2_cor <- 0
  if (sd(mu) > 0) r2_cor <- cor(y, mu)^2

  r2_efron <- 1 - sum((y - mu)^2) / sum((y - mean(y))^2)

  tibble(R2_FCN_nonzero = r2_fcn, n_nonzero = sum(nz),
         R2_cor = r2_cor, R2_efron = r2_efron)
}

# Coefficient of one predictor from a fitted model (NA if model failed)
get_est <- function(m, v) {
  if (!is_ok(m)) return(NA_real_)
  unname(fixef(m)$cond[v])
}

# Readable label for a set of terms
label_terms <- function(vars) {
  if (length(vars) == 0) return("(null)")
  paste(vars, collapse = " + ")
}

# ---- 2. Window selection -----------------------------------------------------
m_null <- fit_ob(character(0))
if (!is_ok(m_null)) stop("Null model failed to converge")
aicc_null <- AICc(m_null)

single_grid <- expand_grid(window = names(windows), base = base_vars) %>%
  mutate(var = map2_chr(base, window, var_name))

single_models <- map(single_grid$var, fit_ob)
names(single_models) <- single_grid$var

single_tab <- single_grid %>%
  mutate(converged  = map_lgl(single_models, is_ok),
         AICc       = map_dbl(single_models, get_aicc),
         estimate   = map2_dbl(single_models, var, get_est),
         dAICc_null = AICc - aicc_null)

all_aicc <- c(null = aicc_null, setNames(single_tab$AICc, single_tab$var))
w_raw    <- exp(-0.5 * (all_aicc - min(all_aicc, na.rm = TRUE)))
w_all    <- w_raw / sum(w_raw, na.rm = TRUE)
single_tab$weight <- unname(w_all[single_tab$var])

window_tab <- single_tab %>%
  group_by(window) %>%
  summarise(best_var        = var[which.min(AICc)],
            best_AICc       = min(AICc, na.rm = TRUE),
            best_dAICc_null = best_AICc - aicc_null,
            n_improve_2     = sum(dAICc_null < -2, na.rm = TRUE),
            cum_weight      = sum(weight, na.rm = TRUE),
            n_failed        = sum(!converged),
            .groups = "drop") %>%
  arrange(best_AICc)

cat("\n==== WINDOW SELECTION ====\n")
cat("Null model AICc:", round(aicc_null, 2), "\n")
print(window_tab)

#chosen_window <- window_tab$window[1]
chosen_window <- "GP"   # manual override if justified
cat("Chosen window:", chosen_window, "\n")

# ---- 3. Predictor screening --------------------------------------------------
lrt_p <- function(v) {
  m <- single_models[[v]]
  if (!is_ok(m)) return(NA_real_)
  a <- anova(m_null, m)
  a[["Pr(>Chisq)"]][2]
}

screen_tab <- single_tab %>%
  filter(window == chosen_window) %>%
  mutate(LRT_p    = map_dbl(var, lrt_p),
         retained = !is.na(LRT_p) & LRT_p < p_screen) %>%
  arrange(LRT_p)

cat("\n==== SCREENING (LRT vs null, p <", p_screen, ") ====\n")
print(screen_tab %>% select(var, estimate, AICc, dAICc_null, LRT_p, retained))

sel_vars <- screen_tab$var[screen_tab$retained]
if (length(sel_vars) == 0) {
  stop("No predictor passed screening in window ", chosen_window,
       ". Consider reporting the null result or revisiting p_screen.")
}

# ---- 4. Global model and admissible subsets ----------------------------------
cor_mat <- cor(dat[, sel_vars, drop = FALSE])

cat("\n==== CORRELATION AMONG RETAINED PREDICTORS ====\n")
print(round(cor_mat, 2))
hp <- which(abs(cor_mat) > r_max & upper.tri(cor_mat), arr.ind = TRUE)
if (nrow(hp)) {
  cat("Pairs with |r| >", r_max, "(never combined):\n")
  print(data.frame(var1 = rownames(cor_mat)[hp[, 1]],
                   var2 = colnames(cor_mat)[hp[, 2]],
                   r    = round(cor_mat[hp], 2)))
}

# Global model (all retained predictors), reported for completeness
global <- fit_ob(sel_vars)
cat("\nGlobal model converged cleanly:", is_ok(global), "\n")
if (is_ok(global)) print(summary(global))

max_terms <- min(length(sel_vars), max(1, floor(nrow(dat) / obs_per_term)))

is_admissible <- function(vars) {
  if (length(vars) < 2) return(TRUE)
  cm <- abs(cor_mat[vars, vars])
  if (any(cm[upper.tri(cm)] > r_max)) return(FALSE)
  # CWB = Prec - Eto: never allow all three together
  !all(c("Prec", "Eto", "CWB") %in% base_of(vars))
}

combos <- c(list(character(0)),
            unlist(lapply(seq_len(max_terms),
                          function(k) combn(sel_vars, k, simplify = FALSE)),
                   recursive = FALSE))
combos <- combos[map_lgl(combos, is_admissible)]
cat("\nAdmissible subsets to fit:", length(combos), "(max terms =", max_terms, ")\n")

subset_models <- map(combos, fit_ob)

subset_tab <- tibble(
  model_id  = seq_along(combos),
  terms     = map_chr(combos, label_terms),
  k_terms   = lengths(combos),
  converged = map_lgl(subset_models, is_ok),
  AICc      = map_dbl(subset_models, get_aicc)
) %>%
  filter(converged) %>%
  arrange(AICc) %>%
  mutate(delta  = AICc - min(AICc),
         weight = exp(-0.5 * delta) / sum(exp(-0.5 * delta)))

cat("\n==== MODEL RANKING (top 15) ====\n")
print(head(subset_tab, 15))

# ---- 5. Diagnostics ----------------------------------------------------------
diagnose <- function(m) {
  res <- simulateResiduals(m, n = n_sim, seed = 123)

  p_unif <- testUniformity(res, plot = FALSE)$p.value
  p_disp <- testDispersion(res, plot = FALSE)$p.value
  p_out  <- testOutliers(res, type = "bootstrap", plot = FALSE)$p.value
  p_zero <- testZeroInflation(res, plot = FALSE)$p.value

  k <- length(fixef(m)$cond) - 1
  vif_max <- if (k >= 2) {
    tryCatch(max(check_collinearity(m, component = "conditional")$VIF),
             error = function(e) NA_real_)
  } else NA_real_

  singular <- if (use_site_re) {
    tryCatch(check_singularity(m), error = function(e) NA)
  } else FALSE

  tibble(p_uniformity = p_unif, p_dispersion = p_disp,
         p_outliers = p_out, p_zeroinfl = p_zero,
         VIF_max = vif_max, singular = singular)
}

cand <- subset_tab %>%
  filter(delta <= delta_cand) %>%
  slice_head(n = max_cand)

diag_tab <- cand %>%
  mutate(diag = map(model_id, ~ diagnose(subset_models[[.x]]))) %>%
  unnest(diag) %>%
  mutate(r2 = map(model_id, ~ get_r2(subset_models[[.x]]))) %>%
  unnest(r2) %>%
  mutate(pass = p_uniformity > alpha_diag &
                p_dispersion > alpha_diag &
                p_outliers   > alpha_diag &
                p_zeroinfl   > alpha_diag &
                (is.na(VIF_max) | VIF_max < 5) &
                !(singular %in% TRUE))

cat("\n==== DIAGNOSTICS OF CANDIDATE MODELS ====\n")
print(diag_tab %>% select(terms, k_terms, delta, weight, p_uniformity, p_dispersion,
                          p_outliers, p_zeroinfl, VIF_max, singular, R2_FCN_nonzero, n_nonzero, R2_cor, R2_efron, pass),
      width = Inf)

# ---- 6. Final model ----------------------------------------------------------
# The null model is a benchmark, not a candidate for the final model
passing <- diag_tab %>% filter(pass, k_terms > 0)

if (nrow(passing) == 0) {
  warning("No predictor model passed all diagnostics. Top-ranked predictor model used; inspect before reporting.")
  final_row <- diag_tab %>% filter(k_terms > 0) %>% slice(1)
} else {
  best_pass <- min(passing$AICc)
  final_row <- passing %>%
    filter(AICc - best_pass <= delta_pref) %>%
    arrange(k_terms, AICc) %>%
    slice(1)
}

final      <- subset_models[[final_row$model_id]]
final_vars <- combos[[final_row$model_id]]

cat("\n==== FINAL MODEL ====\n")
cat("Window:", chosen_window, "\nTerms:", final_row$terms, "\n")
null_row <- subset_tab %>% filter(k_terms == 0)
cat("dAICc vs null:", round(final_row$AICc - null_row$AICc, 2),
    "| Akaike weight final:", round(final_row$weight, 3),
    "| null:", round(null_row$weight, 3), "\n")
if (final_row$AICc - null_row$AICc > -2)
  cat("NOTE: final model is within 2 AICc of the null; evidence for climate effects is weak.\n")
print(summary(final))
ci <- tryCatch(confint(final), error = function(e) NULL)
if (!is.null(ci)) print(ci)
r2 <- get_r2(final)
print(r2)

# DHARMa plot for final model
res_final <- simulateResiduals(final, n = n_sim, seed = 123)
pdf(file.path(out_dir, paste0(model_tag, "_DHARMa_final.pdf")), width = 10, height = 5)
plot(res_final)
dev.off()

# ---- Export ------------------------------------------------------------------
write.csv(single_tab,  file.path(out_dir, paste0(model_tag, "_01_single_models.csv")), row.names = FALSE)
write.csv(window_tab,  file.path(out_dir, paste0(model_tag, "_02_window_selection.csv")), row.names = FALSE)
write.csv(screen_tab,  file.path(out_dir, paste0(model_tag, "_03_screening.csv")), row.names = FALSE)
write.csv(round(cor_mat, 3), file.path(out_dir, paste0(model_tag, "_04_correlation.csv")))
write.csv(subset_tab,  file.path(out_dir, paste0(model_tag, "_05_model_ranking.csv")), row.names = FALSE)
write.csv(diag_tab,    file.path(out_dir, paste0(model_tag, "_06_diagnostics.csv")), row.names = FALSE)
write.csv(diag_tab %>% filter(pass) %>%
            select(model_id, terms, k_terms, AICc, delta, weight,
                   R2_FCN_nonzero, n_nonzero, R2_cor, R2_efron, everything()),
          file.path(out_dir, paste0(model_tag, "_06b_passing_models_R2.csv")), row.names = FALSE)



library(broom.mixed)

# ---- 6c. Coefficients and 95% CIs for the exported (passing) models ----------
# One row per term per model. Predictors are standardised, so each estimate is
# the change in the logit of mortality per 1 SD change in that predictor.
coef_one <- function(id) {
  m <- subset_models[[id]]
  if (!is_ok(m)) return(NULL)
  tt <- tryCatch(
    broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE, conf.level = 0.95),
    error = function(e) NULL
  )
  if (is.null(tt)) return(NULL)
  if ("component" %in% names(tt)) tt <- tt[tt$component == "cond", , drop = FALSE]
  
  # 90% intervals, matching the p < 0.1 marker used at the screening stage
  tt90 <- tryCatch(
    broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE, conf.level = 0.90),
    error = function(e) NULL
  )
  tt$conf.low90  <- NA_real_
  tt$conf.high90 <- NA_real_
  if (!is.null(tt90)) {
    if ("component" %in% names(tt90)) tt90 <- tt90[tt90$component == "cond", , drop = FALSE]
    idx <- match(tt$term, tt90$term)
    tt$conf.low90  <- tt90$conf.low[idx]
    tt$conf.high90 <- tt90$conf.high[idx]
  }
  vars <- combos[[id]]
  lab  <- "(null)"
  if (length(vars) > 0) lab <- paste(vars, collapse = " + ")
  tt$model_id <- id
  tt$terms    <- lab
  tt$k_terms  <- length(vars)
  tt
}

coef_ids <- diag_tab$model_id[diag_tab$pass & diag_tab$k_terms > 0]
coef_tab <- bind_rows(map(coef_ids, coef_one))

if (nrow(coef_tab) > 0) {
  coef_tab <- coef_tab %>%
    mutate(is_final     = model_id == final_row$model_id,
           odds_ratio   = exp(estimate),
           OR_conf.low  = exp(conf.low),
           OR_conf.high = exp(conf.high),
           sig_95       = conf.low > 0 | conf.high < 0,
           sig_90       = conf.low90 > 0 | conf.high90 < 0) %>%
    select(model_id, terms, k_terms, is_final, term, estimate, std.error,
           conf.low, conf.high, conf.low90, conf.high90, statistic, p.value,
           odds_ratio, OR_conf.low, OR_conf.high, sig_95, sig_90) %>%
    arrange(desc(is_final), k_terms, model_id, term)
  
  cat("\n==== COEFFICIENTS AND 95% CIs (passing models with >= 1 predictor) ====\n")
  print(coef_tab, n = Inf, width = Inf)
  
  #write.csv(coef_tab,
            #file.path(out_dir, paste0(model_tag, "_06c_coefficients_CI.csv")),
            #row.names = FALSE)
#} else {
  #cat("\nNo passing model with predictors, so no coefficient table was written.\n")
}



capture.output(
  { cat("Window:", chosen_window, "\nTerms:", final_row$terms,
        "\ndAICc vs null:", round(final_row$AICc - null_row$AICc, 2), "\n\n")
    print(summary(final)); if (!is.null(ci)) print(ci); print(r2) },
  file = file.path(out_dir, paste0(model_tag, "_07_final_model.txt"))
)

# Structure saved so the simulated-mortality script can refit the same model
saveRDS(list(response = response, window = chosen_window,
             vars = final_vars, formula = make_fml(final_vars),
             site_re = use_site_re),
        file.path(out_dir, paste0(model_tag, "_final_structure.rds")))

cat("\nOutputs written to:", out_dir, "\n")















# ---- 8. Leave-one-year-out cross-validation (CH_Visp) ------------------
# Cross-validation of every passing model with at least one predictor.
# Each fold refits the model on the training rows and compares the fold
# estimates with the full-data estimates, so the output shows how far the
# coefficients move when one year at CH_Visp is removed.
# Predictors keep the full-data standardisation, so all folds are on one scale.

cv_fold_var <- dat$Year
cv_folds    <- sort(unique(dat$Year))

# Refit on an arbitrary subset of the data (same family and random structure)
fit_subset <- function(vars, d) {
  rhs <- "1"
  if (length(vars) > 0) rhs <- paste(vars, collapse = " + ")
  f <- as.formula(paste0("y ~ ", rhs, rand_part))
  tryCatch(suppressWarnings(glmmTMB(f, data = d, family = ordbeta())),
           error = function(e) NULL)
}

predict_holdout <- function(m, test) {
  args <- list(m, newdata = test, type = "response", allow.new.levels = TRUE)
  if (use_site_re) args$re.form <- NA          # population level: held-out site is new
  tryCatch(as.numeric(do.call(predict, args)), error = function(e) rep(NA_real_, nrow(test)))
}

cv_ids <- diag_tab$model_id[diag_tab$pass & diag_tab$k_terms > 0]

if (length(cv_ids) == 0) {
  cat("\nLEAVE-ONE-YEAR-OUT CROSS-VALIDATION (CH_VISP): no passing model with predictors, so cross-validation was skipped.\n")
} else {
  
  # Full-data estimates for reference
  cv_full <- bind_rows(lapply(cv_ids, function(id) {
    e <- fixef(subset_models[[id]])$cond
    tibble(model_id = id, term = names(e), full_estimate = as.numeric(e))
  }))
  
  # One fold at a time, for every candidate model
  cv_raw <- bind_rows(lapply(cv_ids, function(id) {
    vars <- combos[[id]]
    lab  <- paste(vars, collapse = " + ")
    bind_rows(lapply(cv_folds, function(f) {
      train <- droplevels(dat[cv_fold_var != f, , drop = FALSE])
      test  <- dat[cv_fold_var == f, , drop = FALSE]
      m     <- fit_subset(vars, train)
      if (!is_ok(m)) {
        return(tibble(model_id = id, terms = lab, fold = as.character(f),
                      term = NA_character_, estimate = NA_real_,
                      converged = FALSE, n_train = nrow(train),
                      n_holdout = nrow(test), holdout_rmse = NA_real_))
      }
      pred <- predict_holdout(m, test)
      e    <- fixef(m)$cond
      tibble(model_id = id, terms = lab, fold = as.character(f),
             term = names(e), estimate = as.numeric(e),
             converged = TRUE, n_train = nrow(train), n_holdout = nrow(test),
             holdout_rmse = sqrt(mean((test$y - pred)^2, na.rm = TRUE)))
    }))
  }))
  
  # ---- Coefficient stability ----
  cv_coef <- cv_raw %>%
    filter(converged) %>%
    left_join(cv_full, by = c("model_id", "term")) %>%
    group_by(model_id, terms, term) %>%
    summarise(full_estimate  = first(full_estimate),
              n_folds_ok     = n(),
              mean_estimate  = mean(estimate),
              sd_estimate    = sd(estimate),
              min_estimate   = min(estimate),
              max_estimate   = max(estimate),
              prop_same_sign = mean(sign(estimate) == sign(full_estimate)),
              max_abs_change = max(abs(estimate - full_estimate)),
              max_rel_change = max(abs(estimate - full_estimate)) /
                max(abs(first(full_estimate)), 1e-8),
              most_infl_fold = fold[which.max(abs(estimate - full_estimate))],
              .groups = "drop") %>%
    mutate(is_intercept = term == "(Intercept)",
           is_final     = model_id == final_row$model_id,
           stable       = prop_same_sign == 1 & max_rel_change < 0.5) %>%
    arrange(desc(is_final), model_id, is_intercept, term)
  
  # ---- Per-model summary (predictor terms only) and hold-out error ----
  cv_model <- cv_coef %>%
    filter(!is_intercept) %>%
    group_by(model_id, terms) %>%
    summarise(min_prop_same_sign = min(prop_same_sign),
              max_rel_change     = max(max_rel_change),
              all_terms_stable   = all(stable),
              .groups = "drop") %>%
    left_join(
      cv_raw %>%
        filter(converged) %>%
        distinct(model_id, fold, holdout_rmse, n_holdout) %>%
        group_by(model_id) %>%
        summarise(n_folds_ok        = n(),
                  mean_holdout_rmse = mean(holdout_rmse, na.rm = TRUE),
                  max_holdout_rmse  = max(holdout_rmse, na.rm = TRUE),
                  worst_fold        = fold[which.max(holdout_rmse)],
                  .groups = "drop"),
      by = "model_id") %>%
    mutate(is_final = model_id == final_row$model_id) %>%
    arrange(desc(is_final), mean_holdout_rmse)
  
  cv_failed <- cv_raw %>% filter(!converged) %>% distinct(model_id, terms, fold)
  
  cat("\n==== LEAVE-ONE-YEAR-OUT CROSS-VALIDATION (CH_VISP) ====\n")
  cat("Folds:", length(cv_folds), "| each fold holds out one year at CH_Visp\n")
  cat("Models cross-validated:", length(cv_ids), "\n")
  if (nrow(cv_failed) > 0) {
    cat("Folds that failed to converge:", nrow(cv_failed), "\n")
    print(cv_failed, n = Inf)
  }
  
  cat("\n-- Coefficient stability (one row per term per model) --\n")
  print(cv_coef %>% select(model_id, terms, term, is_final, full_estimate,
                           mean_estimate, sd_estimate, min_estimate, max_estimate,
                           prop_same_sign, max_abs_change, max_rel_change,
                           most_infl_fold, stable),
        n = Inf, width = Inf)
  
  cat("\n-- Per-model stability and hold-out error (proportion scale) --\n")
  print(cv_model, n = Inf, width = Inf)
  
  write.csv(cv_raw,   file.path(out_dir, paste0(model_tag, "_09_LOYO_folds.csv")), row.names = FALSE)
  write.csv(cv_coef,  file.path(out_dir, paste0(model_tag, "_09_LOYO_coef_stability.csv")), row.names = FALSE)
  write.csv(cv_model, file.path(out_dir, paste0(model_tag, "_09_LOYO_model_summary.csv")), row.names = FALSE)
  
  cat("\nCross-validation outputs written to:", out_dir, "\n")
}
