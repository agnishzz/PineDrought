# =============================================================================
# VPD p95: single-predictor screen under the Tweedie scheme
#
# Same layout as the ordered-beta version: one row per subset x predictor,
# ranked by AICc, with coefficient, LRT against the null, R2 set and
# simulation-based diagnostics. Differences forced by the family:
#   family      ordbeta(logit) -> tweedie(log), power estimated
#   response    ObsMor / 100   -> ObsMor as given (%), zeros kept as zeros
#   R2          Ferrari & Cribari-Neto (beta only) -> Nakagawa for Tweedie
#               (lognormal approximation) plus Nagelkerke vs the null
# Diagnostics use DHARMa when installed, otherwise an equivalent built-in
# simulation (same quantities, same p-values).
# =============================================================================

suppressMessages({
  library(glmmTMB); library(dplyr); library(tidyr); library(tibble)
  library(purrr); library(performance)
})

# ---- Settings ---------------------------------------------------------------
# in_dat: a data frame already in the workspace, or a path to the xlsx
in_dat   <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/Mor_VPD.xlsx"
response <- "ObsMor"
out_dir  <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Models/VPD_single_ObsMor_tweedie"
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
use_dharma <- requireNamespace("DHARMa", quietly = TRUE)

if (is.character(in_dat)) in_dat <- readxl::read_xlsx(in_dat)
stopifnot(is.data.frame(in_dat))
miss <- setdiff(c("Site", "Country", "Year", response, preds), names(in_dat))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

# ---- Helpers ----------------------------------------------------------------
# Tweedie power is estimated per model, unless the null model puts it at a
# boundary (p -> 2 with no zeros, p -> 1), where it is fixed for the subset:
# optimisation at the boundary is slow and unstable, the estimates are the same.
power_fix <- function(d, re) {
  rhs <- if (re) "1 + (1 | Site)" else "1"
  m <- tryCatch(suppressWarnings(glmmTMB(as.formula(paste("y ~", rhs)),
                                         family = tweedie(link = "log"), data = d,
                                         control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  p0 <- unname(family_params(m))
  if (p0 > 1.99) 1.99 else if (p0 < 1.01) 1.01 else NULL
}

fit_tw <- function(d, vars, re, pfix = NULL) {
  rhs <- if (length(vars) > 0) paste(vars, collapse = " + ") else "1"
  if (re) rhs <- paste(rhs, "+ (1 | Site)")
  f <- as.formula(paste("y ~", rhs))
  args <- list(formula = f, data = d, family = tweedie(link = "log"))
  if (!is.null(pfix)) {
    args$start <- list(psi = qlogis(pfix - 1)); args$map <- list(psi = factor(NA))
  }
  # do.call embeds formula and data in the model call, so downstream refits of
  # the null model (e.g. performance) work from inside functions
  m <- tryCatch(suppressWarnings(do.call(glmmTMB, args)), error = function(e) NULL)
  if (!is.null(m) && !isTRUE(m$sdr$pdHess) && is.null(pfix)) {
    args$control <- glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
    m2 <- tryCatch(suppressWarnings(do.call(glmmTMB, args)), error = function(e) NULL)
    if (!is.null(m2) && isTRUE(m2$sdr$pdHess)) m <- m2
  }
  m
}

is_ok <- function(m) {
  !is.null(m) && isTRUE(m$sdr$pdHess) && all(is.finite(sqrt(diag(vcov(m)$cond))))
}

aicc <- function(m) {
  ll <- logLik(m); k <- attr(ll, "df"); n <- nobs(m)
  -2 * as.numeric(ll) + 2 * k + 2 * k * (k + 1) / max(n - k - 1, 1)
}

safe_cor2 <- function(a, b) {
  if (length(a) < 3 || sd(a) == 0 || sd(b) == 0) return(NA_real_)
  cor(a, b)^2
}

# R2 measures
#   R2_marginal / R2_conditional : Nakagawa et al. (2017) for Tweedie (log link),
#     lognormal approximation sigma2_d = log(1 + phi * lambda^(p - 2)).
#     Without a random effect the two are identical.
#   R2_cor        : squared correlation of observed and fitted (response scale)
#   R2_efron      : 1 - SS_res / SS_tot (response scale)
#   R2_Nagelkerke : likelihood-ratio R2 against the intercept-only model (+ RE)
get_r2 <- function(m, m0, d, re) {
  y <- d$y; n <- length(y)
  mu <- as.numeric(predict(m, type = "response"))
  b  <- fixef(m)$cond; X <- getME(m, "X")
  s2f <- var(as.vector(X %*% b))
  s2r <- if (length(VarCorr(m)$cond)) sum(sapply(VarCorr(m)$cond, function(v) v[1])) else 0
  lambda <- unname(exp(b[1] + 0.5 * s2r))
  s2d <- log(1 + sigma(m) * lambda^(unname(family_params(m)) - 2))
  ll <- as.numeric(logLik(m)); ll0 <- as.numeric(logLik(m0))
  cs <- 1 - exp(-2 / n * (ll - ll0))
  tibble(
    R2_marginal    = s2f / (s2f + s2r + s2d),
    R2_conditional = (s2f + s2r) / (s2f + s2r + s2d),
    R2_Nagelkerke  = cs / (1 - exp(2 / n * ll0)),
    R2_cor         = safe_cor2(mu, y),
    R2_efron       = 1 - sum((y - mu)^2) / sum((y - mean(y))^2),
    n_zero         = sum(y == 0)
  )
}

# Within-site temporal autocorrelation (smallest p across sites)
autocorr_p_dharma <- function(res, d) {
  ps <- vapply(unique(d$Site), function(s) {
    idx <- which(d$Site == s)
    if (length(idx) < 4) return(NA_real_)
    r <- DHARMa::recalculateResiduals(res, sel = d$Site == s)
    tryCatch(DHARMa::testTemporalAutocorrelation(r, time = d$Year[idx], plot = FALSE)$p.value,
             error = function(e) NA_real_)
  }, numeric(1))
  if (all(is.na(ps))) NA_real_ else min(ps, na.rm = TRUE)
}

# Built-in equivalent of the DHARMa tests, used when DHARMa is not installed
diag_builtin <- function(m, d, nsim = n_sim) {
  y <- d$y; n <- length(y)
  sims <- as.matrix(simulate(m, nsim = nsim, seed = 123))
  fit  <- as.numeric(predict(m, type = "response"))
  res  <- vapply(seq_len(n), function(i)
    mean(sims[i, ] < y[i]) + runif(1) * mean(sims[i, ] == y[i]), numeric(1))
  p_unif <- suppressWarnings(ks.test(res, "punif"))$p.value
  obs_sd <- sd(y - fit); sim_sd <- apply(sims, 2, function(s) sd(s - fit))
  p_disp <- min(1, 2 * min(mean(sim_sd >= obs_sd), mean(sim_sd <= obs_sd)))
  z_obs <- sum(y == 0); z_sim <- colSums(sims == 0)
  p_zi <- if (z_obs > 0 || mean(z_sim) > 0)
    min(1, 2 * min(mean(z_sim >= z_obs), mean(z_sim <= z_obs))) else NA_real_
  n_out <- sum(res < 1 / (nsim + 1) | res > nsim / (nsim + 1))
  p_out <- binom.test(n_out, n, p = 2 / (nsim + 1))$p.value
  zr <- qnorm(pmin(pmax(res, 1e-4), 1 - 1e-4))
  lag1 <- function(zz) {
    num <- 0; den <- 0
    for (s in unique(d$Site)) {
      idx <- which(d$Site == s); idx <- idx[order(d$Year[idx])]
      v <- zz[idx] - mean(zz[idx])
      if (length(v) > 2) { num <- num + sum(v[-1] * v[-length(v)]); den <- den + sum(v^2) }
    }
    if (den == 0) NA_real_ else num / den
  }
  ac <- lag1(zr)
  p_ac <- if (is.na(ac)) NA_real_ else {
    perm <- replicate(999, { zp <- zr
    for (s in unique(d$Site)) { idx <- which(d$Site == s); zp[idx] <- sample(zr[idx]) }
    lag1(zp) })
    min(1, 2 * min(mean(perm >= ac), mean(perm <= ac)))
  }
  list(p_unif = p_unif, p_disp = p_disp, p_out = p_out, p_zi = p_zi, p_ac = p_ac, res = NULL)
}

diagnose <- function(m, d, re) {
  set.seed(123)
  if (use_dharma) {
    res <- DHARMa::simulateResiduals(m, n = n_sim)
    p_unif <- DHARMa::testUniformity(res, plot = FALSE)$p.value
    p_disp <- DHARMa::testDispersion(res, plot = FALSE)$p.value
    p_out  <- DHARMa::testOutliers(res, type = "bootstrap", plot = FALSE)$p.value
    p_zi   <- DHARMa::testZeroInflation(res, plot = FALSE)$p.value
    p_ac   <- autocorr_p_dharma(res, d)
  } else {
    b <- diag_builtin(m, d)
    p_unif <- b$p_unif; p_disp <- b$p_disp; p_out <- b$p_out; p_zi <- b$p_zi
    p_ac <- b$p_ac; res <- NULL
  }
  sing <- if (re) isTRUE(performance::check_singularity(m)) else NA
  pvals <- c(p_unif, p_disp, p_out, p_zi, p_ac)
  list(
    tab = tibble(p_uniformity = p_unif, p_dispersion = p_disp,
                 p_outliers = p_out, p_zeroinfl = p_zi,
                 p_autocorr_min_site = p_ac, singular = sing,
                 diag_pass = all(pvals[!is.na(pvals)] > alpha_diag) && !isTRUE(sing)),
    res = res
  )
}

get_coef <- function(m, v) {
  cf <- summary(m)$coefficients$cond
  if (!v %in% rownames(cf)) {
    return(tibble(estimate = NA_real_, se = NA_real_, p_wald = NA_real_, ratio_per_SD = NA_real_))
  }
  tibble(estimate = cf[v, "Estimate"], se = cf[v, "Std. Error"],
         p_wald = cf[v, "Pr(>|z|)"], ratio_per_SD = exp(cf[v, "Estimate"]))
}

# ---- Run one subset ---------------------------------------------------------
run_subset <- function(name, cfg) {
  d <- in_dat %>%
    filter(Country %in% cfg$countries) %>%
    mutate(y = .data[[response]]) %>%          # Tweedie: mortality in %, zeros kept
    drop_na(all_of(c("y", preds))) %>%
    mutate(across(all_of(preds), ~ as.numeric(scale(.x))),
           Site = factor(Site))
  
  if (nrow(d) == 0) stop("No rows for subset ", name, ". Check Country labels.")
  
  pfix <- power_fix(d, cfg$re)
  cat("\n==========", name, "==========\n")
  cat("Rows:", nrow(d), "| Sites:", nlevels(d$Site),
      "| Zeros:", sum(d$y == 0), "| Random effect:", cfg$re,
      "| Tweedie power:", if (is.null(pfix)) "estimated" else paste("fixed at", pfix), "\n")
  
  m_null <- fit_tw(d, character(0), cfg$re, pfix)
  if (!is_ok(m_null)) {
    warning("Null model did not converge for ", name, "; subset skipped.")
    return(tibble(subset = name, predictor = "NULL", converged = FALSE))
  }
  aicc_null <- aicc(m_null)
  
  if (use_dharma) {
    pdf(file.path(out_dir, paste0("VPD_", name, "_DHARMa.pdf")), width = 10, height = 5)
    on.exit(dev.off(), add = TRUE)
  }
  
  rows <- map(c("NULL", preds), function(v) {
    vars <- if (v == "NULL") character(0) else v
    m    <- if (v == "NULL") m_null else fit_tw(d, vars, cfg$re, pfix)
    
    if (!is_ok(m)) return(tibble(subset = name, predictor = v, converged = FALSE))
    
    lrt <- if (v == "NULL") NA_real_ else
      pchisq(2 * (as.numeric(logLik(m)) - as.numeric(logLik(m_null))),
             df = attr(logLik(m), "df") - attr(logLik(m_null), "df"), lower.tail = FALSE)
    dg <- diagnose(m, d, cfg$re)
    if (use_dharma && !is.null(dg$res)) plot(dg$res, title = paste(name, "|", v))
    
    bind_cols(
      tibble(subset = name, predictor = v, converged = TRUE, n = nrow(d),
             tweedie_power = unname(family_params(m)),
             AICc = aicc(m), dAICc_null = aicc(m) - aicc_null, LRT_p = lrt),
      if (v == "NULL") tibble(estimate = NA_real_, se = NA_real_, p_wald = NA_real_,
                              ratio_per_SD = NA_real_) else get_coef(m, v),
      get_r2(m, m_null, d, cfg$re),
      dg$tab
    )
  })
  
  tab <- bind_rows(rows) %>%
    arrange(AICc) %>%
    mutate(dAICc_best = AICc - min(AICc, na.rm = TRUE))
  
  if (!cfg$re) tab <- select(tab, -any_of(c("p_autocorr_min_site", "singular")))
  
  print(tab, width = Inf, n = Inf)
  write.csv(tab, file.path(out_dir, paste0("VPD_", name, "_single_predictor_tweedie.csv")),
            row.names = FALSE)
  tab
}

# ---- Run all three ----------------------------------------------------------
results <- imap(subsets, function(cfg, name) run_subset(name, cfg))

all_tab <- bind_rows(results)
write.csv(all_tab, file.path(out_dir, "VPD_all_subsets_single_predictor_tweedie.csv"),
          row.names = FALSE)

# Compact summary: one row per subset x predictor
summary_tab <- all_tab %>%
  filter(predictor != "NULL") %>%
  select(any_of(c("subset", "predictor", "converged", "estimate", "se", "p_wald",
                  "ratio_per_SD", "LRT_p", "AICc", "dAICc_null", "R2_marginal",
                  "R2_conditional", "R2_Nagelkerke", "R2_cor", "R2_efron", "diag_pass")))
print(summary_tab, width = Inf, n = Inf)