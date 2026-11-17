# =============================================================================
# Scots pine mortality drivers: Switzerland SimMor
# Tweedie GLM(M), log link. One of six stand-alone scripts:
#   Pooled / Spain / Switzerland  x  ObsMor / SimMor
#
# Steps
#   0. Timing selection: previous-year (YL_) vs current-year (AN_) anomalies.
#      Univariate models for 7 predictors x 2 windows x 2 timings (+ null).
#      Evidence per timing = summed Akaike weight of its 14 models.
#   1. Window selection (GP, SS) within each selected timing: evidence per
#      window = summed Akaike weight of that timing's 7 models per window.
#      Tie rule (both steps): the best option is selected; another option is
#      ALSO selected when the evidence ratio best/other < TIE_RATIO. Every
#      selected timing x window combination gets its own model.
#   2. Univariate screening in each selected timing x window (p < 0.05; falls back to
#      p < 0.10, then to the single best-AICc predictor, flagged).
#   3. Collinearity filter: greedy by AICc, drop predictor if |r| > 0.7 with
#      an already retained predictor.
#   4. All subsets of retained predictors (max KMAX terms), ranked by AICc.
#      Final = lowest-AICc converged model that passes diagnostics
#      (KS p > 0.05, dispersion p > 0.05, VIF < 5).
#   5. Diagnostics: simulation-based (DHARMa-type) scaled residuals: KS
#      uniformity, dispersion, zero inflation, outliers, within-site lag-1
#      residual autocorrelation; VIF.
#      R2: Nakagawa marginal/conditional (Tweedie lognormal approximation),
#      Nagelkerke vs intercept-only (+RE), r2 observed vs fitted.
#   6. Robustness: log1p(y) LMM/LM with same predictors; leave-one-site-out
#      (multi-site) or leave-one-year-out (single site); AR(1) residual structure.
#
# TIMING_OVERRIDE / WINDOW_OVERRIDE below bypass the automatic selection.
#
# Columns expected in Sheet1: YL_<Pred>_<GP|SS> (previous year) and
# AN_<Pred>_<GP|SS> (current year), e.g. YL_CWB_GP, AN_CWB_GP.
#
# Tweedie power is estimated. If nlminb fails (Hessian not positive definite,
# typically power at the upper bound 2), the model is refitted with BFGS.
# Predictors are z-scored within the dataset: coefficients are the change in
# log(mortality) per 1 SD anomaly; exp(b) is the multiplicative change.
# =============================================================================

rm(list = ls())   # these scripts reuse object names; start from a clean workspace

# ---- settings ------------------------------------------------------------------
DATASET   <- "Switzerland"      # "Pooled", "Spain" or "Switzerland"
RESPONSE  <- "SimMor"     # "ObsMor" or "SimMor"
TIE_RATIO <- 2                  # evidence ratio below which timings/windows are tied
TIMINGS   <- c(YL = "previous year")#, AN = "current year")
TIMING_OVERRIDE <- "YL"   # NULL = automatic; or "YL", "AN", c("YL", "AN")
WINDOW_OVERRIDE <- "SS"   # NULL = automatic; or "GP", "SS", c("GP", "SS")
KMAX      <- 2           # max predictors in all-subsets (~ n / 10-15)
NSIM      <- 1000               # simulations for residual diagnostics
in_file   <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/ObsMor_anomaly.xlsx"
suppressMessages({
  library(readxl); library(glmmTMB); library(performance); library(dplyr)
  library(tidyr); library(lme4); library(ggplot2); library(patchwork); library(openxlsx)
})
set.seed(42)

TAG     <- paste(DATASET, RESPONSE, sep = "_")
out_dir <- file.path("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/results_new1", TAG); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

preds <- c("VPD", "Eto", "CWB", "WDcount", "Prec", "CDD", "SDII")
wins  <- c("SS") #"GP", 

d <- read_excel(in_file, sheet = "Sheet1") %>% filter(!is.na(Site)) %>% arrange(Site, Year)
if (DATASET != "Pooled") d <- filter(d, Country == DATASET)
d$Site <- factor(d$Site); d <- droplevels(d)
RE <- nlevels(d$Site) > 1                # site random effect only if >1 site
cat(TAG, "| n =", nrow(d), "| sites =", nlevels(d$Site), "| RE =", RE, "\n")

# ---- helpers -------------------------------------------------------------------
vname <- function(p, w, t) paste0(t, "_", p, "_", w)       # e.g. AN_CWB_GP
zname <- function(p, w, t) paste0("z_", vname(p, w, t))
lab   <- function(z) sub("^z_", "", z)                      # z_AN_CWB_GP -> AN_CWB_GP
need  <- as.vector(outer(outer(names(TIMINGS), preds, paste, sep = "_"), wins, paste, sep = "_"))
miss  <- setdiff(need, names(d))
if (length(miss)) stop("Missing columns in Sheet1: ", paste(miss, collapse = ", "))
for (t in names(TIMINGS)) for (p in preds) for (w in wins)
  d[[zname(p, w, t)]] <- as.numeric(scale(d[[vname(p, w, t)]]))

converged <- function(m) {
  if (is.null(m)) return(FALSE)
  ok <- isTRUE(m$fit$convergence == 0) && isTRUE(m$sdr$pdHess)
  ok && all(is.finite(sqrt(diag(vcov(m)$cond))))
}

fit_tw <- function(terms, df = d, re = RE, extra = NULL) {
  rhs <- if (length(terms)) paste(terms, collapse = " + ") else "1"
  if (re) rhs <- paste(rhs, "+ (1|Site)")
  if (!is.null(extra)) rhs <- paste(rhs, "+", extra)
  f <- as.formula(paste(RESPONSE, "~", rhs))
  m <- tryCatch(suppressWarnings(glmmTMB(f, family = tweedie(link = "log"), data = df)),
                error = function(e) NULL)
  meth <- "nlminb"
  if (!converged(m)) {
    m2 <- tryCatch(suppressWarnings(glmmTMB(f, family = tweedie(link = "log"), data = df,
                                            control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))),
                   error = function(e) NULL)
    if (converged(m2)) { m <- m2; meth <- "BFGS" }
  }
  if (!is.null(m)) attr(m, "fit_method") <- meth
  m
}

aicc <- function(m) {
  ll <- logLik(m); k <- attr(ll, "df"); n <- nobs(m)
  -2 * as.numeric(ll) + 2 * k + 2 * k * (k + 1) / max(n - k - 1, 1)
}

# Nakagawa R2 for Tweedie (log link), lognormal approximation:
#   sigma2_d = log(1 + phi * lambda^(p - 2)), lambda = exp(b0 + 0.5 * sigma2_RE)
nak_tw <- function(m) {
  b <- fixef(m)$cond; X <- getME(m, "X")
  s2f <- var(as.vector(X %*% b))
  s2r <- if (length(VarCorr(m)$cond)) sum(sapply(VarCorr(m)$cond, function(v) v[1])) else 0
  lambda <- unname(exp(b[1] + 0.5 * s2r))
  s2d <- log(1 + sigma(m) * lambda^(unname(family_params(m)) - 2))
  c(R2m = s2f / (s2f + s2r + s2d), R2c = (s2f + s2r) / (s2f + s2r + s2d))
}

r2_all <- function(m, m0, df = d) {
  y <- df[[RESPONSE]]; n <- length(y); nk <- nak_tw(m)
  ll <- as.numeric(logLik(m)); ll0 <- as.numeric(logLik(m0))
  cs <- 1 - exp(-2 / n * (ll - ll0))
  fc <- predict(m, type = "response"); fm <- predict(m, type = "response", re.form = NA)
  data.frame(R2_Nakagawa_marginal = nk[["R2m"]], R2_Nakagawa_conditional = nk[["R2c"]],
             R2_Nagelkerke_vs_null = cs / (1 - exp(2 / n * ll0)),
             R2_cor_fixed = cor(y, fm)^2, R2_cor_full = cor(y, fc)^2,
             RMSE = sqrt(mean((y - fc)^2)), MAE = mean(abs(y - fc)))
}

vif_max <- function(m) {
  if (length(fixef(m)$cond) <= 2) return(NA_real_)
  cc <- tryCatch(suppressWarnings(check_collinearity(m)), error = function(e) NULL)
  if (is.null(cc)) NA_real_ else max(cc$VIF)
}

sim_diag <- function(m, df = d, nsim = NSIM) {
  y <- df[[RESPONSE]]; n <- length(y)
  sims <- as.matrix(simulate(m, nsim = nsim, seed = 1))
  fit  <- predict(m, type = "response")
  res <- sapply(seq_len(n), function(i) {           # PIT residuals, ties randomised
    mean(sims[i, ] < y[i]) + runif(1) * mean(sims[i, ] == y[i])
  })
  ks <- suppressWarnings(ks.test(res, "punif"))
  obs_sd <- sd(y - fit); sim_sd <- apply(sims, 2, function(s) sd(s - fit))
  disp_p <- min(1, 2 * min(mean(sim_sd >= obs_sd), mean(sim_sd <= obs_sd)))
  z_obs <- sum(y == 0); z_sim <- colSums(sims == 0)
  zi_p <- if (z_obs > 0 || mean(z_sim) > 0) min(1, 2 * min(mean(z_sim >= z_obs), mean(z_sim <= z_obs))) else NA
  n_out <- sum(res < 1 / (nsim + 1) | res > nsim / (nsim + 1))
  zr <- qnorm(pmin(pmax(res, 1e-4), 1 - 1e-4))
  lag1 <- function(zz) {
    num <- 0; den <- 0
    for (s in unique(df$Site)) {
      idx <- which(df$Site == s); idx <- idx[order(df$Year[idx])]; v <- zz[idx] - mean(zz[idx])
      if (length(v) > 2) { num <- num + sum(v[-1] * v[-length(v)]); den <- den + sum(v^2) }
    }
    num / den
  }
  ac <- lag1(zr)
  ac_perm <- replicate(999, { zp <- zr
  for (s in unique(df$Site)) { idx <- which(df$Site == s); zp[idx] <- sample(zr[idx]) }
  lag1(zp) })
  list(res = res, fit = fit,
       tab = data.frame(KS_D = unname(ks$statistic), KS_p = ks$p.value,
                        Dispersion_ratio = obs_sd / mean(sim_sd), Dispersion_p = disp_p,
                        Zeros_obs = z_obs, Zeros_sim_mean = mean(z_sim),
                        ZI_ratio = if (mean(z_sim) > 0) z_obs / mean(z_sim) else NA, ZI_p = zi_p,
                        Outliers = n_out, Outlier_p = binom.test(n_out, n, p = 2 / (nsim + 1))$p.value,
                        Resid_lag1_r = ac, Resid_lag1_p = min(1, 2 * min(mean(ac_perm >= ac), mean(ac_perm <= ac)))))
}

# ---- Steps 0-1: timing (YL vs AN), then window (GP vs SS) ------------------------
akaike_w <- function(aic, a_null) {                 # weights over models + null
  amin <- min(c(aic, a_null)); w <- exp(-0.5 * (aic - amin))
  w / (sum(w) + exp(-0.5 * (a_null - amin)))
}
pick <- function(tab, col) {                        # best + tied options
  tab %>% arrange(desc(summed_Akaike_weight)) %>%
    mutate(Evidence_ratio_best_vs_this = summed_Akaike_weight[1] / summed_Akaike_weight,
           Selected = Evidence_ratio_best_vs_this < TIE_RATIO,
           Status = case_when(Selected & sum(Selected) > 1 ~ "selected (tied)",
                              Selected ~ "selected", TRUE ~ "not selected"))
}
evid <- function(x) summarise(x, n_sig_p05 = sum(p < 0.05), n_sig_p10 = sum(p < 0.10),
                              mean_dAICc_vs_null = mean(dAICc_vs_null),
                              best_predictor = paste(Predictor, Window)[which.min(AICc)],
                              best_dAICc_vs_null = min(dAICc_vs_null), best_p = p[which.min(AICc)],
                              summed_Akaike_weight = sum(wt), .groups = "drop")

cat("Steps 0-1: univariate models (timing x window x predictor)\n")
m0 <- fit_tw(character(0)); a0 <- aicc(m0)
grid <- expand.grid(Predictor = preds, Window = wins, Timing = names(TIMINGS), stringsAsFactors = FALSE)
step1 <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  p <- grid$Predictor[i]; w <- grid$Window[i]; t <- grid$Timing[i]
  tv <- zname(p, w, t); m <- fit_tw(tv); cf <- summary(m)$coefficients$cond
  data.frame(Timing = t, Window = w, Predictor = p, Variable = vname(p, w, t),
             Estimate = cf[tv, 1], SE = cf[tv, 2], z = cf[tv, 3], p = cf[tv, 4],
             AICc = aicc(m), AICc_null = a0, dAICc_vs_null = aicc(m) - a0,
             Tweedie_power = unname(family_params(m)), Fit = attr(m, "fit_method"),
             Converged = converged(m))
}))
step1$Akaike_weight_all <- akaike_w(step1$AICc, a0)

# Step 0: timing
timing_sum <- step1 %>% mutate(wt = Akaike_weight_all) %>% group_by(Timing) %>% evid() %>%
  mutate(Description = TIMINGS[Timing], .after = Timing) %>% pick()
sel_times <- timing_sum$Timing[timing_sum$Selected]
if (!is.null(TIMING_OVERRIDE)) {
  sel_times <- TIMING_OVERRIDE
  timing_sum$Selected <- timing_sum$Timing %in% sel_times
  timing_sum$Status <- ifelse(timing_sum$Selected, "selected (manual)", "not selected (manual)")
}
cat("Selected timing(s):", paste(sel_times, collapse = ", "), "\n"); print(as.data.frame(timing_sum))
#sel_times <- c("YL", "AN")   # uncomment to force previous-year (YL) models; use c("YL", "AN") for both
# Step 1: window, within each selected timing (weights renormalised within timing)
win_sum <- bind_rows(lapply(sel_times, function(t) {
  step1 %>% filter(Timing == t) %>% mutate(wt = akaike_w(AICc, a0)) %>%
    dplyr::group_by(Timing, Window) %>% evid() %>% pick()
}))
sel <- win_sum %>% dplyr::filter(Selected) %>% dplyr::select(Timing, Window)
if (!is.null(WINDOW_OVERRIDE)) {
  win_sum$Selected <- win_sum$Window %in% WINDOW_OVERRIDE
  win_sum$Status <- ifelse(win_sum$Selected, "selected (manual)", "not selected (manual)")
  sel <- win_sum %>% dplyr::filter(Selected) %>% dplyr::select(Timing, Window)
}
cat("Selected timing x window:", paste(sel$Timing, sel$Window, collapse = ", "), "\n")
print(as.data.frame(win_sum))

# ---- Steps 2-6 for each selected timing x window --------------------------------
screen_tab <- list(); cor_tab <- list(); subsets_tab <- list(); coef_tab <- list()
fit_tab <- list(); rob_tab <- list(); loso_tab <- list(); loyo_tab <- list(); ar1_tab <- list()
models <- list(); diags <- list()

for (j in seq_len(nrow(sel))) {
  t <- sel$Timing[j]; w <- sel$Window[j]
  key <- paste(TAG, t, w, sep = "_"); cat("Model:", key, "\n")
  u <- step1 %>% filter(Timing == t, Window == w) %>% arrange(AICc)
  
  # Step 2: screening
  rule <- "p < 0.05"; cand <- u$Predictor[u$p < 0.05 & u$Converged]
  if (!length(cand)) { rule <- "p < 0.10 (none at 0.05)"; cand <- u$Predictor[u$p < 0.10 & u$Converged] }
  if (!length(cand)) { rule <- "best AICc only (none at 0.10)"; cand <- u$Predictor[1] }
  
  # Step 3: collinearity filter
  cm <- cor(d[, sapply(preds, zname, w = w, t = t)]); dimnames(cm) <- list(preds, preds)
  kept <- character(0)
  for (p in cand) if (!length(kept) || all(abs(cm[p, kept]) <= 0.7)) kept <- c(kept, p)
  screen_tab[[key]] <- u %>% mutate(Model_ID = key, Screen_rule = rule,
                                    Passed_screen = Predictor %in% cand, Retained_after_r07 = Predictor %in% kept, .before = 1)
  cor_tab[[key]] <- data.frame(Model_ID = key, Var = preds, round(cm, 3), check.names = FALSE)
  
  # Step 4: all subsets
  combos <- list(character(0))
  for (k in seq_len(min(KMAX, length(kept)))) combos <- c(combos, combn(kept, k, simplify = FALSE))
  fits <- lapply(combos, function(cb) fit_tw(sapply(cb, zname, w = w, t = t)))
  st <- data.frame(Model_ID = key, Timing = t, Window = w,
                   Model = sapply(combos, function(cb) if (length(cb)) paste(cb, collapse = " + ") else "(null)"),
                   k = sapply(combos, length), AICc = sapply(fits, aicc), Converged = sapply(fits, converged),
                   Fit = sapply(fits, function(m) attr(m, "fit_method")))
  st$dAICc <- st$AICc - min(st$AICc[st$Converged])
  st$Akaike_weight <- ifelse(st$Converged, exp(-0.5 * st$dAICc), 0)
  st$Akaike_weight <- st$Akaike_weight / sum(st$Akaike_weight)
  ord <- order(st$AICc); st <- st[ord, ]; fits <- fits[ord]
  
  st$KS_p <- NA; st$Dispersion_p <- NA; st$VIF_max <- NA; st$Pass_diagnostics <- NA; chosen <- NA
  for (i in seq_len(nrow(st))) {
    if (!st$Converged[i] || st$k[i] == 0) next
    dg <- sim_diag(fits[[i]])
    st$KS_p[i] <- dg$tab$KS_p; st$Dispersion_p[i] <- dg$tab$Dispersion_p; st$VIF_max[i] <- vif_max(fits[[i]])
    st$Pass_diagnostics[i] <- dg$tab$KS_p > 0.05 && dg$tab$Dispersion_p > 0.05 &&
      (is.na(st$VIF_max[i]) || st$VIF_max[i] < 5)
    if (st$Pass_diagnostics[i]) { chosen <- i; break }
  }
  flag <- character(0)
  if (is.na(chosen)) { chosen <- which(st$Converged & st$k > 0)[1]
  flag <- c(flag, "no candidate passed all diagnostics; lowest-AICc shown") }
  if (st$Model[1] == "(null)") flag <- c(flag, "null model has lower AICc than best climate model")
  st$Final <- seq_len(nrow(st)) == chosen
  subsets_tab[[key]] <- st
  
  # Step 5: final model, coefficients, fit statistics, diagnostics
  mf <- fits[[chosen]]; models[[key]] <- mf
  dg <- sim_diag(mf); diags[[key]] <- dg
  cf <- summary(mf)$coefficients$cond
  ci <- suppressWarnings(confint(mf, parm = "beta_", method = "wald"))
  coef_tab[[key]] <- data.frame(Model_ID = key, Timing = t, Window = w, Term = lab(rownames(cf)),
                                Estimate = cf[, 1], SE = cf[, 2], CI_low = ci[rownames(cf), 1], CI_high = ci[rownames(cf), 2],
                                z = cf[, 3], p = cf[, 4], Ratio_per_SD = exp(cf[, 1]), row.names = NULL)
  rvi <- sapply(kept, function(p) sum(st$Akaike_weight[sapply(strsplit(st$Model, " \\+ "), function(x) p %in% x)]))
  fit_tab[[key]] <- data.frame(Model_ID = key, Dataset = DATASET, Response = RESPONSE,
                               Timing = t, Timing_status = timing_sum$Status[timing_sum$Timing == t], Window = w,
                               Window_status = win_sum$Status[win_sum$Timing == t & win_sum$Window == w],
                               Final_model = st$Model[chosen], n = nobs(mf), Random_effect = if (RE) "(1|Site)" else "none",
                               Site_RE_variance = if (RE) as.numeric(VarCorr(mf)$cond$Site) else NA,
                               Tweedie_power = unname(family_params(mf)), Dispersion_phi = sigma(mf),
                               Optimizer = attr(mf, "fit_method"),
                               AICc = st$AICc[chosen], AICc_null = a0, dAICc_vs_null = st$AICc[chosen] - a0,
                               LRT_vs_null_p = pchisq(2 * (as.numeric(logLik(mf)) - as.numeric(logLik(m0))),
                                                      df = attr(logLik(mf), "df") - attr(logLik(m0), "df"), lower.tail = FALSE),
                               n_models_dAICc_lt2 = sum(st$dAICc < 2 & st$Converged),
                               RVI = paste(sprintf("%s=%.2f", names(rvi), rvi), collapse = "; "),
                               r2_all(mf, m0), dg$tab, VIF_max = vif_max(mf),
                               Screen_rule = rule, Flag = paste(flag, collapse = "; "))
  
  terms <- names(fixef(mf)$cond)[-1]
  
  # Step 6a: log1p LMM / LM with the same predictors
  dd <- d; dd$ly <- log1p(dd[[RESPONSE]]); f <- paste("ly ~", paste(terms, collapse = " + "))
  if (RE) { ml <- lmer(as.formula(paste(f, "+ (1|Site)")), data = dd, REML = TRUE)
  cl <- summary(ml)$coefficients; pv <- 2 * pnorm(-abs(cl[, 3]))
  r2l <- r2_nakagawa(ml); r2m <- r2l$R2_marginal; r2c <- r2l$R2_conditional
  } else  { ml <- lm(as.formula(f), data = dd); cl <- summary(ml)$coefficients; pv <- cl[, 4]
  r2m <- r2c <- summary(ml)$r.squared }
  rob_tab[[key]] <- data.frame(Model_ID = key, Term = lab(rownames(cl)), Estimate_log1p = cl[, 1],
                               SE = cl[, 2], p = pv, R2_marginal = r2m, R2_conditional = r2c,
                               Shapiro_resid_p = shapiro.test(residuals(ml))$p.value,
                               Same_sign_as_Tweedie = ifelse(rownames(cl) == "(Intercept)", NA,
                                                             sign(cl[, 1]) == sign(fixef(mf)$cond[rownames(cl)])), row.names = NULL)
  
  # Step 6b: leave-one-site-out
  if (RE) loso_tab[[key]] <- bind_rows(lapply(levels(d$Site), function(s_out) {
    dl <- droplevels(filter(d, Site != s_out))
    ml <- fit_tw(terms, df = dl, re = nlevels(dl$Site) > 1); cl <- summary(ml)$coefficients$cond
    data.frame(Model_ID = key, Site_left_out = s_out, Term = lab(terms), Estimate = cl[terms, 1],
               SE = cl[terms, 2], p = cl[terms, 4], Converged = converged(ml), row.names = NULL)
  }))
  # Step 6b (single-site): leave-one-year-out
  if (!RE) loyo_tab[[key]] <- bind_rows(lapply(sort(unique(d$Year)), function(y_out) {
    dy <- filter(d, Year != y_out)
    ml <- fit_tw(terms, df = dy); cl <- summary(ml)$coefficients$cond
    data.frame(Model_ID = key, Year_left_out = y_out, Term = lab(terms), Estimate = cl[terms, 1],
               SE = cl[terms, 2], p = cl[terms, 4], Converged = converged(ml), row.names = NULL)
  }))
  
  # Step 6c: AR(1) residual structure within site
  d$fYear <- factor(d$Year)
  ma <- fit_tw(terms, extra = "ar1(fYear + 0 | Site)")
  ca <- summary(ma)$coefficients$cond
  rho <- tryCatch(attr(VarCorr(ma)$cond[[length(VarCorr(ma)$cond)]], "correlation")[1, 2], error = function(e) NA)
  ar1_tab[[key]] <- data.frame(Model_ID = key, Term = lab(terms),
                               Estimate_base = fixef(mf)$cond[terms], p_base = cf[terms, 4],
                               Estimate_AR1 = ca[terms, 1], SE_AR1 = ca[terms, 2], p_AR1 = ca[terms, 4], AR1_rho = rho,
                               AIC_base = AIC(mf), AIC_AR1 = AIC(ma),
                               LRT_AR1_p = pchisq(2 * (as.numeric(logLik(ma)) - as.numeric(logLik(mf))),
                                                  df = attr(logLik(ma), "df") - attr(logLik(mf), "df"), lower.tail = FALSE),
                               Converged_AR1 = converged(ma), row.names = NULL)
}

bind <- function(l) bind_rows(l)
screen_tab <- bind(screen_tab); cor_tab <- bind(cor_tab); subsets_tab <- bind(subsets_tab)
coef_tab <- bind(coef_tab); fit_tab <- bind(fit_tab); rob_tab <- bind(rob_tab)
loso_tab <- bind(loso_tab); loyo_tab <- bind(loyo_tab); ar1_tab <- bind(ar1_tab)

# ---- figures ---------------------------------------------------------------------
col <- c(ObsMor = "#2a78d6", SimMor = "#eb6834")[[RESPONSE]]
ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"
th <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
        axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink2),
        plot.title = element_text(colour = ink, face = "bold", size = 10),
        plot.subtitle = element_text(colour = ink2, size = 8),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))

for (key in names(models)) {
  mf <- models[[key]]; dg <- diags[[key]]; ft <- fit_tab[fit_tab$Model_ID == key, ]
  terms <- names(fixef(mf)$cond)[-1]; b <- fixef(mf)$cond; V <- vcov(mf)$cond; n <- nrow(d)
  
  # effect of each term (other terms at their mean = 0), population level
  eff <- lapply(terms, function(tv) {
    raw <- sub("^z_", "", tv); zz <- seq(min(d[[tv]]), max(d[[tv]]), length.out = 100)
    X <- matrix(0, 100, length(b), dimnames = list(NULL, names(b))); X[, 1] <- 1; X[, tv] <- zz
    eta <- as.vector(X %*% b); se <- sqrt(rowSums((X %*% V) * X))
    nd <- data.frame(x = zz * sd(d[[raw]]) + mean(d[[raw]]), fit = exp(eta),
                     lo = exp(eta - 1.96 * se), hi = exp(eta + 1.96 * se))
    ggplot() + geom_ribbon(data = nd, aes(x, ymin = lo, ymax = hi), fill = col, alpha = 0.18) +
      geom_point(data = d, aes(.data[[raw]], .data[[RESPONSE]]), colour = col, size = 1.6, alpha = 0.75) +
      geom_line(data = nd, aes(x, fit), colour = col, linewidth = 0.9) +
      labs(title = paste(key, "|", lab(tv)),
           subtitle = sprintf("b = %.3f, p = %s | R2m = %.2f, R2c = %.2f", b[tv],
                              format.pval(summary(mf)$coefficients$cond[tv, 4], digits = 2, eps = 1e-4),
                              ft$R2_Nakagawa_marginal, ft$R2_Nakagawa_conditional),
           x = paste(lab(tv), "anomaly"), y = paste0(RESPONSE, " (%)")) + th
  })
  ggsave(file.path(out_dir, paste0(key, "_effects.png")), wrap_plots(eff, ncol = min(2, length(eff))),
         width = 4.5 * min(2, length(eff)), height = 3.8 * ceiling(length(eff) / 2), dpi = 300)
  
  qq <- data.frame(e = (seq_len(n) - 0.5) / n, o = sort(dg$res))
  p1 <- ggplot(qq, aes(e, o)) + geom_abline(colour = muted, linetype = 2) + geom_point(colour = col, size = 1.4) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = key, subtitle = sprintf("KS p = %.2f | dispersion p = %.2f", ft$KS_p, ft$Dispersion_p),
         x = "Expected (uniform)", y = "Scaled residual") + th
  p2 <- ggplot(data.frame(pr = rank(dg$fit) / n, r = dg$res), aes(pr, r)) +
    geom_hline(yintercept = c(0.25, 0.5, 0.75), colour = muted, linetype = 2, linewidth = 0.3) +
    geom_point(colour = col, size = 1.4) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = ink2, linewidth = 0.6, span = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = " ", subtitle = sprintf("Residual lag-1 r = %.2f (p = %.3f)", ft$Resid_lag1_r, ft$Resid_lag1_p),
         x = "Predicted (rank)", y = "Scaled residual") + th
  lim <- range(c(d[[RESPONSE]], dg$fit))
  p3 <- ggplot(data.frame(f = dg$fit, y = d[[RESPONSE]]), aes(f, y)) + geom_abline(colour = muted, linetype = 2) +
    geom_point(colour = col, size = 1.4) + coord_equal(xlim = lim, ylim = lim) +
    labs(title = " ", subtitle = sprintf("r2 obs vs fitted (incl. RE) = %.2f", ft$R2_cor_full),
         x = "Fitted (%)", y = "Observed (%)") + th
  ggsave(file.path(out_dir, paste0(key, "_diagnostics.png")), p1 | p2 | p3, width = 10, height = 3.6, dpi = 300)
}

# ---- workbook --------------------------------------------------------------------
r4 <- function(df) { df <- as.data.frame(df); df[] <- lapply(df, function(v) if (is.numeric(v)) signif(v, 4) else v); df }
summary_tab <- fit_tab %>%
  left_join(coef_tab %>% filter(Term != "(Intercept)") %>% group_by(Model_ID) %>%
              summarise(Effect = paste(sprintf("%s: b = %.3f (95%% CI %.3f to %.3f), p = %s", Term, Estimate, CI_low, CI_high,
                                               format.pval(p, digits = 2, eps = 1e-4)), collapse = "; ")), by = "Model_ID") %>%
  left_join(ar1_tab %>% group_by(Model_ID) %>%
              summarise(AR1_rho = first(AR1_rho), AR1_LRT_p = first(LRT_AR1_p),
                        Effect_with_AR1 = paste(sprintf("%s: b = %.3f, p = %s", Term, Estimate_AR1,
                                                        format.pval(p_AR1, digits = 2, eps = 1e-4)), collapse = "; ")), by = "Model_ID") %>%
  transmute(Model_ID, Timing, Timing_status, Window, Window_status, Final_model, Effect, n, Random_effect,
            R2_marginal = R2_Nakagawa_marginal, R2_conditional = R2_Nakagawa_conditional,
            R2_Nagelkerke_vs_null, R2_cor_full, RMSE, dAICc_vs_null, LRT_vs_null_p,
            KS_p, Dispersion_p, ZI_p, Outliers, Resid_lag1_r, Resid_lag1_p,
            AR1_rho, AR1_LRT_p, Effect_with_AR1, Tweedie_power, Optimizer, Screen_rule, Flag)

wb <- createWorkbook()
add <- function(name, df) {
  if (is.null(df) || !nrow(df)) return(invisible())
  df <- r4(df); addWorksheet(wb, name); writeData(wb, name, df, withFilter = TRUE)
  freezePane(wb, name, firstRow = TRUE)
  addStyle(wb, name, createStyle(textDecoration = "bold", fgFill = "#e1e0d9"), rows = 1, cols = seq_len(ncol(df)))
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}
add("Summary", summary_tab); add("0_Timing_selection", timing_sum); add("1_Window_selection", win_sum)
add("1_Univariate", step1)
add("2_Screening", screen_tab); add("3_Correlations", cor_tab); add("4_All_subsets", subsets_tab)
add("5_Coefficients", coef_tab); add("5_Fit_diagnostics", fit_tab); add("6_Robust_log1p", rob_tab)
add("6_Robust_LOSO", loso_tab); add("6_Robust_LOYO", loyo_tab); add("6_AR1", ar1_tab)
saveWorkbook(wb, file.path(out_dir, paste0(TAG, "_results.xlsx")), overwrite = TRUE)

saveRDS(list(settings = list(DATASET = DATASET, RESPONSE = RESPONSE, TIE_RATIO = TIE_RATIO, KMAX = KMAX),
             data = d, step1 = step1, timing_sum = timing_sum, win_sum = win_sum, screen = screen_tab, cor = cor_tab,
             subsets = subsets_tab, coef = coef_tab, fit = fit_tab, rob = rob_tab, loso = loso_tab, loyo = loyo_tab,
             ar1 = ar1_tab, models = models, diag = diags),
        file.path(out_dir, paste0(TAG, "_results.rds")))
print(r4(summary_tab[, c("Model_ID", "Timing_status", "Window_status", "Final_model", "Effect", "R2_marginal", "R2_conditional")]))
cat("Done:", out_dir, "\n")

saveRDS(eff, file.path(out_dir, "fig_swiss_sim.rds"))