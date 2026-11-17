## ============================================================================
##  PROLES — Satellite-based redistribution of cumulative tree mortality
## ============================================================================
##
##  PURPOSE
##  -------
##  Converts single cumulative tree-mortality observations (Spanish plots,
##  2007–2015 census interval) into annual mortality time series, using a
##  mortality–vegetation-index (VI) model calibrated on the Swiss Visp site
##  (annual observations 1999–2020; 74 trees in 1999, 20 survivors in 2020).
##  The redistributed plot series are aggregated to sites and evaluated
##  against annual mortality simulated by the forest landscape model iLand.
##
##  PIPELINE (sections below mirror the accompanying report)
##  --------
##    1. Load data (VI_Swiss.xlsx, VI_Spain.xlsx, SiteMor_Spain.xlsx)
##    2. Feature engineering: per-series z-score anomalies of 8 VIs,
##       their 1-year lags, and year-over-year changes (24 features)
##    3. Candidate model zoo fitted on the Swiss site:
##       univariate / bivariate / trivariate OLS, Ridge, Lasso,
##       binomial GLM (logit hazard), random forest
##    4. Predictive power  = R², RMSE, Spearman rho of each model's
##       PREDICTED series vs. the OBSERVED ObsMor series (per project
##       specification).  Leave-one-out CV is computed as a supplementary
##       robustness indicator only.
##    5. Redistribution: model predictions on Spanish plot features
##       (2008–2015) are converted to non-negative weights; each plot's
##       cumulative deaths are allocated proportionally (mass-preserving).
##    6. Aggregation of plots to sites; evaluation against iLand by the
##       OVERALL pooled R² across all site-years (5 sites x 8 years),
##       with per-site Pearson r as supplementary detail.
##    7. Sensitivity: same evaluation with the first comparison year
##       (2008) removed.
##    8. Outputs: results tables (xlsx), redistributed series (xlsx),
##       and two figures in PROLES brand colours.
##
##  EXPECTED INPUT FILES (in DATA_DIR)
##  --------------------
##    VI_Swiss.xlsx      : columns Year, NDVI, NDWI, NBR, NDMI, EVI, NBR2,
##                         TCW, LST, ObsMor (annual mortality, % of live
##                         trees; NA in years without observation)
##    VI_Spain.xlsx      : columns Site, Plot, year + the same 8 VI columns,
##                         one row per plot-year
##    SiteMor_Spain.xlsx : sheet "ObsMor" (Site, Plot, InitialTrees,
##                         CumulativeDead), sheet "SimMor"
##                         (Simulation_site, Year, SimMor)
##
##  USAGE
##  -----
##    Rscript mortality_redistribution.R [DATA_DIR] [OUT_DIR]
##    (both default to the current working directory)
##
##  DEPENDENCIES
##  ------------
##    readxl, glmnet, randomForest, openxlsx   (all on CRAN)
##    install.packages(c("readxl","glmnet","randomForest","openxlsx"))
##
##  NOTES ON REPRODUCIBILITY
##  ------------------------
##  * Ridge is implemented in closed form with an unpenalized intercept,
##    which reproduces scikit-learn's Ridge exactly; Lasso uses glmnet,
##    whose Gaussian objective  (1/2n)*RSS + lambda*||b||_1  equals
##    scikit-learn's Lasso objective, so coefficients match.
##  * The random forest is stochastic and implementation-dependent
##    (randomForest vs. scikit-learn); its scores are approximate.
##  * All other models are deterministic and match the Python reference
##    pipeline to numerical precision.
## ============================================================================

suppressPackageStartupMessages({
  library(readxl)        # read .xlsx input
  library(glmnet)        # Lasso / Ridge paths
  library(randomForest)  # nonparametric benchmark
  library(openxlsx)      # write .xlsx output
})

## ---------------------------------------------------------------------------
## 0.  Configuration
## ---------------------------------------------------------------------------
args     <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- ifelse(length(args) >= 1, args[1], "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/RS/RS_DataForScript")  #"C:\Users\kumar_das\OneDrive - CZU v Praze\R_etal\PineDrought\RS\RS_DataForScript"
OUT_DIR  <- ifelse(length(args) >= 2, args[2], "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/RS/RS_DataForScript/Outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

VIS        <- c("NDVI","NDWI","NBR","NDMI","EVI","NBR2","TCW","LST")
YEARS_ES   <- 2008:2015     # deaths are assigned after the 2007 census,
                            # matching the iLand simulation window
SEED       <- 42            # random forest only
RUN_LOOCV  <- TRUE          # supplementary out-of-sample indicator

## PROLES brand colours 
COL_SMOLA   <- "#2B4334"    # dark green
COL_MECH    <- "darkblue"    # light green
COL_JEHLICI <- "#E19230"    # orange
COL_BRIZA   <- "#FAEED9"    # cream

## ---------------------------------------------------------------------------
## 1.  Load data
## ---------------------------------------------------------------------------
message("Loading data from: ", normalizePath(DATA_DIR))
swiss <- read_excel(file.path(DATA_DIR, "VI_Swiss.xlsx"))
spain <- read_excel(file.path(DATA_DIR, "VI_Spain.xlsx"))
obs   <- read_excel(file.path(DATA_DIR, "SiteMor_Spain.xlsx"), sheet = "ObsMor")
sim   <- read_excel(file.path(DATA_DIR, "SiteMor_Spain.xlsx"), sheet = "SimMor")
sim   <- sim[grepl("^SP", sim$Simulation_site), ]   # Spanish sites only

## ---------------------------------------------------------------------------
## 2.  Feature engineering
## ---------------------------------------------------------------------------
## Absolute VI values do not transfer between the Swiss conifer stand and the
## drier Spanish stands (different baseline reflectance), so every predictor
## is expressed relative to its own series:
##   <VI>_z   : z-score anomaly within the series (population SD)
##   <VI>_zl1 : the anomaly lagged by one year
##   <VI>_dz  : year-over-year change of the anomaly
## -> 8 VIs x 3 variants = 24 candidate features.

#' Add anomaly / lag / delta features to one time-ordered series.
#' @param df data.frame with a year column and the 8 VI columns
#' @return   df with 24 additional feature columns
add_features <- function(df) {
  for (v in VIS) {
    z <- (df[[v]] - mean(df[[v]])) / sd_pop(df[[v]])
    df[[paste0(v, "_z")]]   <- z
    df[[paste0(v, "_zl1")]] <- c(NA, head(z, -1))
    df[[paste0(v, "_dz")]]  <- z - c(NA, head(z, -1))
  }
  df
}

#' Population standard deviation (ddof = 0), matching the reference pipeline.
sd_pop <- function(x) sqrt(mean((x - mean(x))^2))

FEATS <- as.vector(outer(VIS, c("_z", "_zl1", "_dz"), paste0))

## --- Swiss training table ---------------------------------------------------
swiss <- add_features(swiss[order(swiss$Year), ])
swiss$rate <- swiss$ObsMor / 100                       # response as fraction
D <- swiss[!is.na(swiss$rate) & complete.cases(swiss[, FEATS]), ]
row.names(D) <- NULL
message(sprintf("Swiss training rows: %d (%d-%d)", nrow(D), min(D$Year), max(D$Year)))

## Approximate live-tree pool per year (compounding the annual hazard from the
## initial 74 trees).  Used only as prior weights in the binomial GLM.
D$alive <- 74 * cumprod(c(1, head(1 - D$rate, -1)))

y     <- D$rate
X_all <- as.matrix(D[, FEATS])

## --- Spanish plot features --------------------------------------------------
## Anomalies are computed per plot over its full 2000-2015 series, then the
## modelled window 2008-2015 is extracted.
spain_feat <- do.call(rbind, lapply(split(spain, spain$Plot), function(g) {
  add_features(g[order(g$year), ])
}))

## ---------------------------------------------------------------------------
## 3.  Candidate models
## ---------------------------------------------------------------------------
## Each model is described by a list(name, family, feats).  fit_predict()
## fits the model on the full Swiss record and returns predictions for
## arbitrary new rows; predictions used as mortality rates are later clipped
## to [0, 1] (Swiss evaluation) or truncated at 0 (redistribution weights).

#' Exact scikit-learn-equivalent ridge regression (intercept unpenalized).
#' Solves  min ||y - b0 - Xb||^2 + a * ||b||^2  in closed form.
#' @return list(b0 = intercept, b = coefficient vector)
ridge_fit <- function(X, y, a) {
  Xc <- scale(X, center = TRUE, scale = FALSE)
  yc <- y - mean(y)
  b  <- solve(crossprod(Xc) + diag(a, ncol(X)), crossprod(Xc, yc))
  list(b0 = mean(y) - as.numeric(colMeans(X) %*% b), b = as.numeric(b))
}

#' Fit one candidate on training data (Xtr, ytr) and predict for Xnew.
#' @param mod  list(name, family, feats); family one of
#'             "OLS", "Ridge", "Lasso", "RF", "GLM"
#' @param w    prior weights for the binomial GLM (live-tree counts)
fit_predict <- function(mod, Xtr, ytr, Xnew, w = NULL) {
  Xtr  <- Xtr[,  mod$feats, drop = FALSE]
  Xnew <- Xnew[, mod$feats, drop = FALSE]
  switch(mod$family,
    OLS = {
      fit <- lm.fit(cbind(1, Xtr), ytr)
      as.numeric(cbind(1, Xnew) %*% fit$coefficients)
    },
    Ridge = {
      f <- ridge_fit(Xtr, ytr, mod$a)
      as.numeric(f$b0 + Xnew %*% f$b)
    },
    Lasso = {
      ## glmnet Gaussian objective (1/2n)RSS + lambda*||b||_1 == sklearn Lasso
      fit <- glmnet(Xtr, ytr, alpha = 1, lambda = mod$a,
                    standardize = FALSE, intercept = TRUE)
      as.numeric(predict(fit, newx = Xnew))
    },
    RF = {
      set.seed(SEED)
      fit <- randomForest(x = Xtr, y = ytr, ntree = 400,
                          nodesize = 3, maxnodes = 8)   # ~ depth-3 trees
      as.numeric(predict(fit, Xnew))
    },
    GLM = {
      df  <- data.frame(rate = ytr, Xtr, check.names = FALSE)
      fit <- suppressWarnings(
        glm(reformulate(sprintf("`%s`", mod$feats), "rate"),
            data = df, family = binomial(), weights = w))
      as.numeric(predict(fit, newdata = as.data.frame(Xnew), type = "response"))
    })
}

## --- Evaluation helpers ------------------------------------------------------

#' Skill of a predicted vs. an observed series.
#' R2 is the coefficient of determination 1 - SSE/SST (not squared r).
series_skill <- function(obs, pred) {
  pred <- pmin(pmax(pred, 0), 1)
  c(R2   = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
    RMSE = sqrt(mean((obs - pred)^2)) * 100,              # percentage points
    rho  = suppressWarnings(cor(obs, pred, method = "spearman")))
}

#' Supplementary leave-one-out CV R2 (out-of-sample robustness indicator).
loocv_r2 <- function(mod) {
  pred <- vapply(seq_along(y), function(i) {
    fit_predict(mod, X_all[-i, , drop = FALSE], y[-i],
                X_all[i, , drop = FALSE], w = D$alive[-i])
  }, numeric(1))
  pred <- pmin(pmax(pred, 0), 1)
  1 - sum((y - pred)^2) / sum((y - mean(y))^2)
}

## --- Build the model zoo -----------------------------------------------------
message("Screening candidate models ...")
zoo <- list()

## (i) univariate OLS on all 24 features
for (f in FEATS)
  zoo[[length(zoo) + 1]] <- list(name = paste0("OLS[", f, "]"),
                                 family = "OLS", feats = f)

uni_r2 <- vapply(zoo, function(m)
  series_skill(y, fit_predict(m, X_all, y, X_all))["R2"], numeric(1))
top10  <- FEATS[order(-uni_r2)][1:10]

## (ii) bivariate OLS among the top-10 single features
for (cmb in combn(top10, 2, simplify = FALSE))
  zoo[[length(zoo) + 1]] <- list(name = paste0("OLS[", paste(cmb, collapse = "+"), "]"),
                                 family = "OLS", feats = cmb)

## (iii) trivariate OLS among the top-6
for (cmb in combn(top10[1:6], 3, simplify = FALSE))
  zoo[[length(zoo) + 1]] <- list(name = paste0("OLS[", paste(cmb, collapse = "+"), "]"),
                                 family = "OLS", feats = cmb)

## (iv) Ridge and (v) Lasso on all 24 features
for (a in c(0.5, 1, 3, 10))
  zoo[[length(zoo) + 1]] <- list(name = sprintf("Ridge(a=%g)", a),
                                 family = "Ridge", feats = FEATS, a = a)
for (a in c(0.002, 0.005, 0.01, 0.02))
  zoo[[length(zoo) + 1]] <- list(name = sprintf("Lasso(a=%g)", a),
                                 family = "Lasso", feats = FEATS, a = a)

## (vi) random forest benchmark
zoo[[length(zoo) + 1]] <- list(name = "RandomForest", family = "RF", feats = FEATS)

## (vii) binomial GLM (logit hazard) for the 5 strongest single features
##       and the 3 strongest OLS pairs
for (f in top10[1:5])
  zoo[[length(zoo) + 1]] <- list(name = paste0("GLM-binom[", f, "]"),
                                 family = "GLM", feats = f)
is_pair  <- vapply(zoo, function(m) m$family == "OLS" && length(m$feats) == 2, logical(1))
pair_r2  <- vapply(zoo[is_pair], function(m)
  series_skill(y, fit_predict(m, X_all, y, X_all))["R2"], numeric(1))
for (m in zoo[is_pair][order(-pair_r2)][1:3])
  zoo[[length(zoo) + 1]] <- list(name = paste0("GLM-binom[", paste(m$feats, collapse = "+"), "]"),
                                 family = "GLM", feats = m$feats)

## ---------------------------------------------------------------------------
## 4.  Predictive power on the Swiss site (predicted vs. observed series)
## ---------------------------------------------------------------------------
skill <- t(vapply(zoo, function(m) {
  series_skill(y, fit_predict(m, X_all, y, X_all, w = D$alive))
}, numeric(3)))
swiss_res <- data.frame(model  = vapply(zoo, `[[`, "", "name"),
                        family = vapply(zoo, `[[`, "", "family"),
                        fit_R2 = skill[, "R2"], fit_RMSE = skill[, "RMSE"],
                        fit_rho = skill[, "rho"])

if (RUN_LOOCV) {
  message("Computing supplementary LOOCV (this is the slow part) ...")
  swiss_res$loocv_R2 <- vapply(zoo, loocv_r2, numeric(1))
}
swiss_res <- swiss_res[order(-swiss_res$fit_R2), ]

## Carry the 14 best-fitting candidates forward to Spain.
carry_names <- head(swiss_res$model, 14)
carry <- zoo[match(carry_names, vapply(zoo, `[[`, "", "name"))]

## ---------------------------------------------------------------------------
## 5.  Redistribution of Spanish cumulative mortality
## ---------------------------------------------------------------------------
## For each plot: predict a mortality intensity p(t) for 2008-2015, truncate
## negatives at zero, normalise to weights w(t) = p(t)/sum(p), and allocate
##   deaths(t) = CumulativeDead * w(t).
## The allocation is mass-preserving: the satellite signal controls only the
## timing of mortality, never its amount.  If all predictions are <= 0 the
## weights default to uniform.

#' Redistribute all plots under one model.
#' @return data.frame Plot, Site, Year, Deaths, InitialTrees
redistribute <- function(mod) {
  out <- lapply(split(spain_feat, spain_feat$Plot), function(g) {
    g   <- g[g$year %in% YEARS_ES, ]
    g   <- g[order(g$year), ]
    p   <- fit_predict(mod, X_all, y, as.matrix(g[, FEATS]), w = D$alive)
    w   <- pmax(p, 0)
    w   <- if (sum(w) <= 0) rep(1 / length(w), length(w)) else w / sum(w)
    meta <- obs[obs$Plot == g$Plot[1], ]
    data.frame(Plot = g$Plot[1], Site = meta$Site, Year = g$year,
               Deaths = meta$CumulativeDead * w,
               InitialTrees = meta$InitialTrees)
  })
  do.call(rbind, out)
}

#' Aggregate plot allocations to site-level annual mortality (% of the
#' site's initial tree pool).
aggregate_site <- function(pl) {
  agg <- aggregate(cbind(Deaths, InitialTrees) ~ Site + Year, pl, sum)
  agg$RedistMor <- 100 * agg$Deaths / agg$InitialTrees
  agg[order(agg$Site, agg$Year), ]
}

## ---------------------------------------------------------------------------
## 6.  Evaluation against iLand
## ---------------------------------------------------------------------------
## Headline metric: overall R2 (squared Pearson correlation) between the
## redistributed and simulated series across ALL pooled site-years
## (5 sites x 8 years = 40 points).  Per-site Pearson r is supplementary.

#' Pooled and per-site agreement of one site table with iLand.
eval_iland <- function(st, year_min = min(YEARS_ES)) {
  m <- merge(st, sim, by.x = c("Site", "Year"),
             by.y = c("Simulation_site", "Year"))
  m <- m[m$Year >= year_min, ]
  r_pool <- cor(m$RedistMor, m$SimMor)
  list(overall = c(R2 = r_pool^2, r = r_pool,
                   RMSE = sqrt(mean((m$RedistMor - m$SimMor)^2))),
       per_site = do.call(rbind, lapply(split(m, m$Site), function(g)
         data.frame(Site = g$Site[1],
                    r    = cor(g$RedistMor, g$SimMor),
                    rho  = suppressWarnings(cor(g$RedistMor, g$SimMor,
                                                method = "spearman")),
                    RMSE = sqrt(mean((g$RedistMor - g$SimMor)^2))))))
}

message("Redistributing and evaluating ", length(carry), " candidate models ...")
site_tables <- list(); eval_tab <- NULL
for (mod in carry) {
  st <- aggregate_site(redistribute(mod))
  site_tables[[mod$name]] <- st
  ev <- eval_iland(st)
  eval_tab <- rbind(eval_tab, data.frame(model = mod$name, t(ev$overall)))
}
results <- merge(swiss_res, eval_tab, by = "model", all.x = TRUE)
results <- results[order(-results$fit_R2), ]
names(results)[names(results) %in% c("R2","r","RMSE")] <-
  c("iLand_R2", "iLand_r", "iLand_RMSE")

## Recommended model: best overall iLand R2 among the well-performing set.
best_name <- results$model[which.max(results$iLand_R2)]
message("\nRecommended model (best overall iLand R2): ", best_name)
best_st   <- site_tables[[best_name]]
best_eval <- eval_iland(best_st)

## ---------------------------------------------------------------------------
## 7.  Sensitivity: first comparison year (2008) removed
## ---------------------------------------------------------------------------
sens_eval <- eval_iland(best_st, year_min = 2009)

## ---------------------------------------------------------------------------
## 8.  Outputs
## ---------------------------------------------------------------------------
## 8a. Tables -----------------------------------------------------------------
best_pl <- redistribute(carry[[match(best_name,
                                     vapply(carry, `[[`, "", "name"))]])
best_pl$AnnualMortality_pctInitial <- 100 * best_pl$Deaths / best_pl$InitialTrees
best_st_out <- merge(best_st,
                     setNames(sim, c("Site", "Year", "iLand_SimMor")),
                     by = c("Site", "Year"), all.x = TRUE)

wb <- createWorkbook()
hs <- createStyle(fontName = "Mulish", textDecoration = "bold",
                  fontColour = COL_BRIZA, fgFill = COL_SMOLA,
                  halign = "center")
add_sheet <- function(name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hs)
  setColWidths(wb, name, cols = seq_along(df), widths = "auto")
  freezePane(wb, name, firstRow = TRUE)
}
add_sheet("Model_Comparison",  results)
add_sheet("Site_Annual_Mortality", best_st_out)
add_sheet("Plot_Annual_Mortality", best_pl)
add_sheet("BestModel_PerSite_vs_iLand", best_eval$per_site)
add_sheet("Sensitivity_no2008", sens_eval$per_site)
saveWorkbook(wb, file.path(OUT_DIR, "R_Redistributed_Mortality1.xlsx"),
             overwrite = TRUE)

## 8b. Console summary --------------------------------------------------------
fmt <- function(x, d = 3) formatC(x, digits = d, format = "f")
cat("\n================ Swiss predictive power (top 14) ================\n")
print(within(head(results, 14), {
  fit_R2 <- fmt(fit_R2); fit_RMSE <- fmt(fit_RMSE, 2); fit_rho <- fmt(fit_rho, 2)
  iLand_R2 <- fmt(iLand_R2); iLand_r <- fmt(iLand_r, 2); iLand_RMSE <- fmt(iLand_RMSE, 2)
})[, c("model","fit_R2","fit_RMSE","fit_rho","iLand_R2","iLand_r","iLand_RMSE")],
row.names = FALSE)

cat("\n============ Recommended model vs. iLand, per site ==============\n")
print(within(best_eval$per_site,
             { r <- fmt(r, 2); rho <- fmt(rho, 2); RMSE <- fmt(RMSE, 2) }),
      row.names = FALSE)
cat(sprintf("\nPooled 2008-2015 : R2 = %s, r = %s, RMSE = %s %%/yr\n",
            fmt(best_eval$overall["R2"]), fmt(best_eval$overall["r"], 2),
            fmt(best_eval$overall["RMSE"], 2)))
cat(sprintf("Pooled 2009-2015 : R2 = %s, r = %s, RMSE = %s %%/yr  (first year removed)\n",
            fmt(sens_eval$overall["R2"]), fmt(sens_eval$overall["r"], 2),
            fmt(sens_eval$overall["RMSE"], 2)))

## 8c. Figures (base graphics, PROLES brand colours) --------------------------
## Figure A: redistributed vs. iLand per site, recommended model
png(file.path(OUT_DIR, "R_fig_sites_vs_iland1.png"),
    width = 2100, height = 1150, res = 200)
par(mfrow = c(2, 3), mar = c(3.5, 3.8, 2.5, 1), mgp = c(2.2, 0.7, 0),
    col.axis = COL_SMOLA, col.lab = COL_SMOLA, col.main = COL_SMOLA, bty = "l")
for (s in unique(best_st$Site)) {
  g  <- best_st[best_st$Site == s, ]
  sv <- sim[sim$Simulation_site == s, ]
  sv <- sv[order(sv$Year), ]
  ylim <- c(0, max(g$RedistMor, sv$SimMor) * 1.05)
  bp <- barplot(g$RedistMor, names.arg = g$Year, col = COL_MECH, border = NA,
                ylim = ylim, ylab = "Mortality (%/yr)",
                main = sprintf("%s  (r = %.2f)", s,
                               cor(g$RedistMor, sv$SimMor)))
  lines(bp, sv$SimMor, col = COL_JEHLICI, lwd = 2.5, type = "o", pch = 16)
}
plot.new()
legend("center", bty = "n", cex = 1.2, text.col = COL_SMOLA,
       legend = c("Redistributed observation", "iLand simulation"),
       fill = c(COL_MECH, NA), border = NA,
       lty = c(NA, 1), lwd = c(NA, 2.5), pch = c(NA, 16),
       col = c(NA, COL_JEHLICI))
invisible(dev.off())

## Figure B: predicted vs. observed mortality series, Swiss site
best_mod  <- carry[[match(best_name, vapply(carry, `[[`, "", "name"))]]
fit_swiss <- pmin(pmax(fit_predict(best_mod, X_all, y, X_all, w = D$alive), 0), 1)
png(file.path(OUT_DIR, "R_fig_swiss_fit1.png"),
    width = 1800, height = 750, res = 200)
par(mar = c(3.5, 3.8, 2.5, 1), mgp = c(2.2, 0.7, 0),
    col.axis = COL_SMOLA, col.lab = COL_SMOLA, col.main = COL_SMOLA, bty = "l")
bp <- barplot(y * 100, names.arg = D$Year, col = COL_MECH, border = NA,
              ylab = "Mortality (% of live trees)",
              main = sprintf("Swiss site: predicted vs. observed  (R2 = %.2f)",
                             series_skill(y, fit_swiss)["R2"]))
lines(bp, fit_swiss * 100, col = COL_JEHLICI, lwd = 2.5, type = "o", pch = 16)
legend("top", bty = "n", text.col = COL_SMOLA, horiz = TRUE,
       legend = c("Observed (ObsMor)", "Predicted"),
       fill = c(COL_MECH, NA), border = NA,
       lty = c(NA, 1), lwd = c(NA, 2.5), pch = c(NA, 16),
       col = c(NA, COL_JEHLICI))
invisible(dev.off())

message("\nOutputs written to ", normalizePath(OUT_DIR), ":")
message("  R_Redistributed_Mortality1.xlsx")
message("  R_fig_sites_vs_iland.png, R_fig_swiss_fit1.png")



#########################################################################################
#' Coefficient table of one Swiss-fitted candidate (NULL for random forest).
coef_table <- function(mod) {
  Xtr <- X_all[, mod$feats, drop = FALSE]
  b <- switch(mod$family,
              OLS = {
                f <- lm.fit(cbind(1, Xtr), y)
                setNames(f$coefficients, c("(Intercept)", mod$feats))
              },
              Ridge = {
                f <- ridge_fit(Xtr, y, mod$a)
                setNames(c(f$b0, f$b), c("(Intercept)", mod$feats))
              },
              Lasso = {
                f  <- glmnet(Xtr, y, alpha = 1, lambda = mod$a,
                             standardize = FALSE, intercept = TRUE)
                cf <- as.matrix(coef(f))
                setNames(cf[, 1], rownames(cf))
              },
              GLM = {
                df <- data.frame(rate = y, Xtr, check.names = FALSE)
                f  <- suppressWarnings(
                  glm(reformulate(sprintf("`%s`", mod$feats), "rate"),
                      data = df, family = binomial(), weights = D$alive))
                cf <- coef(f)
                names(cf) <- gsub("`", "", names(cf))
                cf
              },
              RF = NULL)
  if (is.null(b)) return(NULL)
  data.frame(model = mod$name, family = mod$family,
             term = names(b), estimate = unname(b),
             row.names = NULL)
}

#' Readable model expression; Lasso terms shrunk to exactly 0 are dropped.
equation <- function(ct) {
  b   <- setNames(ct$estimate, ct$term)
  b   <- b[names(b) == "(Intercept)" | b != 0]
  lhs <- if (ct$family[1] == "GLM") "logit(p)" else "p"
  rhs <- paste0(sprintf("%.5f", b[1]),
                paste0(sprintf(" %s %.5f*%s",
                               ifelse(b[-1] < 0, "-", "+"),
                               abs(b[-1]), names(b)[-1]),
                       collapse = ""))
  paste(lhs, "=", rhs)
}

## --- All 14 carried models ---------------------------------------------------
coef_all <- do.call(rbind, lapply(carry, coef_table))

cat("\n================ Model expressions (Swiss fit) ==================\n")
for (nm in unique(coef_all$model)) {
  cat("\n", nm, "\n  ", equation(coef_all[coef_all$model == nm, ]), "\n", sep = "")
}

## --- Random forest: no coefficients, report variable importance instead ------
rf_imp <- NULL
if (any(vapply(carry, `[[`, "", "family") == "RF")) {
  set.seed(SEED)
  rf_fit <- randomForest(x = X_all, y = y, ntree = 400,
                         nodesize = 3, maxnodes = 8)
  imp    <- importance(rf_fit)
  rf_imp <- data.frame(feature = rownames(imp),
                       IncNodePurity = imp[, 1], row.names = NULL)
  rf_imp <- rf_imp[order(-rf_imp$IncNodePurity), ]
  cat("\nRandomForest has no coefficients; top 10 features by importance:\n")
  print(head(rf_imp, 10), row.names = FALSE)
}

## --- Recommended model: full regression output -------------------------------
cat("\n================ Recommended model: full summary ================\n")
cat(best_mod$name, "\n")
if (best_mod$family == "OLS") {
  best_lm <- lm(reformulate(best_mod$feats, "rate"), data = D)
  print(summary(best_lm))           # estimates, SEs, t, p, R2, adj. R2
  print(confint(best_lm))           # 95% CIs
} else {
  print(coef_table(best_mod))
}

## --- Add to the results workbook ---------------------------------------------
add_sheet("Model_Coefficients", coef_all)
if (!is.null(rf_imp)) add_sheet("RF_Importance", rf_imp)
saveWorkbook(wb, file.path(OUT_DIR, "R_Redistributed_Mortality1.xlsx"),
             overwrite = TRUE)
message("Coefficient sheets added to R_Redistributed_Mortality1.xlsx")
################################################################################################################
## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## 10.  Additional figures
## ---------------------------------------------------------------------------
## Paste at the END of mortality_redistribution.R and run in the same session.
## Uses objects already created: best_mod, best_name, X_all, y, D, best_st,
## fit_predict, series_skill, COL_*, OUT_DIR.  Base graphics only.

fit_swiss <- pmin(pmax(fit_predict(best_mod, X_all, y, X_all, w = D$alive), 0), 1)
sk        <- series_skill(y, fit_swiss)
obs_pct   <- y * 100
fit_pct   <- fit_swiss * 100

par_brand <- function(...)
  par(mgp = c(2.2, 0.7, 0), bty = "l", las = 1,
      col.axis = COL_SMOLA, col.lab = COL_SMOLA, col.main = COL_SMOLA,
      fg = COL_SMOLA, ...)

## ---- Figure C: Swiss time series + fitted-vs-observed scatter --------------
png(file.path(OUT_DIR, "R_fig_swiss_fit_scatter1.png"),
    width = 2600, height = 900, res = 200)
layout(matrix(1:2, nrow = 1), widths = c(2.3, 1))

## (a) time series
par_brand(mar = c(3.5, 4, 2.5, 1))
ymax <- max(obs_pct, fit_pct) * 1.3
bp <- barplot(obs_pct, names.arg = D$Year, col = COL_MECH, border = NA,
              ylim = c(0, ymax), ylab = "Mortality (% of live trees)",
              main = "(a) Swiss site: observed and fitted series")
lines(bp, fit_pct, col = COL_JEHLICI, lwd = 2.5, type = "o", pch = 16)
legend("topleft", bty = "n", text.col = COL_SMOLA, horiz = TRUE,
       legend = c("Observed", paste0("Fitted: ", best_name)),
       fill = c(COL_MECH, NA), border = NA,
       lty = c(NA, 1), lwd = c(NA, 2.5), pch = c(NA, 16),
       col = c(NA, COL_JEHLICI))

## (b) scatter, fitted vs observed, with 1:1 line
par_brand(mar = c(3.5, 4, 2.5, 1.5))
lim <- c(-0.03, 1.05) * max(obs_pct, fit_pct)
plot(obs_pct, fit_pct, xlim = lim, ylim = lim, xaxs = "i", yaxs = "i",
     pch = 21, bg = COL_JEHLICI, col = COL_SMOLA, cex = 1.3,
     xlab = "Observed mortality (%)", ylab = "Fitted mortality (%)",
     main = "(b) Fitted vs. observed")
abline(0, 1, lty = 2, col = COL_SMOLA)
legend("topleft", bty = "n", text.col = COL_SMOLA, inset = 0.02,
       legend = c(as.expression(bquote(R^2 == .(sprintf("%.2f", sk["R2"])))),
                  as.expression(bquote(rho  == .(sprintf("%.2f", sk["rho"])))),
                  "1:1 line"),
       lty = c(NA, NA, 2), col = c(NA, NA, COL_SMOLA))
invisible(dev.off())

## ---- Figure D: Swiss observed + Spanish redistributed, one time axis -------
sites    <- unique(best_st$Site)
site_col <- c(COL_JEHLICI, "#8C5A2B", "#6F9A4E", "#C9A227", "#4F7A8C",
              "#A0522D", "#7B6D8D")[seq_along(sites)]
site_pch <- c(16, 17, 15, 18, 25, 8, 4)[seq_along(sites)]

yrs  <- range(c(D$Year, best_st$Year))
ymax <- max(obs_pct, best_st$RedistMor) * 1.08

png(file.path(OUT_DIR, "R_fig_timeseries_swiss_spain1.png"),
    width = 1700, height = 1000, res = 200)
par_brand(mar = c(3.5, 4, 0.8, 0.8))
plot(NA, xlim = yrs, ylim = c(0, ymax), xaxt = "n",
     xlab = "Year", ylab = "Annual mortality (%)")
axis(1, at = seq(yrs[1], yrs[2], 1),
     labels = ifelse(seq(yrs[1], yrs[2], 1) %% 2 == 0,
                     seq(yrs[1], yrs[2], 1), NA))
## shade the Spanish redistribution window
rect(min(best_st$Year) - 0.5, 0, max(best_st$Year) + 0.5, ymax,
     col = COL_BRIZA, border = NA)
box(bty = "l")

lines(D$Year, obs_pct, col = COL_SMOLA, lwd = 3, type = "o", pch = 21,
      bg = COL_SMOLA)
for (i in seq_along(sites)) {
  g <- best_st[best_st$Site == sites[i], ]
  g <- g[order(g$Year), ]
  lines(g$Year, g$RedistMor, col = site_col[i], lwd = 2, type = "o",
        pch = site_pch[i], bg = site_col[i])
}

## legend inside the plot: two columns, transparent, top-left (in the
## empty space between the 1999 peak and the 2017-2019 peaks).
## Adjust LEG_X (year) if it collides with your data.
LEG_X <- min(D$Year) + 0.7
legend(x = LEG_X, y = ymax, xjust = 0, yjust = 1,
       ncol = 2, bty = "n", text.col = COL_SMOLA, cex = 0.75,
       seg.len = 1.5, x.intersp = 0.5, y.intersp = 1.1,
       legend = c("CH_Visp, observed",  paste0(sites)),
       col = c(COL_SMOLA, site_col, "#D9C49E"),
       lwd = c(3, rep(2, length(sites)), NA),
       pch = c(21, site_pch, 22),
       pt.bg = c(COL_SMOLA, site_col, COL_BRIZA),
       pt.cex = c(1, rep(1, length(sites)), 2))
invisible(dev.off())

message("Figures written: R_fig_swiss_fit_scatter1.png, ",
        "R_fig_timeseries_swiss_spain1.png")
