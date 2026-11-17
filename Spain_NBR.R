# =============================================================================
# Spain: annual mortality allocation (2007-2015) using NBR + NDVI AGREEMENT
# -----------------------------------------------------------------------------
# Same hazard model and the SAME eight signal rules as the Swiss script
# (swiss_consensus_allocation.R), where they were validated against ObsMor:
#
#   NBR, NDVI, AND_NBR, AND_meanz, AND_geomz, ANDlag_NBR, OR_meanz, OR_maxz
#
# Per plot: N0 = TotalTrees (1 Jan 2007), NT = TotalTrees - Dead (end 2015)
#   m_t = 1 - exp(-beta * s_t),  beta = ln(N0/NT) / sum(s_t)  -> plot totals exact
# Scheme A: plot's own NBR/NDVI series.  Scheme B: site-mean NBR/NDVI series,
# applied to each plot. Site series = sum of plots -> site totals exact.
# If a plot has deaths but zero signal under a rule (Scheme A), it falls back to
# the site signal for that rule (flagged).
#
# Annual mortality % = 100 * deaths_t / trees alive at start of year t (as ObsMor)
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
})

# ---- Settings ---------------------------------------------------------------
f_index    <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/RS/NBR_Spain.xlsx"
f_mort     <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/RS/SiteMor_Spain.xlsx"
out_dir    <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/RS/spain_consensus_output/new1"
year_start <- 2007
year_end   <- 2015
window     <- 2
dir.create(out_dir, showWarnings = FALSE)

# ---- Functions (identical to the Swiss script) --------------------------------
drop2 <- function(x, w = 2) {
  s <- rep(NA_real_, length(x))
  for (t in seq_along(x)) if (t > w) s[t] <- max(0, max(x[(t - w):(t - 1)]) - x[t])
  s
}

consensus_signals <- function(nbr, ndmi, w = 2) {
  dn <- drop2(nbr, w);  dv <- drop2(ndmi, w)
  zn <- dn / sd(diff(nbr));  zv <- dv / sd(diff(ndmi))
  dn0 <- ifelse(is.na(dn), 0, dn); dv0 <- ifelse(is.na(dv), 0, dv)
  zn0 <- ifelse(is.na(zn), 0, zn); zv0 <- ifelse(is.na(zv), 0, zv)
  both   <- dn0 > 0 & dv0 > 0
  either <- dn0 > 0 | dv0 > 0
  dv_next <- c(dv0[-1], 0)
  lagconf <- dn0 > 0 & (dv0 > 0 | dv_next > 0)
  list(
    NBR        = dn0,
    NDMI       = dv0,
    AND_NBR    = ifelse(both, dn0, 0),
    AND_meanz  = ifelse(both, (zn0 + zv0) / 2, 0),
    AND_geomz  = ifelse(both, sqrt(zn0 * zv0), 0),
    ANDlag_NBR = ifelse(lagconf, dn0, 0),
    OR_meanz   = ifelse(either, (zn0 + zv0) / 2, 0),
    OR_maxz    = ifelse(either, pmax(zn0, zv0), 0))
}
rules <- c("NBR", "NDMI", "AND_NBR", "AND_meanz", "AND_geomz", "ANDlag_NBR", "OR_meanz", "OR_maxz")

cum_round <- function(d) diff(c(0, round(cumsum(d))))

allocate_hazard <- function(s, N0, dead) {
  s[is.na(s)] <- 0
  NT <- N0 - dead
  if (dead == 0) {
    rate <- rep(0, length(s))
  } else if (NT == 0) {                                # all trees died: proportional
    d <- dead * s / sum(s)
    Ns <- N0 - c(0, head(cumsum(d), -1))
    return(list(N_start = Ns, deaths = d))
  } else {
    rate <- 1 - exp(-log(N0 / NT) / sum(s) * s)
  }
  N_end <- N0 * cumprod(1 - rate)
  N_beg <- c(N0, head(N_end, -1))
  list(N_start = N_beg, deaths = N_beg - N_end)
}

# ---- Data -------------------------------------------------------------------
idx  <- as.data.frame(read_excel(f_index))
mort <- as.data.frame(read_excel(f_mort))
idx  <- idx[order(idx$Plot, idx$year), ]
stopifnot(setequal(idx$Plot, mort$Plot), !anyNA(idx[, c("NBR", "NDMI")]))

yrs_all <- sort(unique(idx$year))
sel     <- yrs_all >= year_start & yrs_all <= year_end
yrs     <- yrs_all[sel]

# Signals: per plot (Scheme A) and per site-mean series (Scheme B)
sig_plot <- lapply(split(idx, idx$Plot), function(g) {
  stopifnot(identical(g$year, yrs_all))
  lapply(consensus_signals(g$NBR, g$NDMI, window), function(s) s[sel])
})
site_mean <- aggregate(cbind(NBR, NDMI) ~ Site + year, idx, mean)
site_mean <- site_mean[order(site_mean$Site, site_mean$year), ]
sig_site <- lapply(split(site_mean, site_mean$Site), function(g)
  lapply(consensus_signals(g$NBR, g$NDMI, window), function(s) s[sel]))

# ---- Allocation: every plot x scheme x rule -----------------------------------
out <- list()
for (i in seq_len(nrow(mort))) {
  p <- mort$Plot[i]; st <- mort$Site[i]; N0 <- mort$TotalTrees[i]; dead <- mort$Dead[i]
  for (rule in rules) for (scheme in c("A", "B")) {
    s <- if (scheme == "A") sig_plot[[p]][[rule]] else sig_site[[st]][[rule]]
    fb <- FALSE
    if (dead > 0 && sum(s) <= 0) {
      s <- sig_site[[st]][[rule]]; fb <- TRUE
      if (sum(s) <= 0) s <- rep(1, length(s))        # last resort: uniform (flagged)
    }
    h <- allocate_hazard(s, N0, dead)
    d_int <- cum_round(h$deaths)
    out[[length(out) + 1]] <- data.frame(
      Rule = rule, Scheme = scheme, Method = paste(scheme, rule, sep = "_"),
      Site = st, Plot = p, Year = yrs, TotalTrees = N0, Dead = dead,
      signal = s, N_start = h$N_start, deaths = h$deaths, deaths_int = d_int,
      N_start_int = N0 - c(0, head(cumsum(d_int), -1)),
      mort_pct = 100 * h$deaths / h$N_start, fallback = fb)
  }
}
plot_long <- do.call(rbind, out)

# Totals preserved?
chk <- aggregate(deaths_int ~ Method + Plot, plot_long, sum)
chk <- merge(chk, mort[, c("Plot", "Dead")], by = "Plot")
stopifnot(all(chk$deaths_int == chk$Dead))
cat("All plot and site totals preserved for every rule and scheme.\n")

fb_tab <- aggregate(fallback ~ Method, unique(plot_long[, c("Method", "Plot", "fallback")]), sum)
names(fb_tab)[2] <- "plots_using_fallback"

# ---- Site series (sum of plots; tree-weighted %) ------------------------------
site_long <- aggregate(cbind(deaths, deaths_int, N_start, N_start_int) ~ Method + Rule + Scheme + Site + Year,
                       plot_long, sum)
site_long$mort_pct <- 100 * site_long$deaths / site_long$N_start
site_long <- site_long[order(site_long$Method, site_long$Site, site_long$Year), ]

# ---- Combined LONG tables: one column per method (A_/B_ x rule), % -------------

method_order <- as.vector(outer(c("A", "B"), rules, paste, sep = "_"))   # A_NBR, B_NBR, A_NDVI, ...

to_methods_wide <- function(df, id) {
  w <- reshape(df[, c(id, "Year", "Method", "mort_pct")], idvar = c(id, "Year"),
               timevar = "Method", direction = "wide")
  names(w) <- sub("^mort_pct\\.", "", names(w))
  w[do.call(order, w[, c(id, "Year")]), c(id, "Year", method_order)]
}

site_census <- aggregate(cbind(TotalTrees, Dead) ~ Site, mort, sum)
site_census$PercentDead_census <- 100 * site_census$Dead / site_census$TotalTrees
site_wide <- merge(site_census, to_methods_wide(site_long, "Site"), by = "Site")
site_wide <- site_wide[order(site_wide$Site, site_wide$Year), ]
write.csv(site_wide, file.path(out_dir, "site_annual_mortpct_CONSENSUS_long.csv"), row.names = FALSE)

plot_census <- mort[, c("Plot", "TotalTrees", "Dead")]
plot_census$PercentDead_census <- 100 * plot_census$Dead / plot_census$TotalTrees
plot_wide <- merge(plot_census, to_methods_wide(plot_long, c("Site", "Plot")), by = "Plot")
plot_wide <- plot_wide[order(plot_wide$Site, plot_wide$Plot, plot_wide$Year),
                       c("Site", setdiff(names(plot_wide), "Site"))]
write.csv(plot_wide, file.path(out_dir, "plot_annual_mortpct_CONSENSUS_long.csv"), row.names = FALSE)

write.csv(plot_long, file.path(out_dir, "plot_annual_CONSENSUS_full.csv"), row.names = FALSE)
write.csv(site_long, file.path(out_dir, "site_annual_CONSENSUS_full.csv"), row.names = FALSE)

# ---- Diagnostics ------------------------------------------------------------
site_of <- setNames(mort$Site, mort$Plot)
haz     <- -log(1 - mort$Dead / mort$TotalTrees)
refA    <- site_long[site_long$Method == "A_NBR", ]
diag <- do.call(rbind, lapply(rules, function(rule) {
  S  <- sapply(mort$Plot, function(p) sig_plot[[p]][[rule]])
  C  <- suppressWarnings(cor(S)); ut <- upper.tri(C)
  same <- outer(site_of[colnames(C)], site_of[colnames(C)], "==")
  A <- site_long[site_long$Method == paste0("A_", rule), ]
  B <- site_long[site_long$Method == paste0("B_", rule), ]
  data.frame(Rule = rule,
             signal_years_per_plot = mean(colSums(S > 0)),
             coherence_within_site = mean(C[ut & same], na.rm = TRUE),
             coherence_between     = mean(C[ut & !same], na.rm = TRUE),
             xplot_spearman        = cor(colSums(S), haz, method = "spearman"),
             site_R2_A_vs_B        = cor(A$deaths, B$deaths)^2,
             site_R2_vs_A_NBR      = cor(A$deaths, refA$deaths)^2)
}))
diag <- merge(diag, data.frame(Rule = sub("^A_", "", fb_tab$Method[startsWith(fb_tab$Method, "A_")]),
                               plots_fallback_A = fb_tab$plots_using_fallback[startsWith(fb_tab$Method, "A_")]),
              by = "Rule")
diag <- diag[match(rules, diag$Rule), ]
write.csv(diag, file.path(out_dir, "spain_consensus_diagnostics.csv"), row.names = FALSE)
cat("\nDiagnostics (no annual observed mortality exists for Spain):\n")
dd <- diag; dd[, 2:7] <- round(dd[, 2:7], 2); print(dd, row.names = FALSE)

cat("\nSite annual mortality (%), Scheme A, ANDlag_NBR:\n")
print(round(xtabs(mort_pct ~ Site + Year, site_long[site_long$Method == "A_ANDlag_NBR", ]), 2))

# ---- Figure -----------------------------------------------------------------
pl <- site_long[site_long$Scheme == "A", ]
pl$Rule <- factor(pl$Rule, levels = rules)
pf <- ggplot(pl, aes(Year, mort_pct, colour = Rule)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1) +
  facet_wrap(~ Site, scales = "free_y") +
  scale_x_continuous(breaks = yrs) +
  labs(x = NULL, y = "Annual mortality (%)", colour = NULL,
       title = "Spain: site annual mortality under NBR/NDMI agreement rules (Scheme A)") +
  theme_bw() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5))
pf
ggsave(file.path(out_dir, "spain_consensus_site_mortpct.png"), pf, width = 11, height = 7, dpi = 300)

cat("\nOutputs written to", normalizePath(out_dir), "\n")
