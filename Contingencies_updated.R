###############################################################################
# Event-based verification of simulated against observed mortality
# Observed series: ObsMor_F
# Criterion: a year is flagged when mortality exceeds the site mean,
#            applied independently to observed and simulated series.
###############################################################################

library(tidyverse)
library(conflicted)
conflicts_prefer(dplyr::select, dplyr::filter, dplyr::lag, .quiet = TRUE)

data <- readxl::read_xlsx(
  "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/ObsMor_anomaly.xlsx"
)

###############################################################################
# 1. Classify every site-year
###############################################################################
ev <- data %>%
  tidyr::drop_na(ObsMor_F, SimMor) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(
    obs_thr = mean(ObsMor_F),
    sim_thr = mean(SimMor),
    obs_ev  = ObsMor_F > obs_thr,
    sim_ev  = SimMor   > sim_thr,
    class   = ifelse( obs_ev &  sim_ev, "Hit",
                      ifelse(!obs_ev &  sim_ev, "False alarm",
                             ifelse( obs_ev & !sim_ev, "Miss", "Correct negative")))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Location = gsub("[0-9]+$", "", Site)) %>%   # Arcalis1-4 share a location
  dplyr::arrange(Country, Site, Year)

# Sanity check: displayed column must match the classified column
stopifnot(all(ev$obs_ev == (ev$ObsMor_F > ev$obs_thr)))

###############################################################################
# 2. Sitewise contingency counts, skill scores and efficiency measures
###############################################################################
skill_sites <- ev %>%
  dplyr::group_by(Country, Location, Site) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    Years     = paste(min(Year), max(Year), sep = "-"),
    obs_thr   = dplyr::first(obs_thr),
    sim_thr   = dplyr::first(sim_thr),
    n_obs_ev  = sum(obs_ev),
    n_sim_ev  = sum(sim_ev),
    
    Hits      = sum( obs_ev &  sim_ev),
    FA        = sum(!obs_ev &  sim_ev),
    Miss      = sum( obs_ev & !sim_ev),
    CN        = sum(!obs_ev & !sim_ev),
    
    NSE       = 1 - sum((ObsMor_F - SimMor)^2) /
      sum((ObsMor_F - mean(ObsMor_F))^2),
    NSE_log   = 1 - sum((log1p(ObsMor_F) - log1p(SimMor))^2) /
      sum((log1p(ObsMor_F) - mean(log1p(ObsMor_F)))^2),
    NSE_z     = 1 - sum(((ObsMor_F - mean(ObsMor_F)) / sd(ObsMor_F) -
                           (SimMor - mean(SimMor)) / sd(SimMor))^2) / (dplyr::n() - 1),
    r_Pear    = cor(ObsMor_F, SimMor),
    R2        = cor(ObsMor_F, SimMor)^2,
    r_Spear   = cor(ObsMor_F, SimMor, method = "spearman"),
    Bias      = mean(SimMor - ObsMor_F),
    KGE_alpha = sd(SimMor) / sd(ObsMor_F),
    KGE_beta  = mean(SimMor) / mean(ObsMor_F),
    .groups   = "drop"
  ) %>%
  dplyr::mutate(
    Exp_hits = (Hits + Miss) * (Hits + FA) / n,
    POD      = Hits / (Hits + Miss),
    FAR      = FA / (Hits + FA),
    POFD     = FA / (FA + CN),
    Accuracy = (Hits + CN) / n,
    PSS      = POD - POFD,
    HSS      = 2 * (Hits * CN - FA * Miss) /
      ((Hits + Miss) * (Miss + CN) + (Hits + FA) * (FA + CN)),
    KGE      = 1 - sqrt((r_Pear - 1)^2 + (KGE_alpha - 1)^2 + (KGE_beta - 1)^2),
    OR_Haldane = ((Hits + 0.5) * (CN + 0.5)) / ((FA + 0.5) * (Miss + 0.5)),
    Fisher_p = purrr::pmap_dbl(
      list(Hits, FA, Miss, CN),
      function(a, b, cc, d)
        fisher.test(matrix(c(a, b, cc, d), nrow = 2, byrow = TRUE),
                    alternative = "greater")$p.value
    )
  ) %>%
  dplyr::relocate(Country, Location, Site, Years, n, obs_thr, sim_thr,
                  n_obs_ev, n_sim_ev, Hits, FA, Miss, CN, Exp_hits)

print(skill_sites, width = Inf)

###############################################################################
# 3. Frequency tables of outcome classes
###############################################################################
freq_sites <- ev %>%
  dplyr::count(Country, Site, class) %>%
  tidyr::pivot_wider(names_from = class, values_from = n, values_fill = 0)

freq_sites <- freq_sites %>%
  dplyr::mutate(Total = rowSums(dplyr::select(., -Country, -Site)))

print(freq_sites, width = Inf)

freq_rel <- ev %>%
  dplyr::count(Country, Site, class) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(freq = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-n) %>%
  tidyr::pivot_wider(names_from = class, values_from = freq, values_fill = 0)

print(freq_rel, width = Inf)

###############################################################################
# 4. Pooled analysis across all sites
###############################################################################
pooled <- skill_sites %>%
  dplyr::summarise(Hits = sum(Hits), FA = sum(FA),
                   Miss = sum(Miss), CN = sum(CN), n = sum(n))

pooled_tab <- matrix(c(pooled$Hits, pooled$Miss, pooled$FA, pooled$CN),
                     nrow = 2, byrow = TRUE,
                     dimnames = list(Observed  = c("Above mean", "Below mean"),
                                     Simulated = c("Above mean", "Below mean")))
pooled_tab
fisher.test(pooled_tab, alternative = "greater")

pooled_scores <- with(pooled, {
  POD  <- Hits / (Hits + Miss)
  POFD <- FA / (FA + CN)
  den  <- (Hits + Miss) * (Miss + CN) + (Hits + FA) * (FA + CN)
  tibble::tibble(
    n = n, Hits = Hits, FA = FA, Miss = Miss, CN = CN,
    Exp_hits = (Hits + Miss) * (Hits + FA) / n,
    POD = POD, FAR = FA / (Hits + FA), POFD = POFD,
    Accuracy = (Hits + CN) / n,
    PSS = POD - POFD,
    HSS = 2 * (Hits * CN - FA * Miss) / den
  )
})
print(pooled_scores, width = Inf)

# Exact Mantel-Haenszel, stratified by site
mh_array <- table(Observed  = ev$obs_ev,
                  Simulated = ev$sim_ev,
                  Site      = ev$Site)
mantelhaen.test(mh_array, exact = TRUE, alternative = "greater")

# Stratified by LOCATION instead, collapsing the four Arcalis plots
mh_array_loc <- table(Observed  = ev$obs_ev,
                      Simulated = ev$sim_ev,
                      Location  = ev$Location)
mantelhaen.test(mh_array_loc, exact = TRUE, alternative = "greater")

###############################################################################
# 4b. Homogeneity of the odds ratio across sites
###############################################################################
if (requireNamespace("DescTools", quietly = TRUE)) {
  print(DescTools::BreslowDayTest(mh_array))
} else {
  message("DescTools not installed; using Woolf test below")
}

# Woolf test, Haldane-corrected (works with zero cells)
woolf_manual <- skill_sites %>%
  dplyr::mutate(
    logOR = log(OR_Haldane),
    w     = 1 / (1/(Hits + 0.5) + 1/(FA + 0.5) + 1/(Miss + 0.5) + 1/(CN + 0.5))
  )

logOR_bar <- with(woolf_manual, sum(w * logOR) / sum(w))
Q_stat    <- with(woolf_manual, sum(w * (logOR - logOR_bar)^2))
df_Q      <- nrow(woolf_manual) - 1

c(pooled_OR = exp(logOR_bar), Q = Q_stat, df = df_Q,
  p_heterogeneity = pchisq(Q_stat, df_Q, lower.tail = FALSE))

###############################################################################
# 4c. Regression alternatives
###############################################################################
glm_fit <- glm(obs_ev ~ scale(SimMor) + Site, data = ev, family = binomial)
summary(glm_fit)

# Mixed model: site nested within location, to respect the Arcalis plots
if (requireNamespace("lme4", quietly = TRUE)) {
  glmer_fit <- lme4::glmer(obs_ev ~ scale(SimMor) + (1 | Location / Site),
                           data = ev, family = binomial)
  print(summary(glmer_fit))
}

###############################################################################
# 5. Circular-shift permutation test (controls serial autocorrelation)
###############################################################################
set.seed(42)
n_perm <- 9999

circ_shift <- function(x, k) {
  if (k == 0) return(x)
  c(x[-seq_len(k)], x[seq_len(k)])
}

perm_hits <- function(dd, n_perm = 9999) {
  sp <- split(dd, dd$Site)
  o  <- lapply(sp, function(s) s$obs_ev)
  s  <- lapply(sp, function(s) s$sim_ev)
  hits_obs <- sum(mapply(function(a, b) sum(a & b), o, s))
  null_hits <- replicate(n_perm, {
    sum(mapply(function(a, b) {
      k <- sample.int(length(b), 1) - 1
      sum(a & circ_shift(b, k))
    }, o, s))
  })
  tibble::tibble(
    Hits      = hits_obs,
    null_mean = mean(null_hits),
    null_sd   = sd(null_hits),
    null_q95  = unname(quantile(null_hits, 0.95)),
    p_shift   = (sum(null_hits >= hits_obs) + 1) / (n_perm + 1)
  )
}

perm_all <- perm_hits(ev, n_perm)
print(perm_all, width = Inf)

###############################################################################
# 6. Country-level analysis
###############################################################################
skill_country <- ev %>%
  dplyr::group_by(Country) %>%
  dplyr::summarise(
    n_sites  = dplyr::n_distinct(Site),
    n        = dplyr::n(),
    Years    = paste(min(Year), max(Year), sep = "-"),
    n_obs_ev = sum(obs_ev),
    n_sim_ev = sum(sim_ev),
    Hits     = sum( obs_ev &  sim_ev),
    FA       = sum(!obs_ev &  sim_ev),
    Miss     = sum( obs_ev & !sim_ev),
    CN       = sum(!obs_ev & !sim_ev),
    NSE_pool = 1 - sum((ObsMor_F - SimMor)^2) /
      sum((ObsMor_F - mean(ObsMor_F))^2),
    R2_pool  = cor(ObsMor_F, SimMor)^2,
    r_Spear  = cor(ObsMor_F, SimMor, method = "spearman"),
    Bias     = mean(SimMor - ObsMor_F),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(
    Exp_hits = (Hits + Miss) * (Hits + FA) / n,
    POD      = Hits / (Hits + Miss),
    FAR      = FA / (Hits + FA),
    POFD     = FA / (FA + CN),
    Accuracy = (Hits + CN) / n,
    PSS      = POD - POFD,
    HSS      = 2 * (Hits * CN - FA * Miss) /
      ((Hits + Miss) * (Miss + CN) + (Hits + FA) * (FA + CN)),
    ft       = purrr::pmap(list(Hits, FA, Miss, CN),
                           function(a, b, cc, d)
                             fisher.test(matrix(c(a, b, cc, d), nrow = 2, byrow = TRUE),
                                         alternative = "greater")),
    OR       = purrr::map_dbl(ft, ~ unname(.x$estimate)),
    OR_lo95  = purrr::map_dbl(ft, ~ .x$conf.int[1]),
    Fisher_p = purrr::map_dbl(ft, ~ .x$p.value)
  ) %>%
  dplyr::select(-ft)

# Sitewise NSE averaged within country, not dominated by the largest-variance site
nse_by_country <- skill_sites %>%
  dplyr::group_by(Country) %>%
  dplyr::summarise(NSE_site_mean   = mean(NSE),
                   NSE_site_median = median(NSE),
                   .groups = "drop")

skill_country <- skill_country %>%
  dplyr::left_join(nse_by_country, by = "Country")

print(skill_country, width = Inf)

# Site-stratified exact test within each country
mh_country <- ev %>%
  dplyr::group_split(Country) %>%
  purrr::map_dfr(function(dd) {
    cname  <- unique(dd$Country)
    nsites <- dplyr::n_distinct(dd$Site)
    res <- tryCatch({
      if (nsites > 1) {
        mantelhaen.test(table(dd$obs_ev, dd$sim_ev, dd$Site),
                        exact = TRUE, alternative = "greater")
      } else {
        fisher.test(table(dd$obs_ev, dd$sim_ev), alternative = "greater")
      }
    }, error = function(e) NULL)
    tibble::tibble(
      Country = cname, n_sites = nsites,
      test    = ifelse(nsites > 1, "Exact Mantel-Haenszel", "Fisher exact"),
      OR      = if (is.null(res)) NA_real_ else unname(res$estimate),
      OR_lo95 = if (is.null(res)) NA_real_ else res$conf.int[1],
      p_value = if (is.null(res)) NA_real_ else res$p.value
    )
  })
print(mh_country, width = Inf)

# Circular-shift permutation within each country
perm_country <- ev %>%
  dplyr::group_split(Country) %>%
  purrr::map_dfr(function(dd)
    dplyr::bind_cols(tibble::tibble(Country = unique(dd$Country)),
                     perm_hits(dd, n_perm)))
print(perm_country, width = Inf)

# Combined country reporting table
report_country <- skill_country %>%
  dplyr::select(Country, n_sites, n, n_obs_ev, Hits, FA, Miss, CN, Exp_hits,
                POD, FAR, POFD, PSS, HSS, Accuracy,
                NSE_pool, NSE_site_mean, R2_pool, Bias) %>%
  dplyr::left_join(dplyr::select(mh_country, Country, test, OR, OR_lo95, p_value),
                   by = "Country") %>%
  dplyr::left_join(dplyr::select(perm_country, Country, p_shift), by = "Country") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))

print(report_country, width = Inf)

###############################################################################
# 7. Year-by-year classification (appendix table)
###############################################################################
event_table <- ev %>%
  dplyr::select(Country, Site, Year, ObsMor_F, SimMor,
                obs_thr, sim_thr, obs_ev, sim_ev, class)

print(event_table, n = Inf, width = Inf)

event_table %>%
  dplyr::filter(class %in% c("Hit", "Miss", "False alarm")) %>%
  print(n = Inf, width = Inf)

###############################################################################
# 8. Figures
###############################################################################
plot_df <- ev %>%
  dplyr::select(Site, Year, ObsMor_F, SimMor, obs_thr, sim_thr, class) %>%
  tidyr::pivot_longer(c(ObsMor_F, SimMor),
                      names_to = "Series", values_to = "Mortality") %>%
  dplyr::mutate(
    thr    = ifelse(Series == "ObsMor_F", obs_thr, sim_thr),
    Series = dplyr::recode(Series, ObsMor_F = "Observed",
                           SimMor   = "Simulated (iLand)")
  )

p_series <- ggplot(plot_df, aes(Year, Mortality)) +
  geom_hline(aes(yintercept = thr), linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  geom_line(colour = "grey60", linewidth = 0.4) +
  geom_point(aes(colour = class), size = 1.9) +
  facet_grid(Series ~ Site, scales = "free") +
  scale_colour_manual(values = c("Hit"              = "#1A9850",
                                 "Miss"             = "#D73027",
                                 "False alarm"      = "#FDAE61",
                                 "Correct negative" = "grey70")) +
  labs(x = NULL, y = "Annual mortality (%)", colour = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"))

print(p_series)

# Forest plot of sitewise odds ratios
forest_df <- skill_sites %>%
  dplyr::transmute(Site, Country, OR = OR_Haldane) %>%
  dplyr::mutate(Site = forcats::fct_reorder(Site, OR))

p_forest <- ggplot(forest_df, aes(OR, Site, colour = Country)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.6) +
  scale_x_log10() +
  labs(x = "Odds ratio (Haldane-corrected, log scale)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

print(p_forest)

###############################################################################
# 9. Export
###############################################################################
# readr::write_csv(skill_sites,     "skill_sites_ObsMorF.csv")
# readr::write_csv(freq_sites,      "class_frequencies_ObsMorF.csv")
# readr::write_csv(report_country,  "verification_by_country_ObsMorF.csv")
# readr::write_csv(event_table,     "event_classification_ObsMorF.csv")
# ggsave("fig_series_classification.png", p_series, width = 12, height = 5, dpi = 300)