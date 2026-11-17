###############################################################################
# Event-based verification of simulated mortality against observations
# Criterion: a year is flagged when mortality exceeds the site mean,
#            applied independently to observed and simulated series.
###############################################################################

library(tidyverse)

# Resolve masking from MASS / plotly
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
    sim_ev  = SimMor > sim_thr,
    class   = ifelse( obs_ev &  sim_ev, "Hit",
                      ifelse(!obs_ev &  sim_ev, "False alarm",
                             ifelse( obs_ev & !sim_ev, "Miss", "Correct negative")))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Country, Site, Year)

###############################################################################
# 2. Sitewise contingency counts, skill scores and efficiency measures
###############################################################################
skill_sites <- ev %>%
  dplyr::group_by(Country, Site) %>%
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
      sum((ObsMor - mean(ObsMor_F))^2),
    NSE_z     = 1 - sum(((ObsMor_F - mean(ObsMor_F)) / sd(ObsMor_F) -
                           (SimMor - mean(SimMor)) / sd(SimMor))^2) / (dplyr::n() - 1),
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
    KGE      = 1 - sqrt((sqrt(R2) - 1)^2 + (KGE_alpha - 1)^2 + (KGE_beta - 1)^2),
    Fisher_p = purrr::pmap_dbl(
      list(Hits, FA, Miss, CN),
      function(a, b, cc, d)
        fisher.test(matrix(c(a, b, cc, d), nrow = 2, byrow = TRUE),
                    alternative = "greater")$p.value
    )
  ) %>%
  dplyr::relocate(Country, Site, Years, n, obs_thr, sim_thr,
                  n_obs_ev, n_sim_ev, Hits, FA, Miss, CN, Exp_hits)

print(skill_sites, width = Inf)

###############################################################################
# 3. Frequency table of outcome classes
###############################################################################
freq_sites <- ev %>%
  dplyr::count(Country, Site, class) %>%
  tidyr::pivot_wider(names_from = class, values_from = n, values_fill = 0) %>%
  dplyr::mutate(Total = rowSums(dplyr::across(where(is.numeric))))

print(freq_sites, width = Inf)

# Relative frequencies
freq_rel <- ev %>%
  dplyr::count(Country, Site, class) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(freq = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-n) %>%
  tidyr::pivot_wider(names_from = class, values_from = freq, values_fill = 0)

print(freq_rel, width = Inf)

###############################################################################
# 4. Pooled test across sites 
###############################################################################
# Pooled 2 x 2 counts
pooled <- skill_sites %>%
  dplyr::summarise(Hits = sum(Hits), FA = sum(FA),
                   Miss = sum(Miss), CN = sum(CN), n = sum(n))

pooled_tab <- matrix(c(pooled$Hits, pooled$FA, pooled$Miss, pooled$CN),
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
    POD = POD, FAR = FA / (Hits + FA), POFD = POFD,
    Accuracy = (Hits + CN) / n,
    PSS = POD - POFD,
    HSS = 2 * (Hits * CN - FA * Miss) / den
  )
})

print(pooled_scores, width = Inf)

# Mantel-Haenszel: stratified by site, respects between-site differences
mh_array <- table(Observed  = ev$obs_ev,
                  Simulated = ev$sim_ev,
                  Site      = ev$Site)

mantelhaen.test(mh_array, exact = TRUE, alternative = "greater")

# Site-stratified logistic regression as a continuous alternative
glm_fit <- glm(obs_ev ~ scale(SimMor) + Site, data = ev, family = binomial)
summary(glm_fit)

###############################################################################
# 5. Year-by-year classification (appendix table)
###############################################################################
event_table <- ev %>%
  dplyr::select(Country, Site, Year, ObsMor, SimMor,
                obs_thr, sim_thr, obs_ev, sim_ev, class)

print(event_table, n = Inf, width = Inf)

# Only the flagged years
event_table %>%
  dplyr::filter(class %in% c("Hit", "Miss", "False alarm")) %>%
  print(n = Inf, width = Inf)

###############################################################################
# 6. Figure: series, site means and classification
###############################################################################
plot_df <- ev %>%
  dplyr::select(Site, Year, ObsMor, SimMor, obs_thr, sim_thr, class) %>%
  tidyr::pivot_longer(c(ObsMor, SimMor),
                      names_to = "Series", values_to = "Mortality") %>%
  dplyr::mutate(
    thr    = ifelse(Series == "ObsMor", obs_thr, sim_thr),
    Series = dplyr::recode(Series, ObsMor = "Observed",
                           SimMor = "Simulated (iLand)")
  )

ggplot(plot_df, aes(Year, Mortality)) +
  geom_hline(aes(yintercept = thr), linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  geom_line(colour = "grey60", linewidth = 0.4) +
  geom_point(aes(colour = class), size = 1.9) +
  facet_grid(Series ~ Site, scales = "free_y") +
  scale_colour_manual(values = c("Hit"              = "#1A9850",
                                 "Miss"             = "#D73027",
                                 "False alarm"      = "#FDAE61",
                                 "Correct negative" = "grey70")) +
  labs(x = NULL, y = "Annual mortality (%)", colour = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"))

###############################################################################
# 7. Export
###############################################################################
# readr::write_csv(skill_sites,  "skill_sites_abovemean.csv")
# readr::write_csv(freq_sites,   "class_frequencies_abovemean.csv")
# readr::write_csv(event_table,  "event_classification_abovemean.csv")