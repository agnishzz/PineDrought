library(dplyr)
library(tidyr)
library(readxl)
library(writexl)
library(ggplot2)
library(corrplot)
library(ggstats)
library(trend)
library(glmmTMB)
library(scales)
library(broom)
library(ggpmisc)


vol <- readxl::read_excel("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/IGA_v2/ExcelSheets2/CumulativeVolume_both.xlsx")
vol <- vol %>%
  mutate(Simulation_site = case_when(
    Simulation_site == "Visp" ~ "CH_Visp",
    Simulation_site == "Simsite0" ~ "SP_Prades",
    Simulation_site == "Simsite1" ~ "SP_Arcalis1",
    Simulation_site == "Simsite2" ~ "SP_Arcalis2",
    Simulation_site == "Simsite3" ~ "SP_Arcalis3",
    Simulation_site == "Simsite4" ~ "SP_Arcalis4",
    TRUE ~ Simulation_site
  ))

head(vol)

vol_f <- vol %>% filter(Year == 2100)
head(vol_f)
mor <- readxl::read_excel("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/IGA_v2/ExcelSheets2/CumulativeMortality_new.xlsx")
head(mor)
mor <- mor %>% filter(Year == 2100)
head(vol_f)
head(mor)
vol_f_new <- vol_f %>%
  left_join(
    mor %>% dplyr::select(Year, Climate_model, Scenario, Country, Simulation_site, Mor),
    by = c("Year", "Climate_model", "Scenario", "Country", "Simulation_site")
  )
head(vol_f_new)

vol_diff_explicit <- vol_f_new %>%
  group_by(Simulation_site) %>%
  mutate(
    # Calculate for the first variable
    perc_diff_cum_VO_tr = ((Mor - Mor[Scenario == "Refclim"]) / Mor[Scenario == "Refclim"]) * 100,
    
    # Repeat the logic for the second variable
    perc_diff_VO = ((VO - VO[Scenario == "Refclim"]) / VO[Scenario == "Refclim"]) * 100
  ) %>%
  filter(Scenario != "Refclim") %>%
  ungroup()

print(vol_diff_explicit)

##Colours for each RCP
rcp_colors <- c(
  "RCP 2.6" = "blue4",
  "RCP 4.5" = "orange",
  "RCP 8.5" = "brown4"
)


# 1. Calculate the mean values for the "Refclim" scenario for each site-model.
# We will append "_ref" to the column names to mark them as the reference values.
clim_sp <- readxl::read_excel("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/ClimateSummaries_Spain.xlsx")
clim_sw <- readxl::read_excel("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/ClimateSummaries_Swiss.xlsx")
head(clim_sp)
head(clim_sw)
clim_sp %>% 
  count(Scenario)


refclim_means_site <- clim_sp %>%
  filter(Scenario == "Refclim") %>%
  group_by(Simulation_site) %>%  # <-- CHANGE HERE
  summarise(
    across(TMean_WY:Wdc, \(x) mean(x, na.rm = TRUE), .names = "{.col}_ref"),
    .groups = "drop"
  )

# 2. Calculate the mean values for the future RCP scenarios (this part is unchanged).
future_means <- clim_sp %>%
  filter(Scenario != "Refclim") %>%
  group_by(Simulation_site, Climate_model, Scenario) %>%
  summarise(
    across(TMean_WY:Wdc, \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

# 3. Join the reference means onto the future means using ONLY Simulation_site.
# This will apply the same site-level reference to all climate models for that site.
delta_values_sp <- future_means %>%
  left_join(refclim_means_site, by = "Simulation_site") %>% # <-- CHANGE HERE
  
  # 4. Calculate the deltas (this part is unchanged and will now work).
  mutate(
    TMean_WY = TMean_WY - TMean_WY_ref,
    Prec_WY = Prec_WY - Prec_WY_ref,
    VPD_WY = VPD_WY - VPD_WY_ref,
    CWB_WY = CWB_WY - CWB_WY_ref,
    TMean_MSEA = TMean_MSEA - TMean_MSEA_ref,
    Prec_MSEA = Prec_MSEA - Prec_MSEA_ref,
    VPD_MSEA = VPD_MSEA - VPD_MSEA_ref,
    CWB_MSEA = CWB_MSEA - CWB_MSEA_ref,
    CDD = CDD - CDD_ref,
    SDII = SDII - SDII_ref,
    Wdc = Wdc - Wdc_ref
  ) %>%
  
  # 5. Clean up (unchanged).
  dplyr:: select(-ends_with("_ref"))

# View the final, correct result
print(delta_values_sp)



refclim_means_site_sw <- clim_sw %>%
  filter(Scenario == "Refclim") %>%
  group_by(Simulation_site) %>%  # <-- CHANGE HERE
  summarise(
    across(TMean_WY:Wdc, \(x) mean(x, na.rm = TRUE), .names = "{.col}_ref"),
    .groups = "drop"
  )

# 2. Calculate the mean values for the future RCP scenarios (this part is unchanged).
future_means_sw <- clim_sw %>%
  filter(Scenario != "Refclim") %>%
  group_by(Simulation_site, Climate_model, Scenario) %>%
  summarise(
    across(TMean_WY:Wdc, \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

# 3. Join the reference means onto the future means using ONLY Simulation_site.
# This will apply the same site-level reference to all climate models for that site.
delta_values_sw <- future_means_sw %>%
  left_join(refclim_means_site_sw, by = "Simulation_site") %>% # <-- CHANGE HERE
  
  # 4. Calculate the deltas (this part is unchanged and will now work).
  mutate(
    TMean_WY = TMean_WY - TMean_WY_ref,
    Prec_WY = Prec_WY - Prec_WY_ref,
    VPD_WY = VPD_WY - VPD_WY_ref,
    CWB_WY = CWB_WY - CWB_WY_ref,
    TMean_MSEA = TMean_MSEA - TMean_MSEA_ref,
    Prec_MSEA = Prec_MSEA - Prec_MSEA_ref,
    VPD_MSEA = VPD_MSEA - VPD_MSEA_ref,
    CWB_MSEA = CWB_MSEA - CWB_MSEA_ref,
    CDD = CDD - CDD_ref,
    SDII = SDII - SDII_ref,
    Wdc = Wdc - Wdc_ref
  ) %>%
  
  # 5. Clean up (unchanged).
  dplyr:: select(-ends_with("_ref"))

# View the final, correct result
print(delta_values_sw)

delta_values_sp
delta_values_sw
delta_values_sp  <- delta_values_sp  %>% mutate(Country = "Spain")
delta_values_sw  <- delta_values_sw  %>% mutate(Country = "Switzerland")

delta_values <- bind_rows(delta_values_sp, delta_values_sw)
head(delta_values)
head(vol_diff_explicit)

vol_diff_explicit_merged <- vol_diff_explicit %>%
  left_join(
    delta_values,
    by = c("Simulation_site", "Climate_model", "Scenario", "Country")
  )

head(vol_diff_explicit_merged)

###################
#Plot for the significant variables, precipitation, CWB, Wet day count#
####################
head(vol_diff_explicit_merged)

# Fit the linear model for the entire dataset
model <- lm(perc_diff_cum_VO_tr ~ Prec_MSEA, data = vol_diff_explicit_merged)

# Extract the values you need
model_summary <- summary(model)
slope <- coef(model)[2]  # This will be a negative number (e.g., -0.51)
r_squared <- model_summary$r.squared
p_value <- model_summary$coefficients[2, 4]

# Format the p-value
p_label <- if (p_value < 0.001) {
  "p < 0.001"
} else if (p_value < 0.01) {
  "p < 0.01"
} else if (p_value < 0.05) {
  "p < 0.05"
} else {
  paste("p =", round(p_value, 3))
}

### --- STEP 2: CREATE THE NEW ANNOTATION TEXT (This is the changed part) ---

# We use abs(slope) to get the positive magnitude of the effect
# and change the text to reflect the "drying effect".
annotation_text <- paste(
  # The text now describes mortality *increase* per mm of *drying*
  #sprintf("Drying Effect: %.2f%% mortality per mm", abs(slope)),
  sprintf("R² = %.2f", r_squared),
  p_label,
  sep = "\n" # Use \n to separate lines
)

# Check the new text
print(annotation_text)
#> [1] "Drying Effect: +0.51% mortality per mm\nR² = 0.43\np < 0.001"


# First, ensure your color palette is defined
rcp_colors <- c(
  "RCP 2.6" = "blue4",
  "RCP 4.5" = "orange",
  "RCP 8.5" = "brown4"
)
head(vol_diff_explicit_merged)
# Now, create the plot
cc_prec_msea_m1 <- ggplot(
  data = vol_diff_explicit_merged, 
  # Base aesthetics for x and y are defined for all layers
  aes(x = Prec_MSEA, y = perc_diff_cum_VO_tr) 
) +
  # --- CHANGE 1: Map the color aesthetic to your data column ---
  # I'm assuming your column is named 'scenario'. Change if necessary.
  geom_point(aes(color = Scenario, shape = Simulation_site), size = 3, alpha = 1) + 
  
  
  # The smooth line remains for the whole dataset as it doesn't inherit the color mapping
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, color = "black") +
  
  
  # Annotation remains the same
  annotate(
    geom = "text",
    x = Inf,
    y = Inf,
    label = annotation_text, # Your custom text string
    hjust = 1.05,
    vjust = 1.2,
    size = 4,
    fontface = "italic"
  ) +
  
  # --- CHANGE 2: Use scale_color_manual to apply your custom colors ---
  scale_color_manual(values = rcp_colors) + 
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    y = "Mortality difference (%)",
    x = "Change in precipitation (mm)",
    # --- CHANGE 3: Add a title for the color legend ---
    color = "RCP Scenario", 
    #caption = "Slope expressed as % mortality change per mm precipitation difference",
    shape = "Site"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y =  element_text(face = "bold", size = 12),
    axis.title.x =  element_text(face = "bold", size = 12),
    # --- CHANGE 4: Make the legend visible ---
    #legend.position = c(0.98, 0.57),
    legend.position = "none",
    
    # Justify the legend so its top-right corner is at the specified position
    legend.justification = c("right"),
    panel.grid.major = element_blank(),
    legend.background = element_blank(),
    panel.grid.minor = element_blank()
  )

# Display the plot
print(cc_prec_msea_m1)




model1 <- lm(perc_diff_VO ~ Prec_MSEA, data = vol_diff_explicit_merged)

# Extract the values you need
model_summary1 <- summary(model1)
slope1 <- coef(model1)[2]  # This will be a negative number (e.g., -0.51)
r_squared1 <- model_summary1$r.squared
p_value1 <- model_summary1$coefficients[2, 4]

# Format the p-value
p_label1 <- if (p_value1 < 0.001) {
  "p < 0.001"
} else if (p_value1 < 0.01) {
  "p < 0.01"
} else if (p_value1 < 0.05) {
  "p < 0.05"
} else {
  paste("p =", round(p_value1, 3))
}

### --- STEP 2: CREATE THE NEW ANNOTATION TEXT (This is the changed part) ---

# We use abs(slope) to get the positive magnitude of the effect
# and change the text to reflect the "drying effect".
annotation_text1 <- paste(
  # The text now describes growing per mm of *drying*
  #sprintf("Drying Effect: -%.2f%% growing stock per mm", slope1),
  sprintf("R² = %.2f", r_squared1),
  p_label1,
  sep = "\n" # Use \n to separate lines
)

# Check the new text
print(annotation_text1)

cc_prec_msea_v1 <- ggplot(
  data = vol_diff_explicit_merged, 
  # Remove color and group from aes() to plot all data together
  aes(x = Prec_MSEA, y = perc_diff_VO) 
) +
  # Use a single color for points, specified outside aes()
  geom_point(aes(color = Scenario, shape = Simulation_site), size = 3, alpha = 1) + 
  # The smooth line will now be for the whole dataset
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, color = "black") +
  
  
  # Add your custom, pre-calculated text using annotate()
  annotate(
    geom = "text",
    x = -30,        # Position at the far left
    y = Inf,         # Position at the very top
    label = annotation_text1, # Your custom text string
    hjust = 0,       # Right-align the text
    vjust = 1.2,     # Nudge it down from the top edge
    size = 4,
    fontface = "italic"
  ) +
  
  scale_color_manual(values = rcp_colors) + 
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    #subtitle = "Effect of climate change quantified by difference in Precipitation in Midsummer to Early Autumn",
    y = "Growing stock difference (%)",
    x = "Change in precipitation (mm)",
    # The color legend title is no longer needed
    # color = "RCP Scenario", 
    #caption = "Slope expressed as % growing stock change per mm precipitation difference"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y =  element_text(face = "bold", size = 12),
    axis.title.x =  element_text(face = "bold", size = 12),
    # --- CHANGE 4: Make the legend visible ---
    #legend.position = c(0.98, 0.57),
    legend.position = "none",
    
    # Justify the legend so its top-right corner is at the specified position
    legend.justification = c("right"),
    panel.grid.major = element_blank(),
    legend.background = element_blank(),
    panel.grid.minor = element_blank()
  )

# Display the plot
print(cc_prec_msea_v1)

ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Figures/Prec_Mor_4.pdf", plot = cc_prec_msea_m1, width = 6, height = 8, dpi = 300)
ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Figures/Prec_GS_4.pdf", plot = cc_prec_msea_v1, width = 6, height = 8, dpi = 300)
######################
#Repeat the above for climatic water balance#

#--- STEP 1: FIT THE LINEAR MODEL AND EXTRACT STATISTICS ---

# Fit the linear model for the entire dataset using CWB_MSEA
model_cwb_m <- lm(perc_diff_cum_VO_tr ~ CWB_MSEA, data = vol_diff_explicit_merged)

# Extract the values you need
model_summary_cwb_m <- summary(model_cwb_m)
slope_cwb_m <- coef(model_cwb_m)[2]
r_squared_cwb_m <- model_summary_cwb_m$r.squared
p_value_cwb_m <- model_summary_cwb_m$coefficients[2, 4]

# Format the p-value
p_label_cwb_m <- if (p_value_cwb_m < 0.001) {
  "p < 0.001"
}  else if (p_label_cwb_m < 0.01) {
  "p < 0.01"
} else if (p_label_cwb_m < 0.05) {
  "p < 0.05"
} else {
  paste("p =", round(p_label_cwb_m, 3))
}

# --- STEP 2: CREATE THE ANNOTATION TEXT ---

# The text still reflects a "drying effect" as a negative CWB is a water deficit.
annotation_text_cwb_m <- paste(
  #sprintf("Drying Effect: %.2f%% mortality per mm", abs(slope_cwb_m)),
  sprintf("R² = %.2f", r_squared_cwb_m),
  p_label_cwb_m,
  sep = "\n" # Use \n to separate lines
)

# Check the new text
print("Annotation for Mortality vs CWB plot:")
print(annotation_text_cwb_m)

max(vol_diff_explicit_merged$CWB_MSEA)
min(vol_diff_explicit_merged$CWB_MSEA)
# --- STEP 3: CREATE THE PLOT ---

# First, ensure your color palette is defined
rcp_colors <- c(
  "RCP 2.6" = "blue4",
  "RCP 4.5" = "orange",
  "RCP 8.5" = "brown4"
)

# Now, create the plot using CWB_MSEA on the x-axis
cc_cwb_msea_m1 <- ggplot(
  data = vol_diff_explicit_merged, 
  aes(x = CWB_MSEA, y = perc_diff_cum_VO_tr) 
) +
  geom_point(aes(color = Scenario, shape = Simulation_site), size = 3, alpha = 1) + 
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, color = "black") +
  scale_x_continuous(breaks = seq(-90, 20, by = 20),
                     limits = c(-90, 20))+
  # Add the new annotation text
  annotate(
    geom = "text",
    x = Inf,
    y = Inf,
    label = annotation_text_cwb_m, # Use the CWB-specific text
    hjust = 1.05,
    vjust = 1.2,
    size = 4,
    fontface = "italic"
  ) +
  
  scale_color_manual(values = rcp_colors) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    y = "Mortality difference (%)",
    # Update the x-axis label for CWB
    x = "Change in climatic water balance (mm)",
    color = "RCP Scenario", 
    # Update the caption for CWB
    #caption = "Slope expressed as % mortality change per mm CWB difference",
    shape = "Site"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y =  element_text(face = "bold", size = 12),
    axis.title.x =  element_text(face = "bold", size = 12),
    # --- CHANGE 4: Make the legend visible ---
    #legend.position = c(0.98, 0.57),
    legend.position = "none",
    
    # Justify the legend so its top-right corner is at the specified position
    legend.justification = c("right"),
    panel.grid.major = element_blank(),
    legend.background = element_blank(),
    panel.grid.minor = element_blank()
  )

# Display the plot
print(cc_cwb_msea_m1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  PART 2: PLOT FOR GROWING STOCK vs. CWB_MSEA                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# --- STEP 1: FIT THE LINEAR MODEL AND EXTRACT STATISTICS ---

# Fit the linear model using CWB_MSEA
model_cwb_v <- lm(perc_diff_VO ~ CWB_MSEA, data = vol_diff_explicit_merged)

# Extract the values you need
model_summary_cwb_v <- summary(model_cwb_v)
slope_cwb_v <- coef(model_cwb_v)[2]
r_squared_cwb_v <- model_summary_cwb_v$r.squared
p_value_cwb_v <- model_summary_cwb_v$coefficients[2, 4]

# Format the p-value
p_label_cwb_v <- if (p_value_cwb_v < 0.001) {
  "p < 0.001"
}  else if (p_label_cwb_v < 0.01) {
  "p < 0.01"
} else if (p_label_cwb_v < 0.05) {
  "p < 0.05"
} else {
  paste("p =", round(p_label_cwb_v, 3))
}

# --- STEP 2: CREATE THE ANNOTATION TEXT ---

annotation_text_cwb_v <- paste(
  #sprintf("Drying Effect: -%.2f%% growing stock per mm", slope_cwb_v),
  sprintf("R² = %.2f", r_squared_cwb_v),
  p_label_cwb_v,
  sep = "\n" # Use \n to separate lines
)

# Check the new text
print("Annotation for Growing Stock vs CWB plot:")
print(annotation_text_cwb_v)


# --- STEP 3: CREATE THE PLOT ---

cc_cwb_msea_v1 <- ggplot(
  data = vol_diff_explicit_merged, 
  aes(x = CWB_MSEA, y = perc_diff_VO) 
) +
  geom_point(aes(color = Scenario, shape = Simulation_site), size = 3, alpha = 1) + 
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, color = "black") +
  scale_x_continuous(breaks = seq(-90, 20, by = 20),
                     limits = c(-90, 20))+
  
  # Add your custom, pre-calculated text
  annotate(
    geom = "text",
    x = -90,      # Position at the far left
    y = Inf,         # Position at the very top
    label = annotation_text_cwb_v, # Use the CWB-specific text
    hjust = 0,       # Left-align the text
    vjust = 1.2,     # Nudge it down from the top edge
    size = 4,
    fontface = "italic"
  ) +
  
  scale_color_manual(values = rcp_colors) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    y = "Growing stock difference (%)",
    # Update the x-axis label for CWB
    x = "Change in climatic water balance (mm)",
    # Update the caption for CWB
    #caption = "Slope expressed as % growing stock change per mm CWB difference"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y =  element_text(face = "bold", size = 12),
    axis.title.x =  element_text(face = "bold", size = 12),
    # --- CHANGE 4: Make the legend visible ---
    legend.position = c(0.98, 0.37),
    #legend.position = "none",
    
    # Justify the legend so its top-right corner is at the specified position
    legend.justification = c("right"),
    panel.grid.major = element_blank(),
    legend.background = element_blank(),
    panel.grid.minor = element_blank()
  )
# Display the plot
print(cc_cwb_msea_v1)

library(patchwork)
library(grid)

# Standardize legend titles in all four plots
cc_prec_msea_v1 <- cc_prec_msea_v1 +
  labs(color = "Scenario", shape = "Simulation site")

cc_prec_msea_m1 <- cc_prec_msea_m1 +
  labs(color = "Scenario", shape = "Simulation site")

cc_cwb_msea_v1 <- cc_cwb_msea_v1 +
  labs(color = "Scenario", shape = "Simulation site")

cc_cwb_msea_m1 <- cc_cwb_msea_m1 +
  labs(color = "Scenario", shape = "Simulation site")

# Combine plots and collect duplicate legends
cc1 <- (
  (cc_prec_msea_v1 + cc_prec_msea_m1) /
    (cc_cwb_msea_v1 + cc_cwb_msea_m1)
) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(
      plot.tag.position = "topleft",
      plot.tag = element_text(
        size = 14,
        face = "bold"
      )
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

cc1
ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Figures/CC_Effect.png", plot = cc1, width = 9, height = 10, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  PART 3: SAVE THE NEW PLOTS                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Update the file paths to save the new CWB plots
ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Figures/CWB_Mor_4.pdf", plot = cc_cwb_msea_m1, width = 8, height = 8, dpi = 300)
ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Figures/CWB_GS_4.pdf", plot = cc_cwb_msea_v1, width = 8, height = 8, dpi = 300)



vol_diff_explicit_merged
############ALso variance decomposition
vd <- vol_diff_explicit_merged %>%
  select(Year, Climate_model, Scenario, Simulation_site, perc_diff_cum_VO_tr, perc_diff_VO)
head(vd)





# =============================================================================
# Variance decomposition of projected change: RCM (climate model) vs RCP (scenario)
# Reviewer request on Figures 6-7
#
# Input : vd  (Year, Climate_model, Scenario, Simulation_site, <responses>)
# Design: per site, 5 RCMs (CC1-CC5) x 3 RCPs, one value per cell
#
# Site-level outputs (OUT_DIR):
#   vd_site_decomposition.csv    two-way ANOVA per site: SS share, omega^2, F-test
#   vd_site_signal_noise.csv     paired RCP contrasts across RCMs: SNR, sign agreement
#   vd_site_loo_rcm.csv          leave-one-RCM-out range of SS shares
#   vd_manuscript_table.csv      compact site table for the supplement
#   vd_variance_fractions.png    stacked bars (and fraction-over-time if >1 Year)
# Pooled outputs (all sites together):
#   vd_pooled_decomposition.csv  Site + RCM + RCP + two-way interactions
#   vd_pooled_signal_noise.csv   paired RCP contrasts across RCMs, sites and cells
#   vd_pooled_table.csv          compact pooled table for the supplement
#   vd_pooled_fractions.png      stacked bars of uncertainty shares
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

OUT_DIR   <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/variance_decomposition"
RESPONSES <- c("perc_diff_cum_VO_tr", "perc_diff_VO")
dir.create(OUT_DIR, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 0. Prepare and check the design
# -----------------------------------------------------------------------------
vd <- vd %>%
  mutate(
    Climate_model   = factor(Climate_model),
    Scenario        = factor(Scenario, levels = sort(unique(Scenario))),
    Simulation_site = factor(Simulation_site)
  )

balance <- vd %>% count(Year, Simulation_site, Climate_model, Scenario)
cat("Cells per Year x Site:", paste(unique(table(balance$Year, balance$Simulation_site)), collapse = ","),
    "| obs per cell:", paste(unique(balance$n), collapse = ","), "\n")
# Expected: 15 cells per Year x Site, 1 obs per cell.
# If n > 1 (iLand stochastic replicates kept), the residual becomes a genuine
# internal-variability term and the RCM x RCP interaction is estimated separately.

yr_final <- max(vd$Year)

# -----------------------------------------------------------------------------
# 1. Per-site two-way ANOVA decomposition
#    y = mu + RCM_i + RCP_j + e_ij
#    Balanced design, so sequential SS are unique and shares sum to 1.
#    omega^2 is reported alongside the raw SS share because raw shares inflate
#    terms with more df (RCM has 4 df, RCP has 2).
# -----------------------------------------------------------------------------
decomp_cell <- function(d, y) {
  d <- droplevels(d[!is.na(d[[y]]), ])
  n_rep <- max(table(d$Climate_model, d$Scenario))
  terms <- c("Climate_model", "Scenario")
  if (n_rep > 1) terms <- c(terms, "Climate_model:Scenario")

  a      <- anova(lm(reformulate(terms, y), data = d))
  ms_res <- a["Residuals", "Mean Sq"]
  ss_tot <- sum(a[["Sum Sq"]])

  out <- tibble(
    term    = trimws(rownames(a)),
    df      = a$Df,
    SS      = a[["Sum Sq"]],
    F       = a[["F value"]],
    p       = a[["Pr(>F)"]],
    frac_SS = SS / ss_tot,
    omega2  = ifelse(term == "Residuals", NA_real_,
                     pmax(0, (SS - df * ms_res) / (ss_tot + ms_res)))
  )
  out$omega2[out$term == "Residuals"] <- 1 - sum(out$omega2, na.rm = TRUE)

  res_label <- if (n_rep > 1) "Internal (iLand replicates)" else "RCM x RCP + residual"
  out %>% mutate(term = recode(term,
                               "Climate_model"          = "RCM (climate model)",
                               "Scenario"               = "RCP (scenario)",
                               "Climate_model:Scenario" = "RCM x RCP",
                               "Residuals"              = res_label))
}

# -----------------------------------------------------------------------------
# 2. Does the RCP signal exceed RCM spread?  Paired contrasts across RCMs
#    Each RCM is run under every RCP, so contrasts are paired within RCM.
#    snr        = mean(delta) / sd(delta) across the 5 RCMs
#    sign_agree = number of RCMs agreeing on the sign of the RCP effect
#    ratio_rcp_range_to_rcm_sd = range of RCP means / mean within-RCP RCM sd
# -----------------------------------------------------------------------------
signal_noise <- function(d, y) {
  d <- d[!is.na(d[[y]]), ]
  cell <- d %>%
    group_by(Climate_model, Scenario) %>%
    summarise(v = mean(.data[[y]]), .groups = "drop") %>%
    droplevels()
  rcp_stats <- cell %>%
    group_by(Scenario) %>%
    summarise(m = mean(v), s = sd(v), .groups = "drop")
  wide  <- pivot_wider(cell, names_from = Scenario, values_from = v)
  pairs <- combn(levels(cell$Scenario), 2, simplify = FALSE)

  map_dfr(pairs, function(pr) {
    delta <- wide[[pr[2]]] - wide[[pr[1]]]
    tibble(
      contrast   = paste(pr[2], "minus", pr[1]),
      mean_delta = mean(delta),
      sd_delta   = sd(delta),
      snr        = mean(delta) / sd(delta),
      sign_agree = max(sum(delta > 0), sum(delta < 0)),
      n_rcm      = length(delta),
      p_paired_t = tryCatch(t.test(delta)$p.value, error = function(e) NA_real_)
    )
  }) %>%
    mutate(
      rcp_range                 = diff(range(rcp_stats$m)),
      rcm_sd_within             = mean(rcp_stats$s),
      ratio_rcp_range_to_rcm_sd = rcp_range / rcm_sd_within,
      p_friedman = tryCatch(friedman.test(v ~ Scenario | Climate_model, data = cell)$p.value,
                            error = function(e) NA_real_)
    )
}

# -----------------------------------------------------------------------------
# 3. Leave-one-RCM-out sensitivity (5 RCMs: is one model driving the result?)
# -----------------------------------------------------------------------------
loo_rcm <- function(d, y) {
  mods <- levels(droplevels(d$Climate_model))
  map_dfr(mods, function(m) decomp_cell(d[d$Climate_model != m, ], y) %>% mutate(dropped = m)) %>%
    group_by(term) %>%
    summarise(frac_SS_min = min(frac_SS), frac_SS_max = max(frac_SS),
              rcm_max_effect = dropped[which.max(abs(frac_SS - median(frac_SS)))],
              .groups = "drop")
}

# -----------------------------------------------------------------------------
# 4. Run site-level analysis
# -----------------------------------------------------------------------------
by_site <- function(fun, y) {
  vd %>% group_by(Year, Simulation_site) %>%
    group_modify(~ fun(.x, y)) %>% ungroup() %>% mutate(response = y)
}

site_dec <- map_dfr(RESPONSES, ~ by_site(decomp_cell, .x))
site_sn  <- map_dfr(RESPONSES, ~ by_site(signal_noise, .x))
site_loo <- map_dfr(RESPONSES, ~ by_site(loo_rcm, .x))

write.csv(site_dec, file.path(OUT_DIR, "vd_site_decomposition.csv"), row.names = FALSE)
write.csv(site_sn,  file.path(OUT_DIR, "vd_site_signal_noise.csv"),  row.names = FALSE)
write.csv(site_loo, file.path(OUT_DIR, "vd_site_loo_rcm.csv"),       row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Compact site table (final Year, RCP 8.5 vs RCP 2.6 contrast)
# -----------------------------------------------------------------------------
tab <- site_dec %>%
  filter(Year == yr_final) %>%
  mutate(term_short = case_when(grepl("^RCP", term) ~ "RCP",
                                grepl("^RCM \\(", term) ~ "RCM",
                                TRUE ~ "Resid")) %>%
  select(response, Simulation_site, term_short, frac_SS, omega2, p) %>%
  pivot_wider(names_from = term_short, values_from = c(frac_SS, omega2, p)) %>%
  left_join(
    site_sn %>%
      filter(Year == yr_final, grepl("8.5 minus .*2.6", contrast)) %>%
      select(response, Simulation_site, mean_delta, snr, sign_agree, n_rcm,
             ratio_rcp_range_to_rcm_sd, p_friedman),
    by = c("response", "Simulation_site")
  ) %>%
  select(-any_of("p_Resid")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write.csv(tab, file.path(OUT_DIR, "vd_manuscript_table.csv"), row.names = FALSE)
print(tab, n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# 6. Site-level figures
# -----------------------------------------------------------------------------
term_levels <- c("RCP (scenario)", "RCM (climate model)", "RCM x RCP",
                 "RCM x RCP + residual", "Internal (iLand replicates)")
term_cols <- c("RCP (scenario)" = "#C0504D", "RCM (climate model)" = "#4F81BD",
               "RCM x RCP" = "#9BBB59", "RCM x RCP + residual" = "#BFBFBF",
               "Internal (iLand replicates)" = "#7F7F7F")

p_bar <- site_dec %>%
  filter(Year == yr_final) %>%
  mutate(term = factor(term, levels = rev(term_levels))) %>%
  ggplot(aes(Simulation_site, frac_SS, fill = term)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  facet_wrap(~ response) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = term_cols, breaks = term_levels, drop = TRUE) +
  labs(x = NULL, y = "Share of ensemble variance (sum of squares)", fill = NULL,
       title = paste0("Variance decomposition of projected change, ", yr_final)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

if (n_distinct(vd$Year) > 1) {
  p_time <- site_dec %>%
    mutate(term = factor(term, levels = rev(term_levels))) %>%
    ggplot(aes(Year, frac_SS, fill = term)) +
    geom_area(colour = "white", linewidth = 0.2) +
    facet_grid(response ~ Simulation_site) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    scale_fill_manual(values = term_cols, breaks = term_levels, drop = TRUE) +
    labs(y = "Share of ensemble variance", fill = NULL) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_bar <- patchwork::wrap_plots(p_bar, p_time, ncol = 1, heights = c(1, 1.3))
  }
}

ggsave(file.path(OUT_DIR, "vd_variance_fractions.png"), p_bar,
       width = 9, height = if (n_distinct(vd$Year) > 1) 10 else 4.5, dpi = 300)

# =============================================================================
# 7. POOLED ANALYSIS: all sites together
#
#    y = mu + Site + RCM + RCP + Site:RCM + Site:RCP + RCM:RCP + e
#    One value per Site x RCM x RCP cell, so the residual is the three-way
#    interaction (RCM x RCP x Site). Balanced design, so shares are unique.
#
#    Site main effect = spatial difference in level, not uncertainty.
#    frac_uncert therefore expresses each term relative to the uncertainty
#    terms only (total SS minus Site SS).
#
#    Interpretation of the terms:
#      RCP        scenario signal common to all sites
#      RCM        climate model offset common to all sites
#      Site x RCP how much the scenario signal differs between sites
#      Site x RCM how much the model offset differs between sites
#      RCM x RCP  models disagree on the scenario signal, consistently across sites
#      Residual   model-specific scenario response that also varies by site
# =============================================================================

pooled_labels <- c(
  "Simulation_site"               = "Site (spatial)",
  "Climate_model"                 = "RCM",
  "Scenario"                      = "RCP",
  "Simulation_site:Climate_model" = "Site x RCM",
  "Simulation_site:Scenario"      = "Site x RCP",
  "Climate_model:Scenario"        = "RCM x RCP",
  "Residuals"                     = "RCM x RCP x Site (residual)"
)

decomp_pooled <- function(d, y) {
  d <- droplevels(d[!is.na(d[[y]]), ])
  f <- reformulate(c("Simulation_site", "Climate_model", "Scenario",
                     "Simulation_site:Climate_model", "Simulation_site:Scenario",
                     "Climate_model:Scenario"), y)
  a      <- anova(lm(f, data = d))
  ms_res <- a["Residuals", "Mean Sq"]
  ss_tot <- sum(a[["Sum Sq"]])
  ss_unc <- ss_tot - a["Simulation_site", "Sum Sq"]

  tibble(term = trimws(rownames(a)), df = a$Df, SS = a[["Sum Sq"]],
         F = a[["F value"]], p = a[["Pr(>F)"]]) %>%
    mutate(
      frac_total  = SS / ss_tot,
      frac_uncert = ifelse(term == "Simulation_site", NA_real_, SS / ss_unc),
      omega2      = ifelse(term == "Residuals", NA_real_,
                           pmax(0, (SS - df * ms_res) / (ss_tot + ms_res))),
      term        = recode(term, !!!pooled_labels)
    )
}

# Paired RCP contrasts pooled over sites.
#   *_cell : each Site x RCM pair is one paired contrast (e.g. 6 x 5 = 30)
#   *_rcm  : delta averaged over sites, one value per RCM (n = 5)
#   p_friedman_blocks : Friedman test with Site x RCM as blocks
pooled_signal <- function(d, y) {
  cell <- d %>%
    filter(!is.na(.data[[y]])) %>%
    group_by(Simulation_site, Climate_model, Scenario) %>%
    summarise(v = mean(.data[[y]]), .groups = "drop") %>%
    droplevels()
  wide  <- pivot_wider(cell, names_from = Scenario, values_from = v)
  pairs <- combn(levels(cell$Scenario), 2, simplify = FALSE)

  p_fr <- tryCatch(
    friedman.test(v ~ Scenario | block,
                  data = mutate(cell, block = interaction(Simulation_site, Climate_model,
                                                          drop = TRUE)))$p.value,
    error = function(e) NA_real_)

  map_dfr(pairs, function(pr) {
    delta     <- wide[[pr[2]]] - wide[[pr[1]]]
    rcm_mean  <- tapply(delta, wide$Climate_model, mean)
    tibble(
      contrast         = paste(pr[2], "minus", pr[1]),
      mean_delta       = mean(delta),
      sign_agree_cell  = max(sum(delta > 0), sum(delta < 0)),
      n_cell           = length(delta),
      snr_rcm          = mean(rcm_mean) / sd(rcm_mean),
      sign_agree_rcm   = max(sum(rcm_mean > 0), sum(rcm_mean < 0)),
      n_rcm            = length(rcm_mean),
      p_paired_t_rcm   = tryCatch(t.test(rcm_mean)$p.value, error = function(e) NA_real_),
      p_friedman_blocks = p_fr
    )
  })
}

run_pooled <- function(fun) {
  map_dfr(RESPONSES, function(y) {
    vd %>% group_by(Year) %>% group_modify(~ fun(.x, y)) %>%
      ungroup() %>% mutate(response = y)
  })
}

pooled_dec <- run_pooled(decomp_pooled)
pooled_sn  <- run_pooled(pooled_signal)

write.csv(pooled_dec, file.path(OUT_DIR, "vd_pooled_decomposition.csv"), row.names = FALSE)
write.csv(pooled_sn,  file.path(OUT_DIR, "vd_pooled_signal_noise.csv"),  row.names = FALSE)

# Compact pooled table (final Year, RCP 8.5 vs RCP 2.6)
pooled_tab <- pooled_dec %>%
  filter(Year == yr_final) %>%
  select(response, term, frac_uncert, p) %>%
  mutate(term = factor(term, levels = pooled_labels)) %>%
  arrange(response, term) %>%
  left_join(
    pooled_sn %>%
      filter(Year == yr_final, grepl("8.5 minus .*2.6", contrast)) %>%
      select(response, mean_delta, sign_agree_cell, n_cell,
             snr_rcm, sign_agree_rcm, n_rcm, p_friedman_blocks),
    by = "response"
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write.csv(pooled_tab, file.path(OUT_DIR, "vd_pooled_table.csv"), row.names = FALSE)
print(pooled_tab, n = Inf, width = Inf)

# Pooled figure: shares of uncertainty (Site main effect excluded)
pooled_cols <- c("RCP" = "#C0504D", "RCM" = "#4F81BD",
                 "Site x RCP" = "#E6A0A0", "Site x RCM" = "#A9C4E8",
                 "RCM x RCP" = "#9BBB59", "RCM x RCP x Site (residual)" = "#BFBFBF")

p_pooled <- pooled_dec %>%
  filter(Year == yr_final, term != "Site (spatial)") %>%
  mutate(term = factor(term, levels = rev(names(pooled_cols)))) %>%
  ggplot(aes(response, frac_uncert, fill = term)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.2) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = pooled_cols, breaks = names(pooled_cols)) +
  labs(x = NULL, y = "Share of ensemble uncertainty (site main effect excluded)", fill = NULL,
       title = paste0("Pooled variance decomposition, ", yr_final)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(nrow = 2))

ggsave(file.path(OUT_DIR, "vd_pooled_fractions.png"), p_pooled,
       width = 9, height = 3.8, dpi = 300)









# =============================================================================
# Variance decomposition, following the Baldo et al. (RESONATE WP4) approach
#   a) Discriminatory power = SD of RCP x RCM cell means
#   b) Two-way ANOVA with replication, % of total SS per source
#
#   Value ~ Scenario * Climate_model, all sites together
#   Replicates within each RCP x RCM cell = the 6 simulation sites
#   -> "Residual" = between-site variation (site differences, which are not
#      uncertainty, plus site-specific responses to RCM and RCP)
#
# Input: vd (Year, Climate_model, Scenario, Simulation_site, responses)
# =============================================================================

# library(tidyverse)
# library(patchwork)
# 
# OUT_DIR <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/variance_decomposition"
# dir.create(OUT_DIR, showWarnings = FALSE)
# 
# var_labels <- c(
#   perc_diff_cum_VO_tr = "Cumulative mortality (% change)",
#   perc_diff_VO        = "Growing stock (% change)"
# )
# 
# # ── 1. Long format, final year, z-score per variable ─────────────────────────
# df_long <- vd |>
#   filter(Year == max(Year)) |>
#   select(Simulation_site, Climate_model, Scenario, all_of(names(var_labels))) |>
#   pivot_longer(cols = all_of(names(var_labels)),
#                names_to = "Variable", values_to = "Value") |>
#   mutate(
#     Variable = recode(Variable, !!!var_labels),
#     Scenario = as.factor(Scenario),
#     RCM      = as.factor(Climate_model)
#   ) |>
#   group_by(Variable) |>
#   mutate(Value = scale(Value)[, 1]) |>   # z-score: makes disc power comparable across variables
#   ungroup()                              # (ANOVA % shares and p-values are unaffected)
# 
# # Check: 3 RCP x 5 RCM = 15 cells, 6 sites per cell (balanced -> Type I SS are unique)
# print(count(df_long, Variable, Scenario, RCM) |> count(Variable, n, name = "cells"))
# 
# # ── 2. Discriminatory power: SD of RCP x RCM cell means ──────────────────────
# disc_power <- df_long |>
#   group_by(Variable, Scenario, RCM) |>
#   summarise(cell_mean = mean(Value), .groups = "drop") |>
#   group_by(Variable) |>
#   summarise(disc = sd(cell_mean), .groups = "drop")
# 
# # ── 3. Two-way ANOVA with replication ────────────────────────────────────────
# anova_res <- df_long |>
#   group_by(Variable) |>
#   group_modify(~ {
#     fit <- aov(Value ~ Scenario * RCM, data = .x)
#     s   <- summary(fit)[[1]]
#     rownames(s) <- trimws(rownames(s))
#     ss_total <- sum(s$`Sum Sq`)
#     tibble(
#       pct_rcp   = s["Scenario",     "Sum Sq"] / ss_total * 100,
#       pct_rcm   = s["RCM",          "Sum Sq"] / ss_total * 100,
#       pct_inter = s["Scenario:RCM", "Sum Sq"] / ss_total * 100,
#       pct_resid = s["Residuals",    "Sum Sq"] / ss_total * 100,
#       F_rcp     = s["Scenario",     "F value"],
#       F_rcm     = s["RCM",          "F value"],
#       F_inter   = s["Scenario:RCM", "F value"],
#       p_rcp     = s["Scenario",     "Pr(>F)"],
#       p_rcm     = s["RCM",          "Pr(>F)"],
#       p_inter   = s["Scenario:RCM", "Pr(>F)"]
#     )
#   }) |>
#   ungroup()
# 
# # Optional residual diagnostics (run manually)
# # test <- df_long |> filter(Variable == "Growing stock (% change)")
# # fit  <- aov(Value ~ Scenario * RCM, data = test)
# # plot(fit); shapiro.test(residuals(fit))
# # bartlett.test(Value ~ interaction(Scenario, RCM), data = test)
# 
# # ── 4. Combine & annotate significance (stars = RCP effect) ──────────────────
# stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
# 
# results <- disc_power |>
#   left_join(anova_res, by = "Variable") |>
#   mutate(Variable = fct_reorder(Variable, disc),
#          sig      = stars(p_rcp))
# 
# print(results, width = Inf)
# write.csv(results, file.path(OUT_DIR, "vd_baldo_style_results.csv"), row.names = FALSE)
# 
# max_sd <- max(results$disc)
# 
# # ── 5. Colours ───────────────────────────────────────────────────────────────
# var_colors <- c(
#   "Cumulative mortality (% change)" = "#D85A30",
#   "Growing stock (% change)"        = "#185FA5"
# )
# 
# src_colors <- c(
#   "RCP (scenario)"            = "#C0504D",
#   "RCM (climate model)"       = "#4F81BD",
#   "Interaction"               = "#D0CECA",
#   "Residual (between sites)"  = "#787670"
# )
# 
# # ── 6. Panel A: discriminatory power ─────────────────────────────────────────
# pA <- ggplot(results, aes(x = disc, y = Variable, fill = Variable)) +
#   geom_col(width = 0.6) +
#   geom_text(aes(label = sig), hjust = -0.3, size = 3.5) +
#   scale_fill_manual(values = var_colors) +
#   scale_x_continuous(expand = c(0, 0), limits = c(0, max_sd * 1.15)) +
#   labs(x = "Discriminatory power (SD of cell means, z-scored)",
#        y = NULL, title = "a") +
#   theme_minimal(base_size = 11) +
#   theme(panel.grid.major.y = element_blank(),
#         panel.grid.minor   = element_blank(),
#         legend.position    = "none")
# 
# # ── 7. Panel B: variance decomposition ───────────────────────────────────────
# plot_decomp <- results |>
#   select(Variable, pct_rcp, pct_rcm, pct_inter, pct_resid) |>
#   pivot_longer(cols = starts_with("pct_"), names_to = "source", values_to = "pct") |>
#   mutate(
#     source = recode(source,
#                     pct_rcp   = "RCP (scenario)",
#                     pct_rcm   = "RCM (climate model)",
#                     pct_inter = "Interaction",
#                     pct_resid = "Residual (between sites)"),
#     source = factor(source, levels = rev(names(src_colors)))
#   )
# 
# pB <- ggplot(plot_decomp, aes(x = pct, y = Variable, fill = source)) +
#   geom_col(position = "stack") +
#   scale_fill_manual(values = src_colors, breaks = names(src_colors),
#                     name = "Variance source") +
#   scale_x_continuous(expand = c(0, 0), limits = c(0, 101),
#                      labels = function(x) paste0(x, "%")) +
#   labs(x = "Variance decomposition (%)", y = NULL, title = "b") +
#   theme_minimal(base_size = 11) +
#   theme(panel.grid.major.y    = element_blank(),
#         panel.grid.minor      = element_blank(),
#         axis.text.y           = element_blank(),
#         legend.position       = "bottom",
#         legend.title          = element_text(face = "bold"),
#         legend.box.background = element_rect(colour = "black"),
#         legend.margin         = margin(4, 4, 4, 4)) +
#   guides(fill = guide_legend(nrow = 2))
# 
# # ── 8. Combine & save ────────────────────────────────────────────────────────
# final_panel <- pA + pB
# final_panel
# 
# ggsave(file.path(OUT_DIR, "Fig_VarDec_baldo_style.pdf"), final_panel,
#        width = 28, height = 8, units = "cm")