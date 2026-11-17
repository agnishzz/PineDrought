#Trend analysis of mortality dynamics
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

spain<- readxl::read_excel("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/WaterDynamics_Spain.xlsx")
swiss<- readxl::read_excel("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/WaterDynamics_Switzerland.xlsx")



head(spain)
head(swiss)
tail(spain)
tail(swiss)

water <- bind_rows(spain, swiss)
head(water)
tail(water)
writexl::write_xlsx(water, "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/WaterDynamics_Both.xlsx")


plot_data1 <- water %>%
  group_by(Simulation_site, Scenario, Year) %>%
  summarise(
    median_water_mm = median(mean_swc_gs_mm, na.rm = TRUE),
    min_water_mm = min(mean_swc_gs_mm, na.rm = TRUE),
    max_water_mm = max(mean_swc_gs_mm, na.rm = TRUE),
    .groups = 'drop' # Drop grouping after summarising
  )

# View the prepared data
head(plot_data1)

# trend_analysis1 <- plot_data1 %>%
#   # Group by the variables that define each trend line
#   group_by(Simulation_site, Scenario) %>%
#   # Perform a linear model (lm) and tidy the output
#   do(tidy(lm(median_gwl_m3 ~ Year, data = .))) %>%
#   # We only need the slope of the 'Year' term
#   filter(term == "Year") %>%
#   # Select relevant columns and format the p-value
#   select(Simulation_site, Scenario, estimate, p.value) %>%
#   rename(slope = estimate) %>%
#   mutate(
#     significance = case_when(
#       p.value < 0.001 ~ "***",
#       p.value < 0.01  ~ "**",
#       p.value < 0.05  ~ "*",
#       TRUE            ~ "ns" # not significant
#     ),
#     # Create the text label for the plot
#     trend_label = paste0(
#       "Trend: ", round(slope, 2), " ", significance
#     )
#   )
# 
# # To position the labels on the plot, we'll join this back
# # with a summary of the y-axis range for each facet.
# # We'll place the labels at the top of each panel.
# label_positions1 <- plot_data1 %>%
#   group_by(Simulation_site) %>%
#   summarise(y_max = max(max_gwl_m3, na.rm = TRUE)) %>%
#   # Join with the trend labels
#   left_join(trend_analysis, by = "Simulation_site") %>%
#   # Create distinct y-positions for each label within a facet to avoid overlap
#   group_by(Simulation_site) %>%
#   mutate(y_position = y_max * (1 - 0.1 * (row_number() - 1)))
# 
# 
# # View the final labels data
# print(label_positions1)
unique(plot_data1$Scenario)

my_colors <- c(
  "Refclim" = "black",
  "RCP 2.6" = "blue4", # A nice blue
  "RCP 4.5" = "orange", # A warm amber/orange
  "RCP 8.5" = "brown4"  # A strong red
)

ggplot(plot_data1, aes(x = Year, group = Scenario)) +
  
  # 1. Ribbon showing the min-max range across climate models
  geom_ribbon(aes(ymin = min_water_mm, ymax = max_water_mm, fill = Scenario), alpha = 0.3) +
  
  # 2. Median line
  geom_line(aes(y = median_water_mm, color = Scenario), linewidth = 1.5) +
  
  # 3. Add the trend line from the linear model
  # geom_smooth(aes(y = median_vo, color = Scenario),
  #             method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.8) +
  # 
  # 4. Add the trend analysis text labels
  # geom_text(data = label_positions1,
  #           aes(label = trend_label, color = Scenario, y = y_position),
  #           x = -Inf, 
  #           hjust = -0.1,
  #           vjust = 1.2,
  #           fontface = "bold") +
  # 
  # 5. Facet by simulation site with free y-scales
  facet_wrap(~ Simulation_site, scales = "free") +
  
  # 6. Apply the custom color palette to both fill and color
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  
  # 7. Apply a clean theme and add informative labels
  theme_bw() +
  labs(
    title = "Water dynamics Over Time by Site and Scenario",
    #subtitle = "Ribbon = Min/Max of Climate Models, Line = Median, Dashed = Linear Trend",
    y = "Mean soil water content in growing season (mm)",
    x = "Year",
    color = "Scenario",
    fill = "Scenario"
  ) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"))


head(water)


#############Testing###################
library(ggpubr)
water$Scenario <- factor(water$Scenario, levels = c("Refclim", "RCP 2.6", "RCP 4.5", "RCP 8.5"))
Scenario_colors <- c(
  "Refclim"     = "skyblue3",
  "RCP 2.6"     = "springgreen4",
  "RCP 4.5"     = "gold2",
  "RCP 8.5"     = "maroon"
)

# Define which groups to compare
# A list of character vectors, where each vector is a pair to compare
my_comparisons <- list(
  c("Refclim", "RCP 2.6"),
  c("Refclim", "RCP 4.5"),
  c("Refclim", "RCP 8.5")
)


# Create the plot with adjustments for visibility
p_brackets_fixed <- ggplot(water, aes(x = Scenario, y = mean_swc_gs_mm, fill = Scenario)) +
  geom_boxplot(alpha = 0.8) +
  
  # Add the pairwise comparisons using brackets
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",
    label = "p.signif",
    step.increase = 0.1  # MODIFICATION: Increase spacing between stacked brackets
  ) +
  
  # Facet by site
  facet_wrap(~ Simulation_site, scales = "free_y") +
  
  # *** KEY ADDITION: Add extra padding to the top of the y-axis ***
  # The second value in mult=c(bottom, top) adds 25% padding to the top.
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  
  # Apply colors and labels
  scale_fill_manual(values = Scenario_colors) +
  labs(
    title = "Growing Season Soil Water Content Across Scenarios",
    subtitle = "Faceted by Site with Pairwise T-tests vs. Refclim",
    x = "Climate Scenario",
    y = "Mean Growing Season Soil Water Content (mm)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "bottom"
  )

# Display the final, corrected plot
print(p_brackets_fixed)
ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/IGA_v2/Plots_outputs/SoilWaterContent.tif", plot = p_brackets_fixed, width = 10, height = 10, dpi = 300)



#########################
head(plot_data1)


head(water)
library(tidyverse)
library(patchwork)

# --- Prepare observed-period, Refclim data at landscape-aggregate level ---
water_obs <- water %>%
  filter(
    Scenario == "Refclim",
    ru == -1, rid == -1,                       # keep aggregate rows only
    (Simulation_site == "CH_Visp"  & Year >= 1998 & Year <= 2020) |
      (Simulation_site != "CH_Visp"  & Year >= 2007 & Year <= 2015)
  ) %>%
  mutate(
    Simulation_site = factor(
      Simulation_site,
      levels = c("SP_Arcalis1", "SP_Arcalis2", "SP_Arcalis3", "SP_Arcalis4", "SP_Prades", "CH_Visp")
    )
  )

site_cols <- c(
  "SP_Arcalis1" = "#D95F02", "SP_Arcalis2" = "#E6AB02", "SP_Arcalis3" = "#A6761D",
  "SP_Arcalis4" = "#B2182B", "SP_Prades"    = "#7570B3", "CH_Visp"   = "#1B7837"
)

# --- 1. Time series ---
p_ts <- ggplot(water_obs, aes(Year, mean_swc_gs_mm,
                              colour = Simulation_site, group = Simulation_site)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = site_cols, name = "Site") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 8)) +
  labs(
    x = "Year",
    y = expression("Mean growing-season SWC (mm)"),
    title = "Simulated growing-season soil water content, observed period (Refclim)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 11, face = "plain")
  )

# --- 2. Boxplots ---
p_box <- ggplot(water_obs, aes(Simulation_site, mean_swc_gs_mm, fill = Simulation_site)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5, width = 0.6, alpha = 0.75) +
  geom_jitter(width = 0.12, size = 1.1, alpha = 0.45, colour = "grey25") +
  scale_fill_manual(values = site_cols, guide = "none") +
  labs(
    x = NULL,
    y = expression("Mean growing-season SWC (mm)"),
    title = "Distribution by site, observed period (Refclim)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(size = 11, face = "plain")
  )

p_ts / p_box + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

ggsave("swc_gs_observed.png", width = 9, height = 8, dpi = 400)