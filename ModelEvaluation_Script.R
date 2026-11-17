##Model evaluation scripts

library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(ggstats)
library(glmmTMB)
library(mgcv)
library(DHARMa)
library(GGally)
library(performance)
library(glmmTMB)
library(lmerTest)
library(lme4)
library(ggpmisc)
library(ggplot2)
library(ggExtra)

##Read the data

data <- read_xlsx("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/IGA_v2/ExcelSheets/MeanTreeVolumes.xlsx")
# Convert Country and Site to factors
data <- data %>%
  mutate(Country = as.factor(Country),
         Site = as.factor(Site))

head(data)
data <- data %>%
  mutate(Site = case_when(
    Site == "Visp" ~ "CH_Visp",
    Site == "Simsite0" ~ "SP_Prades",
    Site == "Simsite1" ~ "SP_Arcalis1",
    Site == "Simsite2" ~ "SP_Arcalis2",
    Site == "Simsite3" ~ "SP_Arcalis3",
    Site == "Simsite4" ~ "SP_Arcalis4",
    TRUE ~ Site
  ))
head(data)
filtered_data <- data %>%
  filter(
    (Country == "Spain" & Year %in% c(2007, 2015)) |
      (Country == "Switzerland" & Year %in% c(1999, 2005, 2009, 2014, 2017))
  )


##1: First we will visualize the regression of the mean tree volume with the simulated mean tree volume
mtv <- ggplot(data, aes(x = ObsMTV, y = SimMTV, color = Site)) + 
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm",formula = y ~ x,se = TRUE,aes(fill = Site),
              alpha = 0.2, level = 0.95)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")+
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, parse = TRUE, size = 5, label.x = "left", label.y = "top"
  ) +
  #scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(
    #title = paste("Regression of Simulated Mean Tree Volume against Observed"),
    x = "Observed",
    y = "Simulated",
    color = "Site"
  ) +
  theme_bw() +
  theme(
    # Axis label styling
    axis.title = element_text(size = 14, face = "bold"), # Increases size for both X and Y titles
    
    # Legend positioning
    legend.position = c(0.95, 0.30), # Position legend: 95% to the right, 5% from the bottom
    legend.justification = c("right", "bottom"), # Anchor the legend's bottom-right corner
    legend.background = element_rect(fill = "transparent"), # Make legend background transparent
    legend.box.background = element_rect(fill = "transparent", color = NA), # Remove box background
    
    # Panel styling (as you had before)
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
print(mtv)

mtv1 <- ggplot(filtered_data, aes(x = ObsMTV, y = SimMTV, color = Country)) + 
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm",formula = y ~ x,se = TRUE,aes(fill = Country),
              alpha = 0.2, level = 0.95)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")+
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, parse = TRUE, size = 5, label.x = "left", label.y = "top"
  ) +
  #scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(
    #title = paste("Regression of Simulated Mean Tree Volume against Observed"),
    x = "Observed",
    y = "Simulated",
    color = "Country"
  ) +
  theme_bw() +
  theme(
    # Axis label styling
    axis.title = element_text(size = 14, face = "bold"), # Increases size for both X and Y titles
    
    # Legend positioning
    legend.position = c(0.95, 0.30), # Position legend: 95% to the right, 5% from the bottom
    legend.justification = c("right", "bottom"), # Anchor the legend's bottom-right corner
    legend.background = element_rect(fill = "transparent"), # Make legend background transparent
    legend.box.background = element_rect(fill = "transparent", color = NA), # Remove box background
    
    # Panel styling (as you had before)
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
print(mtv1)
df_increments <- filtered_data %>%
  arrange(Site, Year) %>%
  group_by(Site) %>%
  mutate(
    ObsIncrement = ObsMTV - lag(ObsMTV),
    SimIncrement = SimMTV - lag(SimMTV),
    YearDiff = Year - lag(Year),
    ObsAnnualIncrement = ObsIncrement / YearDiff,
    SimAnnualIncrement = SimIncrement / YearDiff
  ) %>%
  ungroup()
print(df_increments)

df_plot <- df_increments %>%
  filter(!is.na(ObsAnnualIncrement) & !is.na(SimAnnualIncrement))
print(df_plot)


# Filter your data for convenience
spain_data <- subset(df_plot, Country == "Spain")
swiss_data <- subset(df_plot, Country == "Switzerland")


mti <- ggplot(
  df_plot,
  aes(
    x = ObsAnnualIncrement,
    y = SimAnnualIncrement,
    color = Country
  )
) +
  
  # Scatter points
  geom_point(
    aes(shape = Site),
    size = 2,
    alpha = 0.8
  ) +
  
  # Linear regression lines with 95% confidence intervals
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    alpha = 0.2,
    level = 0.95
  ) +
  
  # 1:1 reference line
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "black",
    linetype = "dashed"
  ) +
  
  # R-squared label for Switzerland
  stat_poly_eq(
    data = swiss_data,
    aes(
      x = ObsAnnualIncrement,
      y = SimAnnualIncrement,
      color = Country,
      label = after_stat(rr.label)
    ),
    formula = y ~ x,
    parse = TRUE,
    size = 4,
    label.x = "right",
    label.y = 0.98,
    inherit.aes = FALSE
  ) +
  
  # R-squared label for Spain
  stat_poly_eq(
    data = spain_data,
    aes(
      x = ObsAnnualIncrement,
      y = SimAnnualIncrement,
      color = Country,
      label = after_stat(rr.label)
    ),
    formula = y ~ x,
    parse = TRUE,
    size = 4,
    label.x = "right",
    label.y = 0.93,
    inherit.aes = FALSE
  ) +
  
  # Axis labels
  labs(
    x = expression(bold("Observed annual increment (" * m^3 *
                          " tree"^{-1} * " yr"^{-1} * ")")),
    y = expression(bold("Simulated annual increment (" * m^3 *
                          " tree"^{-1} * " yr"^{-1} * ")")),
    color = "Country",
    shape = "Site"
  ) +
  
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "none",
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(mti)

####Mortality######
data2 <- read_excel("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/ObsMor_anomaly.xlsx")
head(data2)
data2 <- data2 %>%
  mutate(Site = case_when(
    Site == "Visp" ~ "CH_Visp",
    Site == "Prades_Simsite0" ~ "SP_Prades",
    Site == "Arcalis_Simsite1" ~ "SP_Arcalis1",
    Site == "Arcalis_Simsite2" ~ "SP_Arcalis2",
    Site == "Arcalis_Simsite3" ~ "SP_Arcalis3",
    Site == "Arcalis_Simsite4" ~ "SP_Arcalis4",
    TRUE ~ Site
  ))

spain_data1 <- subset(data2, Country == "Spain")
swiss_data1 <- subset(data2, Country == "Switzerland")



dat <- data2 %>%
  filter(Site == "CH_Visp")

cor.test(dat$ObsMor_F, dat$SimMor, method = "spearman", exact = FALSE)

# stats_lab <- data2 %>%
#   group_by(Country) %>%
#   summarise(
#     rho = cor(ObsMor_F, SimMor, method = "spearman", use = "complete.obs"),
#     r2  = summary(lm(SimMor ~ ObsMor_F))$r.squared,
#     #p   = summary(lm(SimMor ~ ObsMor_F))$coefficients[2, 4],
#     .groups = "drop"
#   ) %>%
#   #mutate(
#     #p_lab = case_when(
#       #p < 0.001 ~ "italic(p) < 0.001",
#       #p < 0.01  ~ "italic(p) < 0.01",
#       #p < 0.05  ~ "italic(p) < 0.05",
#      # TRUE      ~ sprintf("italic(p) == \"%.3f\"", p)
#     #),
#     label = sprintf("rho == \"%.2f\"~~~italic(R)^2 == \"%.2f\"~~~%s", rho, r2)#, p_lab),
#     npcx  = 0.97,
#     npcy  = ifelse(Country == "Switzerland", 0.98, 0.93)
#   )

stats_lab <- data2 %>%
  group_by(Country) %>%
  summarise(
    rho = cor(ObsMor_F, SimMor, method = "spearman", use = "complete.obs"),
    r2  = summary(lm(SimMor ~ ObsMor_F))$r.squared,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("rho == \"%.2f\"~~~italic(R)^2 == \"%.2f\"", rho, r2),
    npcx  = 0.97,
    npcy  = ifelse(Country == "Switzerland", 0.98, 0.93)
  )

mor <- ggplot(data2, aes(x = ObsMor_F, y = SimMor, color = Country)) + 
  
  # The shape aesthetic is now LOCAL to geom_point
  geom_point(aes(shape = Site), size = 2, alpha = 0.8) + 
  
  # This geom_smooth now only inherits 'color = Country' for grouping, which is correct
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              alpha = 0.2, level = 0.95) +
  
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  
  geom_text_npc(
    data = stats_lab,
    aes(npcx = npcx, npcy = npcy, label = label, color = Country),
    parse = TRUE, hjust = 1, size = 4,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  
  labs(
    x = "Observed annual mortality (% of trees per year)",
    y = "Simulated annual mortality (% of trees per year) ",
    color = "Country",
    shape = "Site" # This correctly labels the shape legend
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = c(0.95, 0.25),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(mor)

library(patchwork)
eval <- mti + mor
print(eval)
ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Figures/ModelEval.png", plot = eval, width = 10, height = 6, dpi = 300)



# Simple and robust: Use mean absolute bias
absolute_bias_results <- data2 %>%
  group_by(Country) %>%
  summarize(
    # For each row, calculate (Sim - Obs), then find the mean
    mean_absolute_bias = mean(SimMor - ObsMor)
  )

# Print the results
print(absolute_bias_results)
