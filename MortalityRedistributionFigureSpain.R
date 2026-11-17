# Packages
library(readxl)
library(ggplot2)

# Path to the file (change to where it sits on your machine)
f <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Climate/ObsMor_anomaly.xlsx"

# Read the first sheet; skip the title row so row 2 becomes the header
mort <- read_excel(f, sheet = 1)

# Check the data
str(mort)
head(mort)
# Make sure the types are right
mort$Year <- as.integer(mort$Year)
mort$ObsMor_F <- as.numeric(mort$ObsMor_F)
mort$Site <- factor(mort$Site)
mort <- mort %>% filter(!Country == "Switzerland")
unique(mort$Country)

# Time series plot, one coloured line per site
p <- ggplot(mort, aes(x = Year, y = ObsMor_F, colour = Site)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = sort(unique(mort$Year))) +
  labs(
    x = "Year",
    y = "Redistributed Annual mortality (%)",
    colour = "Site",
    #caption = "Redistributed annual mortality (%)"
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),   # x, y in 0–1 panel units
    legend.justification = c(1, 1),           # anchor the legend's top-right corner there
    #legend.background = element_rect(fill = alpha("white", 0.8), colour = "grey70"),
    #legend.key = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
print(p)

# Optional: save
ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Figures/Redistributed_annual_mortality.png", width = 8, height = 5, dpi = 300)
