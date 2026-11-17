library(here)
library(ggplot2)
library(dplyr)

library(gridExtra)
library(grid)

path0 <- "E:/Outputs_Spain/S0_Obs_ref.sqlite"
path1 <- "E:/Outputs_Spain/S1_Obs_ref.sqlite"
path2 <- "E:/Outputs_Spain/S2_Obs_ref.sqlite"
path3 <- "E:/Outputs_Spain/S3_Obs_ref.sqlite"
path4 <- "E:/Outputs_Spain/S4_Obs_ref.sqlite"
path5 <- "E:/Outputs_Switzerland/Visp_Obs_ref.sqlite"

file0 <- here::here(path0)
file1 <- here::here(path1)
file2 <- here::here(path2)
file3 <- here::here(path3)
file4 <- here::here(path4)
file5 <- here::here(path5)


name0 <- "SP_Prades0"
name1 <- "SP_Arcalis1"
name2 <- "SP_Arcalis2"
name3 <- "SP_Arcalis3"
name4 <- "SP_Arcalis4"
name5 <- "CH_Visp"

# Read in data using the RSQLite package
sqlite.driver <- RSQLite::dbDriver("SQLite")
db0 <-
  RSQLite::dbConnect(sqlite.driver, dbname = file0) # connect to the file
tables.in.the.file0 <-
  RSQLite::dbListTables(db0) # explore the tables in the file1
print(tables.in.the.file0)

# We will work with "tree" table and tree-scale data:
tree0 <- RSQLite::dbReadTable(db0, "tree")
RSQLite::dbDisconnect(db0) # disconnect to the file

# Read in data using the RSQLite package
sqlite.driver <- RSQLite::dbDriver("SQLite")
db1 <-
  RSQLite::dbConnect(sqlite.driver, dbname = file1) # connect to the file
tables.in.the.file1 <-
  RSQLite::dbListTables(db1) # explore the tables in the file1
print(tables.in.the.file1)

# We will work with "tree" table and tree-scale data:
tree1 <- RSQLite::dbReadTable(db1, "tree")
RSQLite::dbDisconnect(db1) # disconnect to the file


# Read in data using the RSQLite package
sqlite.driver <- RSQLite::dbDriver("SQLite")
db2 <-
  RSQLite::dbConnect(sqlite.driver, dbname = file2) # connect to the file
tables.in.the.file2 <-
  RSQLite::dbListTables(db2) # explore the tables in the file1
print(tables.in.the.file2)

# We will work with "tree" table and tree-scale data:
tree2 <- RSQLite::dbReadTable(db2, "tree")
RSQLite::dbDisconnect(db2) # disconnect to the file

# Read in data using the RSQLite package
sqlite.driver <- RSQLite::dbDriver("SQLite")
db3 <-
  RSQLite::dbConnect(sqlite.driver, dbname = file3) # connect to the file
tables.in.the.file3 <-
  RSQLite::dbListTables(db3) # explore the tables in the file1
print(tables.in.the.file3)

# We will work with "tree" table and tree-scale data:
tree3 <- RSQLite::dbReadTable(db3, "tree")
RSQLite::dbDisconnect(db3) # disconnect to the file

# Read in data using the RSQLite package
sqlite.driver <- RSQLite::dbDriver("SQLite")
db4 <-
  RSQLite::dbConnect(sqlite.driver, dbname = file4) # connect to the file
tables.in.the.file4 <-
  RSQLite::dbListTables(db4) # explore the tables in the file1
print(tables.in.the.file4)

# We will work with "tree" table and tree-scale data:
tree4 <- RSQLite::dbReadTable(db4, "tree")
RSQLite::dbDisconnect(db4) # disconnect to the file

# Read in data using the RSQLite package
sqlite.driver <- RSQLite::dbDriver("SQLite")
db5 <-
  RSQLite::dbConnect(sqlite.driver, dbname = file5) # connect to the file
tables.in.the.file5 <-
  RSQLite::dbListTables(db5) # explore the tables in the file1
print(tables.in.the.file5)

# We will work with "tree" table and tree-scale data:
tree5 <- RSQLite::dbReadTable(db5, "tree")
RSQLite::dbDisconnect(db5) # disconnect to the file



tree <- rbind(tree0 |> dplyr::mutate(Site = name0),
              tree1 |> dplyr::mutate(Site = name1),
              tree2 |> dplyr::mutate(Site = name2),
              tree3 |> dplyr::mutate(Site = name3),
              tree4 |> dplyr::mutate(Site = name4),
              tree5 |> dplyr::mutate(Site = name5))
head(tree)

print(tree)

tree <- tree %>%
  mutate(species = case_when(
    species == "acps" ~ "Acer pseudoplatanus",
    species == "frex" ~ "Fraxinus excelsior",
    species == "qupu" ~ "Quercus pubescens",
    species == "bepe" ~ "Betula pendula",
    species == "pisy" ~ "Pinus sylvestris",
    species == "soar" ~ "Sorbus aria",
    TRUE ~ species
  ))

tree_year0 <- tree |> dplyr::filter(year == 0)


head(tree_year0)
species_expr <- setNames(
  lapply(unique(tree_year0$species), function(s) bquote(italic(.(s)))),
  unique(tree_year0$species)
)

g1 <- ggplot(tree_year0, aes(x = x, y = y, color = species)) +
  geom_point(size = tree_year0$dbh / 15, alpha = 1) +
  facet_wrap(~ Site) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 14, face = "bold"),
    legend.text     = element_text(size = 13),
    legend.key.size = unit(1, "cm")          # space each key takes up
  ) +
  labs(
    y = "y (in metres)",
    x = "x (in metres)",
    color = "Species"
  ) +
  coord_equal() +
  scale_color_manual(
    values = scales::hue_pal()(length(unique(tree_year0$species))),
    labels = species_expr
  )

print(g1)

ggsave("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/Figures/TreePositions_3.tif", plot = g1, width = 10, height = 6, dpi = 300)

head(tree5)
visp_prop <- tree5 %>%
  filter(year == 0) %>%
  group_by(species) %>%
  summarise(total_basalArea = sum(basalArea, na.rm = TRUE)) %>%
  mutate(perc_basalArea = 100 * total_basalArea / sum(total_basalArea))
visp_prop


#
tree_sum <- tree_year0 %>% group_by(Site, species) %>% summarise(NT = n(),
                                                                (Mean_dbh = mean(dbh)), (Total_BasalArea = sum (basalArea)), MH = mean(height),
                                                                 TV = sum(volume_m3))
tree_sum

writexl::write_xlsx(tree_sum, "C:/R_etal/IGA_v2/ExcelSheets/tree_initial.xlsx")
