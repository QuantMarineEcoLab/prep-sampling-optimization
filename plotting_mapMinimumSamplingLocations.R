# Create a map of study locations

# Load plyr before dplyr to avoid package conflict issues
library(dplyr)
library(tidyverse)
library(here)

# Load data
data <- read.csv(here("data/Great Bay Environmental Data/clean_imputedCART_Great Bay_environmental data.csv"))

# Source functions needed for this script
source(here("sample_optimization/scripts/minimumTime_functions.R"))

# Clean data ----

# Clean site names
data$sitename <- clean_vec(data$sitename)

# Clean variable_units
data$variable_units <- clean_vec(data$variable_units)

# Filter only to these variables
data.filtered <- data %>%
  dplyr::filter(variable_units == "dissolved_oxygen_mg_l" |
    variable_units == "dissolved_oxygen_saturation_percent" |
    variable_units == "nitrogen_nitrite_no2_nitrate_no3_as_n_diss_mg_l" |
    variable_units == "solids_suspended_total_mg_l" |
    variable_units == "temperature_water_deg_c") %>%
  dplyr::filter(sitename != "rte_108_bridge_mill_pond") %>% # No sampling at this location for water temperature (but have for the 4 other variables
  dplyr::select(site = sitename, variable_units) %>%
  unique()

# Bring in site coordinates
site.loc <- read.csv("./data/metadata/PREP_sampling_locations.csv")

# Merge minimum time with site map data
site.map <- left_join(data.filtered, site.loc)

# Check # of variables sampled for each site
site_variable_counts <- aggregate(variable_units ~ site, site.map, FUN = function(x) length(unique(x)))
colnames(site_variable_counts) <- c("site", "variable_count")

# Set up the map ----
library(sf)

# Extract shapefile for USA
USA <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))

# Download New Hampshire spatial data for rivers, streams and wetlands
# https://www.nhgeodata.unh.edu/search?tags=hydrography
stream <- st_read(here("data/spatial_data/New_Hampshire_Stream_Order_Dataset/New_Hampshire_Stream_Order_Dataset.shp"))
river <- st_read(here("data/spatial_data/NH_Designated_Rivers_24K/NH_Designated_Rivers_24K.shp"))
wetland <- st_read(here("data/spatial_data/National_Wetlands_Inventory_Plus__(NWI_Plus)/National_Wetlands_Inventory_Plus__(NWI_Plus).shp"))

# Convert objects to sf objects and points to coordinates
state.sf <- st_as_sf(USA)
stream.sf <- st_transform(stream, crs = 4326) # Transform to a latitude and longitude coordinate system
river.sf <- st_transform(river, crs = 4326) # Transform to a latitude and longitude coordinate system
wetland.sf <- st_transform(wetland, crs = 4326) %>% # Transform to a latitude and longitude coordinate system
  filter(Landscape == "LK") # Restrict data only to lakes

# Bounding box for cropping the map to a specific region
box <- c(ymin = 43, ymax = 43.27, xmin = -71, xmax = -70.7)

# Limit spatial data to the bounding box
state_crop <- st_crop(state.sf, box)
river_crop <- st_crop(river.sf, box)
stream_crop <- st_crop(stream.sf, box)
# wetland_crop <- st_crop(wetland.sf, box)


# Create a map ----
library(ggspatial)
library(ggrepel)

## Great Bay Estuary map ----
location.map <-
  ggplot() +
  geom_sf(data = state_crop, fill = "gray97") +
  geom_sf(data = stream_crop, col = "blue", size = 1, alpha = 0.2) +
  geom_sf(data = river_crop, col = "blue", size = 1, alpha = 0.2) +
  geom_sf(data = wetland.sf, fill = "blue", alpha = 0.2) +
  geom_point(
    data = site.map,
    aes(x = longitude, y = latitude), fill = "gold1",
    size = 4, pch = 21, inherit.aes = FALSE
  ) +
  scale_x_continuous(limits = c(-71, -70.7), expand = c(0, 0)) +
  scale_y_continuous(limits = c(43, 43.27), expand = c(0, 0)) +
  annotation_scale(location = "tr", width_hint = 0.4, bar_cols = c("grey60", "white")) +
  annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(0.05, "in"), pad_y = unit(0.4, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  labs(y = "Latitude", x = "Longitude") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.position = "none",
    panel.background = element_rect(fill = alpha("blue", 0.2))
  )


## Inset map of New England ----
library(ggthemes)

inset <- ggplot() +
  geom_sf(data = usa, fill = "gray97") +
  xlim(c(74, 69)) +
  ylim(c(40.7, 45.5)) +
  geom_rect(aes(xmin = 71, xmax = 70.7, ymin = 43, ymax = 43.7),
    color = "red", fill = NA
  ) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey20", linetype = 1, linewidth = 1),
    plot.margin = unit(c(0, 0, -1, -1), "mm")
  )


## Plot Great Great Bay with New England inset
library(cowplot)

final.map <- ggdraw() +
  draw_plot(location.map) +
  draw_plot(inset, height = 0.2, x = -0.18, y = 0.77)

# Save plot
ggsave(
  filename = paste0("sampling_location_map.jpeg"),
  plot = final.map,
  path = here("sample_optimization/figures/Minimum Effort, Maximum Insight Manuscript")
)
