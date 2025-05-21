# Plot model dissimilarities for each environmental parameter in Great Bay Estuary

pacman::p_load(
  plyr,
  dplyr,
  here,
  tidyverse,
  vegan,
  sf,
  ggspatial,
  ggrepel,
  viridis,
  patchwork
)

# Source scripts to read specific lines ----
source2 <- function(file, start, end, ...) {
  file.lines <- scan(file, what = character(), skip = start - 1, nlines = end - start + 1, sep = "\n")
  file.lines.collapsed <- paste(file.lines, collapse = "\n")
  source(textConnection(file.lines.collapsed), ...)
}

# Load data ----

# Load data and minimum site functions
source2(here("sample_optimization/scripts/modeling_minimumSite_trendMagnitude.R"), start = 1, end = 40)

# Load site usability data
source2(here("sample_optimization/scripts/modeling_minimumSite_trendMagnitude.R"), start = 81, end = 94)

# Pull only site usability data at the 80% representativeness threshold
all_slope_diff.80 <- Data.Usability %>%
  mutate(all_slope_differences = map(site_contribution.80, ~ .x$all_slope_differences)) %>%
  unnest(all_slope_differences) %>%
  select(variable_units, all_slope_differences)

# Filter to individual sites and their slope differences
all_slope_diff.80$all_slope_differences <- all_slope_diff.80$all_slope_differences %>%
  map(~ .x %>%
    filter(N_sites == 9) %>%
    select(site = site_removed, slope_diff) %>%
    mutate(site = clean_vec(site)))

# Standardize slopes to make them unitless for comparison across environmental parameters
all_slope_diff.80 <- all_slope_diff.80 %>%
  mutate(slope_diff_scaled = map(all_slope_differences, ~ {
    .x %>%
      mutate(
        slope_diff_z = (slope_diff - mean(slope_diff, na.rm = TRUE)) / sd(slope_diff, na.rm = TRUE)
      )
  }))

# Re-format variable names
all_slope_diff.80$variable_units_formatted <- fct_recode(all_slope_diff.80$variable_units,
  "Dissolved Oxygen (mg/L)" = "Dissolved Oxygen_mg/l",
  "Dissolved Oxygen Saturation (%)" = "Dissolved Oxygen Saturation_%",
  "Nitrite + Nitrate, dissolved (mg/L)" = "Nitrogen, Nitrite (No2) + Nitrate (No3) As N Diss_mg/l",
  "Suspended Solids (mg/L)" = "Solids, Suspended Total_mg/l",
  "Water Temperature (°C)" = "Temperature Water_deg c"
)

# Load site coordinates
site.loc <- read_csv(here("data/metadata/PREP_sampling_locations.csv"))

# Filter to Great Bay Estuary sites
site.map <- site.loc %>%
  filter(site %in% c(
    "great_bay_estuary_great_bay",
    "great_bay_estuary_adams_point",
    "winnicut_river_at_rt33_bridge",
    "great_bay_estuary_squamscott_river",
    "great_bay_estuary_lamprey_river",
    "rte_108",
    "great_bay_estuary_oyster_river",
    "rte_108_bridge",
    "rte_9_bridge_central_ave",
    "rte_4_bridge"
  ))

# Calculate the lower bound z-score (slope difference = 0)
z_score_lower_bounds <- all_slope_diff.80 %>%
  mutate(lower_bound_z = map_dbl(all_slope_differences, ~ {
    mean_diff <- mean(.x$slope_diff, na.rm = TRUE)
    sd_diff <- sd(.x$slope_diff, na.rm = TRUE)
    z_min <- -mean_diff / sd_diff
    return(z_min)
  })) %>%
  select(variable_units_formatted, lower_bound_z)

# View results
print(z_score_lower_bounds)

# Merge lower bounds into main object
all_slope_diff.80 <- all_slope_diff.80 %>%
  left_join(z_score_lower_bounds, by = "variable_units_formatted")

# Set up the map ----

# Extract shapefile for USA
states <- readRDS(here("data/spatial_data/gadm36_USA_0_sp.rds"))

# Download New Hampshire spatial data for rivers, streams and wetlands
# https://www.nhgeodata.unh.edu/search?tags=hydrography
stream <- st_read(here("data/spatial_data/New_Hampshire_Stream_Order_Dataset/New_Hampshire_Stream_Order_Dataset.shp"))
wetland <- st_read(here("data/spatial_data/National_Wetlands_Inventory_Plus__(NWI_Plus)/National_Wetlands_Inventory_Plus__(NWI_Plus).shp"))

# Convert objects to sf objects and points to coordinates
state.sf <- st_as_sf(states)
stream.sf <- st_transform(stream, crs = 4326) # Transform to a latitude and longitude coordinate system
wetland.sf <- st_transform(wetland, crs = 4326) %>% # Transform to a latitude and longitude coordinate system
  filter(Landscape == "LK") # Restrict data only to lakes

# Bounding box for cropping the map to a specific region
box <- c(ymin = 43, ymax = 43.27, xmin = -71, xmax = -70.7)

# Limit spatial data to the bounding box
state_crop <- st_crop(state.sf, box)
stream_crop <- st_crop(stream.sf, box)

# Create a water mask layer (inverse of land) to hide stream lines inside major water bodies

# Define bounding box for the study area
bbox <- st_as_sf(st_sfc(st_polygon(list(rbind(
  c(-71, 43), c(-71, 43.27), c(-70.7, 43.27), c(-70.7, 43), c(-71, 43)
))), crs = st_crs(state.sf)))

# Create a water mask by subtracting land from the bounding box
water_mask <- st_difference(bbox, st_union(state.sf))

# Create a spatial map ----

# Define map plotting function using precomputed dissimilarity
create_slopeDiff_map <- function(slope_diff_z, param_label, lower_bound_z) {

  # Join slope differences to site.map
  site.map.diss <- left_join(site.map, slope_diff_z, by = "site")

  # Plot
  p <-
    ggplot() +
    scale_x_continuous(limits = c(-71, -70.7), expand = c(0, 0)) +
    scale_y_continuous(limits = c(43, 43.27), expand = c(0, 0)) +
    geom_sf(data = state.sf, fill = "gray97") +
    geom_sf(data = wetland.sf, fill = "blue", alpha = 0.2) +
    geom_sf(data = stream_crop, col = "blue", size = 1, alpha = 0.2) +
    geom_sf(data = water_mask, fill = "#CBCBFF") +
    geom_point(
      data = site.map.diss,
      aes(x = longitude, y = latitude, fill = slope_diff_z),
      size = 4, pch = 21, inherit.aes = FALSE
    ) +
    scale_fill_viridis_c(
      option = "plasma",
      limits = c(-1.5, 3),
      name = paste0("Slope Difference\n(System z = ", round(lower_bound_z, 2), ")")
    ) +
    annotation_scale(location = "tr", width_hint = 0.4) +
    annotation_north_arrow(
      location = "tr", which_north = "true",
      pad_x = unit(0.05, "in"), pad_y = unit(0.4, "in"),
      style = north_arrow_fancy_orienteering
    ) +
    labs(x = "Longitude", y = "Latitude", title = param_label) +
    ggtitle(param_label) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, vjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.key = element_blank(),
      legend.background = element_blank()
    )

  return(p)
}

# Generate maps with added annotation
plot_list <- pmap(
  list(
    slope_diff_z = all_slope_diff.80$slope_diff_scaled,
    param_label = all_slope_diff.80$variable_units_formatted,
    lower_bound_z = all_slope_diff.80$lower_bound_z
  ),
  create_slopeDiff_map
)

# Assemble into a grid with labels A–E
grid_plot <- wrap_plots(plotlist = plot_list, ncol = 2) +
  plot_annotation(tag_levels = "A")

# Save the combined grid
ggsave(
  filename = "model-dissimilarity-maps.png",
  plot = grid_plot,
  path = here("sample_optimization/figures/Minimum Effort, Maximum Insight Manuscript"),
  dpi = 600,
  width = 16,
  height = 18,
  units = "in"
)
