# Create figures for minimum site across Great Bay for each environmental variable

# Load required packages
pacman::p_load(
  plyr,
  dplyr,
  tidyr,
  here,
  purrr,
  ggplot2,
  ggspatial,
  sf,
  viridis,
  cowplot,
  forcats,
  htmltools,
  webshot,
  stringr,
  formattable
)

# Minimum Site Table (80% and 100% correct) ----

# Load data
source(here("sample_optimization/scripts/modeling_minimumSite_trendMagnitude.R"))

# Extract minimum site (100% correct) and label threshold
variable.min.100 <- minSite_summary.100 %>%
  select(variable, `100% correct` = min_site_mag)

# Extract minimum site (80% correct) and label threshold
variable.min.80 <- minSite_summary.80 %>%
  select(variable, `80% correct` = min_site_mag)

# Combine both datasets into one long table
variable.min <- full_join(
  variable.min.80,
  variable.min.100
)

# Re-name threshold columns
variable.min.2 <- variable.min %>%
  select(
    Variable = variable,
    "Minimum Sites (80% correct)" = `80% correct`,
    "Minimum Sites (100% correct)" = `100% correct`
  )

# For publication, only extract 5 variables
variable.min.pub <- variable.min.2 %>%
  filter(Variable == "Dissolved Oxygen Saturation_%" |
           Variable == "Dissolved Oxygen_mg/l" |
           Variable == "Nitrogen, Nitrite (No2) + Nitrate (No3) As N Diss_mg/l" |
           Variable == "Solids, Suspended Total_mg/l" |
           Variable == "Temperature Water_deg c")

# Clean up names
variable.min.pub$Variable <- fct_recode(variable.min.pub$Variable,
                                        "Dissolved Oxygen Saturation (%)" = "Dissolved Oxygen Saturation_%",
                                        "Dissolved Oxygen (mg/L)" = "Dissolved Oxygen_mg/l",
                                        "Nitrite + Nitrate, dissolved (mg/L)" = "Nitrogen, Nitrite (No2) + Nitrate (No3) As N Diss_mg/l",
                                        "Suspended Solids (mg/L)" = "Solids, Suspended Total_mg/l",
                                        "Water Temperature (°C)" = "Temperature Water_deg c"
)

# Create a table
min_site.table <- formattable(variable.min.pub,
                              align = c("l", "l", "l"),
                              list(area(col = c("Minimum Sites (80% correct)", "Minimum Sites (100% correct)"))
                                   ~ normalize_bar("pink", 0.2))
)

# Function to export minimum site table
export_formattable <- function(f, file, width = "100%", height = NULL,
                               background = "white", delay = 0.2) {
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay,
          zoom = 5
  )
}

## Export minimum site table ----
export_formattable(min_site.table, here("sample_optimization/figures/Minimum Effort, Maximum Insight Manuscript/minimumSite_variable_acrossYears.png"))


# Site Usability Functions ----

# Function to clean up site names
clean_site_name <- function(x) {
  x <- as.character(x)
  
  x_cleaned <- x %>%
    gsub("\\s*,\\s*", ", ", .) %>%     # Normalize commas
    gsub(",\\s*$", "", .) %>%          # Remove trailing comma
    gsub("[()]", "", .) %>%            # Remove parentheses
    gsub("[^a-zA-Z0-9,]+", "_", .) %>% # Replace non-alphanum (except commas) with _
    gsub("_+", "_", .) %>%             # Collapse multiple underscores
    gsub("_?,_?", ", ", .) %>%         # Ensure commas have space after
    gsub("^_|_$", "", .) %>%           # Trim leading/trailing underscores
    tolower()                          # Convert to lowercase
  
  return(x_cleaned)
}

# Function to export site usability table
export_formattable_plain <- function(df, file, width = "100%", height = NULL,
                                     background = "white", delay = 0.2) {
  library(formattable)
  library(htmltools)
  library(webshot)
  library(htmlwidgets)
  
  # Create formattable object
  f <- formattable(df, list())
  
  # Convert to HTML widget
  w <- as.htmlwidget(f, width = width, height = height)
  
  # Custom CSS: wider columns 1 and 4, narrower column 2 with wrap (no centering)
  style_tag <- tags$style(HTML("
    .formattable_widget table th:nth-child(1),
    .formattable_widget table td:nth-child(1),
    .formattable_widget table th:nth-child(4),
    .formattable_widget table td:nth-child(4) {
      min-width: 200px !important;
      white-space: nowrap !important;
    }

    .formattable_widget table th:nth-child(2),
    .formattable_widget table td:nth-child(2) {
      min-width: 100px !important;
      white-space: normal !important;
    }
  "))
  
  # Combine widget and styles
  w <- htmltools::browsable(
    tagList(style_tag, w)
  )
  
  # Render HTML to temp file
  path <- html_print(w, background = background, viewer = NULL)
  
  # Take screenshot as PNG
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay,
          zoom = 5)
}


# Site Usability Table (80% correct) ----

# Generate a summary table
site_usability_summary.80 <- Data.Usability %>%
  
  transmute(
    variable_units,
    
    # Extract: Minimum sites needed
    minimum_sites_required = map_dbl(site_contribution.80, ~ .x$minimum_sites_required),
    
    # Extract: Best combo at minimum site count
    smallest_difference_min_site = map(site_contribution.80, ~ .x$smallest_difference_min_site[[1]]),
    
    # Extract: Full model (slope = baseline)
    full_model = map(site_contribution.80, ~ {
      diffs <- .x$all_slope_differences[[1]]
      diffs %>%
        filter(N_sites == 10)
    })
  ) %>%
  
  # Flatten everything into rows
  unnest(cols = c(smallest_difference_min_site, full_model), names_sep = "_", keep_empty = TRUE) %>%
  
  # Select and rename only the needed columns
  select(
    'Variable' = variable_units,
    'Minimum # of Sites' = minimum_sites_required,
    'Total Sampled Sites' = full_model_N_sites,
    'Sites Removed' = smallest_difference_min_site_site_removed,
    'Subsample Slope' = smallest_difference_min_site_year_slope,
    'Full Sample Slope' = full_model_year_slope,
    'Slope Difference' = smallest_difference_min_site_slope_diff
  ) %>%
  
  # Round numeric columns to 2 decimal places
  mutate(
    `Subsample Slope` = round(`Subsample Slope`, 6),
    `Full Sample Slope` = round(`Full Sample Slope`, 6),
    `Slope Difference` = round(`Slope Difference`, 6)
  )

# Bring in site names
site_names <- read.csv("data/metadata/PREP_sampling_locations.csv")

# Replace comma and everything after with parentheses
site_names <- site_names %>%
  mutate(
    site_formatted = if_else(
      str_detect(site_formatted, ","),
      str_replace(site_formatted, "\\s*,\\s*", " ("),  # Replace comma with opening parenthesis
      site_formatted
    ),
    site_formatted = if_else(
      str_detect(site_formatted, "\\("),
      paste0(site_formatted, ")"),  # Add closing parenthesis if opening one exists
      site_formatted
    )
  )

# Clean the site names in the "Sites Removed" column
site_usability_summary.80$`Sites Removed` <- clean_site_name(site_usability_summary.80$`Sites Removed`)

# Create lookup table from raw site names to formatted labels
site_lookup <- site_names %>%
  select(site, site_formatted) %>%
  mutate(
    site = clean_site_name(site),           # Ensure names match cleaned format
    site_formatted = as.character(site_formatted)
  )

# Create named vector for replacing cleaned codes with display names
site_map <- setNames(site_lookup$site_formatted, site_lookup$site)

# Replace all site codes in 'Sites Removed' with readable names
site_usability_summary.80 <- site_usability_summary.80 %>%
  mutate(
    `Sites Removed` = str_replace_all(`Sites Removed`, site_map)
  )

# Clean up variable names
site_usability_summary.80$Variable <- fct_recode(site_usability_summary.80$Variable,
                                        "Dissolved Oxygen Saturation (%)" = "Dissolved Oxygen Saturation_%",
                                        "Dissolved Oxygen (mg/L)" = "Dissolved Oxygen_mg/l",
                                        "Nitrite + Nitrate, dissolved (mg/L)" = "Nitrogen, Nitrite (No2) + Nitrate (No3) As N Diss_mg/l",
                                        "Suspended Solids (mg/L)" = "Solids, Suspended Total_mg/l",
                                        "Water Temperature (°C)" = "Temperature Water_deg c"
)

# Create a line break for "Minimum # of Sites" so "Of Sites" is on the second row
names(site_usability_summary.80)[names(site_usability_summary.80) == "Minimum # of Sites"] <- "Minimum #<br>of Sites"

## Export site usability table ----
export_formattable_plain(
  site_usability_summary.80,
  here("sample_optimization/figures/Minimum Effort, Maximum Insight Manuscript/site_usability_summary_80.png")
)


# One-Site Removal Difference Table (80% correct) ----

one_site_diff_summary.80 <- Data.Usability %>%
  
  transmute(
    variable_units,
    
    # Extract slope difference info for largest & smallest one-site removal
    largest_diff_one_site = map(site_contribution.80, ~ .x$largest_difference_one_site_removed[[1]]),
    smallest_diff_one_site = map(site_contribution.80, ~ .x$smallest_difference_one_site_removed[[1]])
  ) %>%
  
  # Flatten
  unnest(cols = c(largest_diff_one_site, smallest_diff_one_site), names_sep = "_", keep_empty = TRUE) %>%
  
  # Select and rename
  select(
    'Variable' = variable_units,
    'Most Useful Site' = largest_diff_one_site_site_removed,
    'Least Useful Site' = smallest_diff_one_site_site_removed
  )

## Clean site names for 1-site removal ----

# Bring in site names
site_names <- read.csv("data/metadata/PREP_sampling_locations.csv")

# Replace comma and everything after with parentheses
site_names <- site_names %>%
  mutate(
    site_formatted = if_else(
      str_detect(site_formatted, ","),
      str_replace(site_formatted, "\\s*,\\s*", " ("),  # Replace comma with opening parenthesis
      site_formatted
    ),
    site_formatted = if_else(
      str_detect(site_formatted, "\\("),
      paste0(site_formatted, ")"),  # Add closing parenthesis if opening one exists
      site_formatted
    )
  )

# Clean site names
one_site_diff_summary.80$`Most Useful Site` <- clean_site_name(one_site_diff_summary.80$`Most Useful Site`)
one_site_diff_summary.80$`Least Useful Site` <- clean_site_name(one_site_diff_summary.80$`Least Useful Site`)

# Create lookup table from raw site names to formatted labels
site_lookup <- site_names %>%
  select(site, site_formatted) %>%
  mutate(
    site = clean_site_name(site),           # Ensure names match cleaned format
    site_formatted = as.character(site_formatted)
  )

# Create named vector for replacing cleaned codes with display names
site_map <- setNames(site_lookup$site_formatted, site_lookup$site)

# Replace codes with human-readable names
one_site_diff_summary.80 <- one_site_diff_summary.80 %>%
  mutate(
    `Most Useful Site` = str_replace_all(`Most Useful Site`, site_map),
    `Least Useful Site` = str_replace_all(`Least Useful Site`, site_map)
  )

# Clean up variable names
one_site_diff_summary.80$Variable <- fct_recode(one_site_diff_summary.80$Variable,
                                                "Dissolved Oxygen Saturation (%)" = "Dissolved Oxygen Saturation_%",
                                                "Dissolved Oxygen (mg/L)" = "Dissolved Oxygen_mg/l",
                                                "Nitrite + Nitrate, dissolved (mg/L)" = "Nitrogen, Nitrite (No2) + Nitrate (No3) As N Diss_mg/l",
                                                "Suspended Solids (mg/L)" = "Solids, Suspended Total_mg/l",
                                                "Water Temperature (°C)" = "Temperature Water_deg c"
)

## Export as formattable ----

export_formattable_plain(
  one_site_diff_summary.80,
  here("sample_optimization/figures/Minimum Effort, Maximum Insight Manuscript/one_site_diff_summary_80.png")
)




# Create spatial maps ----

# Data processing for calculating z-scores for each site
Data <- Data %>%
  mutate(
    # Clean site names for each dataset
    data = map(data, function(df) {
      df %>%
        mutate(site = clean_vec(site))  # Clean site names using clean_vec function
    }),
    
    # Calculate z-scores for each site based on response variable differences
    z_stats = map(data, function(df) {
      df %>%
        group_by(site) %>%
        arrange(year) %>%
        mutate(
          difference = response - lag(response),  # Calculate difference between consecutive years
        ) %>%
        filter(!is.na(difference)) %>%
        mutate(
          median_diff = median(difference, na.rm = TRUE),  # Median of differences
          sd_diff = sd(difference, na.rm = TRUE),           # Standard deviation of differences
          z_score = abs((difference - median_diff) / sd_diff)  # Calculate z-scores
        )
    }),
    
    # Calculate the median z-score for each site across datasets
    site_median_z_score = map(z_stats, function(df) {
      df %>%
        group_by(site) %>%
        summarize(site_median_z_score = median(z_score, na.rm = TRUE), .groups = 'drop')
    }),
    
    # Calculate overall median z-score across sites for each dataset
    overall_median_z_score = map_dbl(site_median_z_score, function(df) {
      median(df$site_median_z_score, na.rm = TRUE)  # Calculate the median z-score across sites
    })
  )

# Clean up variable names to standardize units
Data$variable_units <- fct_recode(Data$variable_units,
                                  "Dissolved Oxygen (mg/L)" = "Dissolved Oxygen_mg/l",
                                  "Dissolved Oxygen Saturation (%)" = "Dissolved Oxygen Saturation_%",
                                  "Nitrite + Nitrate, dissolved (mg/L)" = "Nitrogen, Nitrite (No2) + Nitrate (No3) As N Diss_mg/l",
                                  "Suspended Solids (mg/L)" = "Solids, Suspended Total_mg/l",
                                  "Water Temperature (°C)" = "Temperature Water_deg c"
)

# Unnest the dataframe to extract site_median_z_score from the nested structure
unnested_data <- Data %>%
  unnest(cols = c(site_median_z_score))

# Calculate the global min and max values for site_median_z_score across all variables
z_min <- min(unnested_data$site_median_z_score, na.rm = TRUE)
z_max <- max(unnested_data$site_median_z_score, na.rm = TRUE)

# Prepare spatial data for mapping
Data <- Data %>%
  mutate(
    spatial.map = map(site_median_z_score, function(data) {
      
      # Bring in site coordinates from metadata
      site.loc <- read.csv(here("data/metadata/PREP_sampling_locations.csv"))
      
      # Join site data with location coordinates
      site.map <- left_join(data, site.loc)
      
      # Ensure longitude and latitude are numeric
      site.map$longitude <- as.numeric(as.character(site.map$longitude))
      site.map$latitude <- as.numeric(as.character(site.map$latitude))
      
      # Extract shapefile for the USA and NH data
      USA <- readRDS(here("data/spatial_data/gadm36_USA_0_sp.rds"))
      stream <- st_read(here("data/spatial_data/New_Hampshire_Stream_Order_Dataset/New_Hampshire_Stream_Order_Dataset.shp"))
      river <- st_read(here("data/spatial_data/NH_Designated_Rivers_24K/NH_Designated_Rivers_24K.shp"))
      wetland <- st_read(here("data/spatial_data/National_Wetlands_Inventory_Plus__(NWI_Plus)/National_Wetlands_Inventory_Plus__(NWI_Plus).shp"))
      
      # Transform spatial data to correct coordinate reference system (CRS)
      state.sf <- st_as_sf(USA)
      stream.sf <- st_transform(stream, crs = 4326)  # Transform to lat/lon
      river.sf <- st_transform(river, crs = 4326)    # Transform to lat/lon
      wetland.sf <- st_transform(wetland, crs = 4326) %>%
        filter(Landscape == "LK")  # Filter wetland data to lakes only
      
      # Define bounding box for cropping
      box <- c(ymin = 43, ymax = 43.27, xmin = -71, xmax = -70.7)
      
      # Crop spatial data to the region of interest
      state_crop <- st_crop(state.sf, box)
      river_crop <- st_crop(river.sf, box)
      stream_crop <- st_crop(stream.sf, box)
      
      # Plot spatial map with spatial data layers
      spatial.plot <- ggplot() +
        geom_sf(data = state_crop, fill = "gray97") +
        geom_sf(data = stream_crop, col = "blue", size = 1, alpha = 0.2) +
        geom_sf(data = river_crop, col = "blue", size = 1, alpha = 0.2) +
        geom_sf(data = wetland.sf, col = "blue", alpha = 0.2) +
        geom_point(data = site.map, aes(x = longitude, y = latitude, fill = site_median_z_score),
                   size = 4, pch = 21, inherit.aes = FALSE) +
        scale_fill_viridis_c(limits = c(z_min, z_max), option = "viridis") +  # Use global min/max for all variables
        scale_x_continuous(limits = c(-71, -70.7), expand = c(0, 0)) +
        scale_y_continuous(limits = c(43, 43.27), expand = c(0, 0)) +
        ggtitle(variable_units) +
        annotation_scale(location = "tr", width_hint = 0.4, bar_cols = c("grey60", "white")) +
        annotation_north_arrow(location = "tr", which_north = "true",
                               pad_x = unit(0.05, "in"), pad_y = unit(0.4, "in"),
                               style = north_arrow_fancy_orienteering) +
        labs(y = "Latitude", x = "Longitude", fill = bquote(tilde(bold(z)))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 14, face = "bold", hjust = 0),
              axis.title = element_text(size = 14, face = "bold"),
              axis.text = element_text(size = 12),
              axis.text.x = element_text(angle = 45, vjust = 0.5),
              panel.background = element_rect(fill = alpha("blue", 0.2)))
      
      return(spatial.plot)  # Return the map plot
    })
  )

# Check one map to make sure it looks good
Data$spatial.map[[5]]

# Filter and store individual plots for each variable
DO_mgL <- Data %>%
  filter(variable_units == 'Dissolved Oxygen (mg/L)')

DO_sat <- Data %>%
  filter(variable_units == 'Dissolved Oxygen Saturation (%)')

nitrogen <- Data %>%
  filter(variable_units == 'Nitrite + Nitrate, dissolved (mg/L)')

solids <- Data %>%
  filter(variable_units == 'Suspended Solids (mg/L)')

w_temp <- Data %>%
  filter(variable_units == 'Water Temperature (°C)')

# Combine the plots into a grid layout
grid.plot <- plot_grid(
  DO_mgL$spatial.map[[1]], DO_sat$spatial.map[[1]], 
  nitrogen$spatial.map[[1]], solids$spatial.map[[1]], 
  w_temp$spatial.map[[1]],
  ncol = 2,  # 2 columns for the grid layout
  align = 'hv',  # Align horizontally and vertically
  labels = c("A", "B", "C", "D", "E"),  # Custom labels for each plot
  label_size = 15  # Adjust label size
)

# Save the combined plot as an image file
ggsave(filename = paste0("map_variable_Zscores.jpeg"),
       plot = grid.plot,
       path = here("sample_optimization/figures/Minimum Effort, Maximum Insight Manuscript/"),
       width = 11,
       height = 13,
       dpi = 800
)
