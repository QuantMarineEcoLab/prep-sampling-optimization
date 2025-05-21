# Run spatial reduction on Great Bay environmental variables to determine the minimum # of sites needed to detect trends and which site yields information

# Load required packages
pacman::p_load(
  plyr,
  dplyr,
  purrr,
  here,
  ggplot2,
  tidyr,
  forcats,
  cowplot,
  car
)

# Clean data ----

# Load data
CleanData <- read.csv(here("data/Great Bay Environmental Data/clean_imputedPMM_Great Bay_environmental data.csv"))

# Ensure year is numeric
CleanData$year <- as.numeric(CleanData$year)

# Filter data to remove any variables that have less than 10 sites
Data <- CleanData %>%
  select(site = sitename, variable_units, year, response = datavalue) %>% # Rename and select relevant columns
  group_by(variable_units) %>% # Group by variable_units
  filter(n_distinct(site) >= 10) %>% # Filter to retain only variables with at least 10 sites
  filter(site != "RTE 108 BRIDGE/MILL POND") # Remove Rt. 108 Bridge Mill Pond because water temperature was not sampled at this location

# Organize each site's data into one place by creating a nested data frame
Data <- Data %>%
  group_by(variable_units) %>%
  tidyr::nest()

# Site Minimum ----

# Load site minimum functions
source(here("sample_optimization/scripts/minimumSite_functions.R"))

# Run minimum site analyses (magnitude only) for 80% of the subsets matching the true slope
Data.80 <- Data %>%
  mutate(multiple_breakups = map(
    data,
    function(data) multiple_breakups(data)
  )) %>%
  mutate(min_site_mag = map(
    data,
    function(data) min_site_magnitude(data, min_percent = 80, error_multiplier = 1)
  ))

# Create a summary table with minimum sites for trend magnitude
minSite_summary.80 <- Data.80 %>%
  select(variable = variable_units, min_site_mag) %>%
  mutate(
    min_site_mag = map_dbl(min_site_mag, ~ unique(.$minimum_site))
  ) %>%
  select(variable, min_site_mag)


# Run minimum site analyses (magnitude only) for 100% of the subsets matching the true slope
Data.100 <- Data %>%
  mutate(multiple_breakups = map(
    data,
    function(data) multiple_breakups(data)
  )) %>%
  mutate(min_site_mag = map(
    data,
    function(data) min_site_magnitude(data, min_percent = 100, error_multiplier = 1)
  ))

# Create a summary table with minimum sites for trend magnitude
minSite_summary.100 <- Data.100 %>%
  select(variable = variable_units, min_site_mag) %>%
  mutate(
    min_site_mag = map_dbl(min_site_mag, ~ unique(.$minimum_site))
  ) %>%
  select(variable, min_site_mag)


# Site Usability ----

Data.Usability <- Data %>%
  mutate(
    site_contribution.80 = map(
      data,
      ~ usable_sites(.x, min_percent = 80, error_multiplier = 1)
    ),
    site_contribution.100 = map(
      data,
      ~ usable_sites(.x, min_percent = 100, error_multiplier = 1)
    )
  )

# Variance Partitioning ----

# Calculate variance components
variable.nest <- Data %>%
  mutate(
    min_site_mag = map(
      data,
      function(data) min_site_magnitude(data, min_percent = 80, error_multiplier = 1)
    ),
    slope = map_dbl(
      data,
      function(df) {
        stand.data <- standardize(df)
        linefit_output <- linefit_with_site(stand.data)
        slope <- as.numeric(linefit_output$slope) # safer than hardcoding [4]
        return(slope)
      }
    ),
    slope_se = map_dbl(
      data,
      function(df) {
        stand.data <- standardize(df)
        linefit_output <- linefit_with_site(stand.data)
        slope_se <- as.numeric(linefit_output$slope_se)
        return(slope_se)
      }
    ),
    SD = map_dbl(
      data,
      function(df) {
        sd_value <- sd(df$response, na.rm = TRUE)
        return(sd_value)
      }
    ),
    autocorrelation = map_dbl(
      data,
      function(df) {
        # Ensure data is sorted by year
        df[order(df$year), ]

        # Calculate autocorrelation without plotting
        acf_results <- acf(df$response, plot = FALSE)

        # Extract the lag 1 autocorrelation value
        lag_1_autocorr <- acf_results$acf[2] # Lag 1 is the second element (Lag 0 is the first)

        # Print the autocorrelation value
        return(lag_1_autocorr)
      }
    )
  )

# Unnest the magnitude stability results and combine with relevant statistics
variable.min.mag <- variable.nest %>%
  select(variable_units, min_site_mag, slope, slope_se, SD, autocorrelation) %>%
  unnest(cols = c(min_site_mag))

# Identify the first window length that meets 80% and 100% magnitude stability criteria
variable.min.mag.1 <- variable.min.mag %>%
  group_by(variable_units) %>%
  mutate(stability_time_80_test = if_else(proportion_correct >= 80, minimum_site, NA_integer_)) %>%
  mutate(stability_time_80 = min(stability_time_80_test, na.rm = TRUE)) %>%
  mutate(stability_time_100_test = if_else(proportion_correct >= 100, minimum_site, NA_integer_)) %>%
  mutate(stability_time_100 = min(stability_time_100_test, na.rm = TRUE))

# Clean and relabel output table for clarity
variable.min.mag.2 <- variable.min.mag.1 %>%
  select(
    Variable = variable_units,
    "Minimum Site (80% correct)" = stability_time_80,
    "Minimum Site (100% correct)" = stability_time_100,
    slope,
    slope_se,
    SD,
    autocorrelation
  ) %>%
  unique()

# Subset to five key variables for publication/presentation
variable.min.mag.pub <- variable.min.mag.2 %>%
  filter(Variable == "Dissolved Oxygen_mg/l" |
    Variable == "Dissolved Oxygen Saturation_%" |
    Variable == "Nitrogen, Nitrite (No2) + Nitrate (No3) As N Diss_mg/l" |
    Variable == "Solids, Suspended Total_mg/l" |
    Variable == "Temperature Water_deg c")

# Recode variable names for clean labels in publication
variable.min.mag.pub$Variable <- fct_recode(variable.min.mag.pub$Variable,
  "Dissolved Oxygen (mg/L)" = "Dissolved Oxygen_mg/l",
  "Dissolved Oxygen Saturation (%)" = "Dissolved Oxygen Saturation_%",
  "Nitrite + Nitrate, dissolved (mg/L)" = "Nitrogen, Nitrite (No2) + Nitrate (No3) As N Diss_mg/l",
  "Suspended Solids (mg/L)" = "Solids, Suspended Total_mg/l",
  "Water Temperature (°C)" = "Temperature Water_deg c"
)

## Run regression models ----

### All Variables ----

# Fit a linear model for 100% correct based on predictors
var_min_mag_model <- lm(`Minimum Site (100% correct)` ~ autocorrelation + slope_se + SD, data = variable.min.mag.pub)
summary(var_min_mag_model)

# Compute ANOVA and percentage of variance explained
anova_results <- Anova(var_min_mag_model)
total_ss <- sum(anova_results$`Sum Sq`)
variance_explained <- (anova_results$`Sum Sq` / total_ss) * 100

# Summarize into a clean table
minSite_mag_table <- data.frame(
  Term = rownames(anova_results),
  Variance_Explained_Percentage = variance_explained
)


# Summary Table: Variance Partitioning ---

# Extract variance partitioning and format
variance_partition_table <- minSite_mag_table %>%
  select(Metric = Term, `% Variance Explained` = Variance_Explained_Percentage)

# Rename and format metrics
variance_partition_table$Metric <- fct_recode(variance_partition_table$Metric,
  "Slope SE" = "slope_se",
  "Data SD" = "SD",
  "Autocorrelation" = "autocorrelation"
)

# Round variance %s
variance_partition_table <- variance_partition_table %>%
  dplyr::mutate(across(-Metric, ~ paste0(round(.x, 1), "%")))
