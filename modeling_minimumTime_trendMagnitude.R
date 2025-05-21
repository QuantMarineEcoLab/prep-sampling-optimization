# Determine minimum years needed to detect trend magnitude in Great Bay environmental data.

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


# Minimum Time ----

# Load minimum time functions
source(here("sample_optimization/scripts/minimumTime_functions.R"))

# Run minimum time analyses (magnitude) for 80% of the subsets matching the true slope
Data.80 <- Data %>%
  mutate(multiple_breakups = map(
    data,
    function(data) multiple_breakups(data)
  )) %>%
  mutate(min_time_mag = map(
    data,
    function(data) min_time_magnitude(data, min_percent = 80, error_multiplier = 1)
  ))

# Create a summary table with minimum times for trend magnitude
minTime_summary.80 <- Data.80 %>%
  select(variable = variable_units, min_time_mag) %>%
  mutate(
    min_time_mag = map_dbl(min_time_mag, ~ unique(.$stability_time))
  ) %>%
  select(variable, min_time_mag)

# Run minimum time analyses (magnitude) for 100% of the subsets matching the true slope
Data.100 <- Data %>%
  mutate(multiple_breakups = map(
    data,
    function(data) multiple_breakups(data)
  )) %>%
  mutate(min_time_mag = map(
    data,
    function(data) min_time_magnitude(data, min_percent = 100, error_multiplier = 1)
  ))

# Create summary table with stability times for magnitude
minTime_summary.100 <- Data.100 %>%
  select(variable = variable_units, min_time_mag) %>%
  mutate(
    min_time_mag = map_dbl(min_time_mag, ~ unique(.$stability_time))
  ) %>%
  select(variable, min_time_mag)

# Variance Partitioning ----

# Rearrange data
variable.nest <- CleanData %>%
  select(site = sitename, variable_units, year, response = datavalue) %>%
  group_by(variable_units, site) %>%
  tidyr::nest()


# Calculate variance components
variable.nest <- variable.nest %>%
  mutate(
    min_time_mag = map(
      data,
      function(data) min_time_magnitude(data, min_percent = 80, error_multiplier = 1)
    ),
    slope = map_dbl(
      data,
      function(df) {
        stand.data <- standardize(df)
        linefit_output <- linefit_wrapper(stand.data)
        slope <- as.numeric(linefit_output$slope) # safer than hardcoding [4]
        return(slope)
      }
    ),
    slope_se = map_dbl(
      data,
      function(df) {
        stand.data <- standardize(df)
        linefit_output <- linefit_wrapper(stand.data)
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
  dplyr::select(variable_units, site, min_time_mag, slope, slope_se, SD, autocorrelation) %>%
  unnest(cols = c(min_time_mag))

# Make window_length and stability_time numeric
variable.min.mag$window_length <- as.numeric(variable.min.mag$window_length)
variable.min.mag$stability_time <- as.numeric(variable.min.mag$stability_time)

# Identify the first window length that meets 80% and 100% magnitude stability criteria
variable.min.mag.1 <- variable.min.mag %>%
  group_by(variable_units, site) %>%
  mutate(stability_time_80_test = if_else(proportion_correct >= 80, window_length, NA_real_)) %>%
  mutate(stability_time_80 = min(stability_time_80_test, na.rm = TRUE)) %>%
  mutate(stability_time_100_test = if_else(proportion_correct >= 100, window_length, NA_real_)) %>%
  mutate(stability_time_100 = min(stability_time_100_test, na.rm = TRUE))

# Clean and relabel output table for clarity
variable.min.mag.2 <- variable.min.mag.1 %>%
  dplyr::select(
    Variable = variable_units,
    Site = site,
    "Minimum Years (80% correct)" = stability_time_80,
    "Minimum Years (100% correct)" = stability_time_100,
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

### Visualization: Trend Magnitude vs Metrics ----

# Plot: 100% correct vs slope
min100_slope <- ggplot(variable.min.mag.pub, aes(x = slope, y = `Minimum Years (100% correct)`, color = Variable)) +
  geom_point(size = 3) +
  geom_smooth(alpha = 0.2, se = FALSE) +
  labs(x = "Slope") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.title = element_text(size = 12, face = "bold")
  )

# Plot: 80% correct vs slope
min80_slope <- ggplot(variable.min.mag.pub, aes(x = slope, y = `Minimum Years (80% correct)`)) +
  geom_point(size = 3, aes(color = Variable)) +
  geom_smooth(color = "black", alpha = 0.2) +
  labs(x = "Slope") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.title = element_text(size = 12, face = "bold")
  )

# Plot: 100% correct vs slope SE
min100_slope_se <- ggplot(variable.min.mag.pub, aes(x = slope_se, y = `Minimum Years (100% correct)`, color = Variable)) +
  geom_point(size = 3) +
  geom_smooth(se = FALSE, alpha = 0.2, method = "lm") +
  labs(x = "Slope SE") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.title = element_text(size = 12, face = "bold")
  )

# Plot: 80% correct vs slope SE
min80_slope_se <- ggplot(variable.min.mag.pub, aes(x = slope_se, y = `Minimum Years (80% correct)`)) +
  geom_point(size = 3, aes(color = Variable)) +
  geom_smooth(color = "black", alpha = 0.2) +
  labs(x = "Slope SE") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.title = element_text(size = 12, face = "bold")
  )

# Plot: 100% correct vs standard deviation
min100_SD <- ggplot(variable.min.mag.pub, aes(x = SD, y = `Minimum Years (100% correct)`, color = Variable)) +
  geom_point(size = 3) +
  geom_smooth(se = FALSE, alpha = 0.2, method = "lm") +
  labs(x = "Variability (SD)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.title = element_text(size = 12, face = "bold")
  )

# Plot: 80% correct vs standard deviation
min80_SD <- ggplot(variable.min.mag.pub, aes(x = SD, y = `Minimum Years (80% correct)`)) +
  geom_point(size = 3, aes(color = Variable)) +
  geom_smooth(color = "black", alpha = 0.2) +
  labs(x = "Variability (SD)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.title = element_text(size = 12, face = "bold")
  )

# Plot: 100% correct vs autocorrelation
min100_autoCorr <- ggplot(variable.min.mag.pub, aes(x = autocorrelation, y = `Minimum Years (100% correct)`, color = Variable)) +
  geom_point(size = 3) +
  geom_smooth(se = FALSE, alpha = 0.2, method = "lm") +
  labs(x = "Lag-1 Autocorrelation") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.title = element_text(size = 12, face = "bold")
  )

# Plot: 80% correct vs autocorrelation
min80_autoCorr <- ggplot(variable.min.mag.pub, aes(x = autocorrelation, y = `Minimum Years (80% correct)`)) +
  geom_point(size = 3, aes(color = Variable)) +
  geom_smooth(color = "black", alpha = 0.2) +
  labs(x = "Lag-1 Autocorrelation") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.title = element_text(size = 12, face = "bold")
  )

# Combine all plots into a single grid for publication
grid.plot <- plot_grid(
  min100_slope, min100_slope_se, min100_SD, min100_autoCorr,
  min80_slope, min80_slope_se, min80_SD, min80_autoCorr,
  align = "hv",
  labels = "AUTO", ncol = 4
)

### Run regression models ----

#### All Variables ----

# Fit a linear model for 100% correct based on predictors
var_min_mag_model <- lm(`Minimum Years (100% correct)` ~ slope + slope_se + SD + autocorrelation, data = variable.min.mag.pub)
summary(var_min_mag_model)

# Compute ANOVA and percentage of variance explained
anova_results <- Anova(var_min_mag_model)
total_ss <- sum(anova_results$`Sum Sq`)
variance_explained <- (anova_results$`Sum Sq` / total_ss) * 100

# Summarize into a clean table
minTime_mag_table <- data.frame(
  Term = rownames(anova_results),
  Variance_Explained_Percentage = variance_explained
)

#### Dissolved Oxygen ----

# Filter data
variable.min.pub.DO <- variable.min.mag.pub %>%
  filter(Variable == "Dissolved Oxygen (mg/L)")

# Fit a linear model for 100% correct based on predictors
var_min_model_DO <- lm(`Minimum Years (100% correct)` ~ slope + slope_se + SD + autocorrelation, data = variable.min.pub.DO)
summary(var_min_model_DO)

#### Dissolved Oxygen Saturation ----

# Filter data
variable.min.pub.DOS <- variable.min.mag.pub %>%
  filter(Variable == "Dissolved Oxygen Saturation (%)")

# Fit a linear model for 100% correct based on predictors
var_min_model_DOS <- lm(`Minimum Years (100% correct)` ~ slope + slope_se + SD + autocorrelation, data = variable.min.pub.DOS)
summary(var_min_model_DOS)

#### Water Temperature ----

# Filter data
variable.min.pub.wtemp <- variable.min.mag.pub %>%
  filter(Variable == "Water Temperature (°C)")

# Fit a linear model for 100% correct based on predictors
var_min_model_wtemp <- lm(`Minimum Years (100% correct)` ~ slope + slope_se + SD + autocorrelation, data = variable.min.pub.wtemp)
summary(var_min_model_wtemp)

#### Nitrite + Nitrate ----

# Filter data
variable.min.pub.nitrogen <- variable.min.mag.pub %>%
  filter(Variable == "Nitrite + Nitrate, dissolved (mg/L)")

# Fit a linear model for 100% correct based on predictors
var_min_model_nitrogen <- lm(`Minimum Years (100% correct)` ~ slope + slope_se + SD + autocorrelation, data = variable.min.pub.nitrogen)
summary(var_min_model_nitrogen)

#### Total Suspended Solids ----

# Filter data
variable.min.pub.tss <- variable.min.mag.pub %>%
  filter(Variable == "Suspended Solids (mg/L)")

# Fit a linear model for 100% correct based on predictors
var_min_model_tss <- lm(`Minimum Years (100% correct)` ~ slope + slope_se + SD + autocorrelation, data = variable.min.pub.tss)
summary(var_min_model_tss)

# Functions for Variance Partitioning ----

# Extract variance explained from each model
get_variance_partition <- function(model) {
  anova_table <- car::Anova(model, type = "II")
  total_ss <- sum(anova_table$"Sum Sq")
  percent_explained <- (anova_table$"Sum Sq" / total_ss) * 100
  
  # Make it a one-row data frame with proper column names
  result <- as.data.frame(t(percent_explained))
  colnames(result) <- rownames(anova_table)  # Safely set names from ANOVA table rownames
  return(result)
}

# Summary Table: Variance Partitioning ----

# List of models
models <- list(
  DO = var_min_model_DO,
  DOS = var_min_model_DOS,
  WTemp = var_min_model_wtemp,
  Nitrogen = var_min_model_nitrogen,
  TSS = var_min_model_tss
)

# Extract variance partitioning and format
variance_partition_table <- purrr::map_dfr(models, get_variance_partition, .id = "Variable") %>%
  tibble::column_to_rownames("Variable") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Metric")

# Rename and format metrics
variance_partition_table$Metric <- fct_recode(variance_partition_table$Metric,
                                                       "Slope" = "slope",
                                                       "Slope SE" = "slope_se",
                                                       "Data SD" = "SD",
                                                       "Autocorrelation" = "autocorrelation"
)

# Round variance %s
variance_partition_table <- variance_partition_table %>%
  mutate(across(-Metric, ~ paste0(round(.x, 1), "%")))

# Extract minimum time values (from 80% and 100%) ---

# Extract minimum time (100% correct) and label threshold
variable.min.100 <- minTime_summary.100 %>%
  select(`Variable` = variable,
         `Minimum Years (100%)` = min_time_mag)

# Extract minimum time (80% correct) and label threshold
variable.min.80 <- minTime_summary.80 %>%
  select(`Variable` = variable,
         `Minimum Years (80%)` = min_time_mag)

# Combine both datasets into one long table
variable.min <- full_join(
  variable.min.80,
  variable.min.100
)

# Re-shape the dataset to put both thresholds under one column
min_time_table <- variable.min %>%
  pivot_longer(
    cols = c(`Minimum Years (80%)`, `Minimum Years (100%)`),
    names_to = "Metric",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = Variable, values_from = value) %>%
  mutate(across(-Metric, as.character))

# Ensure column order matches variance_partition_table
min_time_table <- min_time_table %>%
  select(
    Metric,
    DO = `Dissolved Oxygen_mg/l`,
    DOS = `Dissolved Oxygen Saturation_%`,
    WTemp = `Temperature Water_deg c`,
    Nitrogen = `Nitrogen, Nitrite (No2) + Nitrate (No3) As N Diss_mg/l`,
    TSS = `Solids, Suspended Total_mg/l`
  )

# Combine both tables
final_summary_table <- bind_rows(
  variance_partition_table,
  min_time_table
)

# Re-name columns
final_varPart_Time_table <- final_summary_table %>%
  rename(
    "Dissolved Oxygen (mg/L)" = DO,
    "Dissolved Oxygen Saturation (%)" = `DOS`,
    "Water Temperature (°C)" = WTemp,
    "Nitrite + Nitrate, dissolved (mg/L)" = Nitrogen,
    "Suspended Solids (mg/L)" = TSS
  )

# Coefficient of Variation ----

# Function to calculate coefficient of variation (CV) to compare across different units and scales
cv <- function(x) {
  (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100
}

# Run the CV analysis
Data.cv <- Data %>%
  mutate(cv = map_dbl(data, ~ cv(.x$response)))

# Rename and format metrics
Data.cv$variable_units <- fct_recode(Data.cv$variable_units,
                           "Dissolved Oxygen (mg/L)" = "Dissolved Oxygen_mg/l",
                           "Dissolved Oxygen Saturation (%)" = "Dissolved Oxygen Saturation_%",
                           "Nitrite + Nitrate, dissolved (mg/L)" = "Nitrogen, Nitrite (No2) + Nitrate (No3) As N Diss_mg/l",
                           "Suspended Solids (mg/L)" = "Solids, Suspended Total_mg/l",
                           "Water Temperature (°C)" = "Temperature Water_deg c"
)

# Round CVs
Data.cv.1 <- Data.cv %>%
  mutate(cv = paste0(round(cv, 1), "%"))

# Reshape the data from long to wide
Data.cv.1 <- Data.cv.1 %>%
  select(cv) %>%
  pivot_wider(names_from = variable_units, values_from = cv)

# Append CV to the final partitioning table
final_varPart_Time_table <- final_varPart_Time_table %>%
  bind_rows(Data.cv.1) %>%
  mutate(Metric = ifelse(is.na(Metric), "Coefficient of Variation", Metric))
