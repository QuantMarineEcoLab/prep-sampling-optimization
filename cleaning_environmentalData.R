# Clean Great Bay environmental data

# Collate Great Bay environmental data ----
# into one main data frame and select appropriate variables for analyses

# Load plyr before dplyr to avoid package conflict issues
pacman::p_load(
  plyr,
  dplyr,
  tidyverse,
  here,
  Rmisc,
  mice,
  ggplot2,
  cowplot
)

# Create a function that reads a CSV file and adds a new column 'filename' to the resulting data frame
# A good way to keep track of the filename associated with each row in the data frame.
read_file <- function(flnm) {
  readr::read_csv(flnm) %>%
    mutate(filename = flnm)
}

# Create a list of files, excluding those that start with 'clean'
gb_env_list <- list.files(
  path = here("data/Great Bay Environmental Data/"),
  pattern = "\\.csv$",      # Only CSV files
  full.names = TRUE         # Return full paths
) %>%
  # Filter out files that start with "clean"
  .[!grepl("^clean", basename(.))] %>%
  # Read each file
  map(~ read_file(.))

# ^^ Comments for above code
# list.files: Get a list of file names in the current directory that match the specified pattern.

# The pattern "\\.csv$" ensures that only files with a CSV extension are included. The full.names = TRUE argument ensures that the full path to each file is returned.

# Pattern breakdown:
# \\: The double backslash (\\) is an escape sequence in regular expressions. In R, you need to escape the dot (.) with a backslash to match a literal dot. So, \\. matches the dot in ".csv".
# .csv: Matches the literal characters ".csv".
# $: The dollar sign ($) asserts the position at the end of the string. In this context, it ensures that the pattern ".csv" is matched only at the end of the file name.

# map(~read_file(.)): uses the map function to apply the read_file function to each element of the list. The ~read_file(.) is a shorthand for an anonymous function that takes a single argument and applies read_file to it.

# Combine 'variable' and 'units' into 1 column and remove the original columns
gb_env_clean <- gb_env_list %>%
  map_dfr(~ unite(., variable_units, variable, units, remove = T)) %>% # Combine 'variable' and 'units' into 1 column
  na.omit() # Remove any rows that contain missing values

# Quickly look at environmental variables of interest present in the data
env_variable_list <- gb_env_clean %>%
  select(., variable_units) %>%
  unique(.) %>% # Remove redundant variable names
  mutate(oxygen = case_when( # Create a new column based on criteria
    str_detect(variable_units, "(?i)oxygen") ~ variable_units, # Identify variable_units that contain the pre-defined word. Note (?!) ensures the search is case-insensitive.
    TRUE ~ NA_character_
  )) %>% # Otherwise, set it as NA
  mutate(temperature = case_when(
    str_detect(variable_units, "(?i)temperature") ~ variable_units,
    TRUE ~ NA_character_
  )) %>%
  mutate(nitrogen = case_when(
    str_detect(variable_units, "(?i)nitrogen") ~ variable_units,
    TRUE ~ NA_character_
  )) %>%
  mutate(phosphorus = case_when(
    str_detect(variable_units, "(?i)phosphorus") ~ variable_units,
    TRUE ~ NA_character_
  )) %>%
  mutate(turbidity = case_when(
    str_detect(variable_units, "(?i)turbidity") ~ variable_units,
    TRUE ~ NA_character_
  )) %>%
  mutate(carbon = case_when(
    str_detect(variable_units, "(?i)carbon") ~ variable_units,
    TRUE ~ NA_character_
  )) %>%
  mutate(solid = case_when(
    str_detect(variable_units, "(?i)solid") ~ variable_units,
    TRUE ~ NA_character_
  )) %>%
  mutate(pH = case_when(
    str_detect(variable_units, "(?i)ph_none|ph total_none") ~ variable_units,
    TRUE ~ NA_character_
  )) %>%
  mutate(chlorophyll = case_when(
    str_detect(variable_units, "(?i)chlorophyll") ~ variable_units,
    TRUE ~ NA_character_
  )) %>%
  mutate(cdom_fdom = case_when(
    str_detect(variable_units, "(?i)cdom|fdom") ~ variable_units,
    TRUE ~ NA_character_
  ))

# Remove 'variable_units' and remove any redundant variable names
env_variable_list_condensed <- env_variable_list %>%
  select(-variable_units) %>%
  filter_all(any_vars(!is.na(.)))

# See summary data by site
gb_env_summary_site <- gb_env_clean %>%
  Rmisc::summarySE(measurevar = "datavalue", groupvars = c("sitecode", "sitename", "variable_units"))

# See summary data
gb_env_summary <- gb_env_clean %>%
  Rmisc::summarySE(measurevar = "datavalue", groupvars = c("variable_units"))

# Choose variables of interest
# Barely any Chl-A, CDOM, fDOM data -- remove these from the environmental data
gb_env_core_data <- gb_env_clean %>%
  filter(variable_units %in% c("Temperature Water_deg c", "Dissolved Oxygen Saturation_%", "Dissolved Oxygen_mg/l", "Turbidity_ntu", "Solids, Suspended Total_mg/l", "pH_none", "Carbon, Organic Dissolved_mg/l", "Phosphorus As P Total_mg/l", "Phosphorus, Orthophosphate As P Dissolved_mg/l", "Nitrogen Total_mg/l", "Nitrogen, Organic Dissolved_mg/l", "Nitrogen, Nitrite (No2) + Nitrate (No3) As N Diss_mg/l", "Nitrogen, Dissolved Total_mg/l", "Nitrogen, Ammonia As N Dissolved_mg/l")) %>%
  filter(dataqualitycode == "Good")

# Add date columns
gb_env_core_data$year <- as.numeric(format(gb_env_core_data$valuedatetime, "%Y"))
gb_env_core_data$month <- as.numeric(format(gb_env_core_data$valuedatetime, "%m"))
gb_env_core_data$day <- as.numeric(format(gb_env_core_data$valuedatetime, "%d"))

# Calculates the sampling frequency by calculating the time difference between each sampling event
gb_env_core_data_timediff <- gb_env_core_data %>%
  group_by(variable_units, sitename) %>%
  arrange(valuedatetime) %>%
  mutate(time_diff = c(NA, difftime(valuedatetime[-1], valuedatetime[-n()], units = "min"))) # Creates a vector of time differences. NA is added as the first element because the time difference for the first row is not defined. 'difftime' calculates the time difference in minutes between each row and the previous row within each group. [-1] and [-n()] exclude the first and last elements, respectively, to align the differences correctly.


# Run sampling summary stats to determine which variables to narrow down once more.
# Total samples = NA should be 0
sampling_summary <- gb_env_core_data_timediff %>%
  group_by(variable_units, sitename) %>%
  summarize(
    start_year = min(year),
    end_year = max(year),
    total_samples = n(),
    total_samples_NA = sum(is.na(datavalue)),
    sampling_freq_day = ((mean(time_diff, na.rm = T) / 60) / 24)
  )


# Yearly data ----

## Data preparation ----

# Average data by year
yearly_var <- gb_env_core_data %>%
  dplyr::group_by(sitename, year, variable_units) %>%
  summarize(datavalue = mean(datavalue))

# Subset only those sites with >= 15 yearly data points
yearly_var.1 <- yearly_var %>%
  group_by(sitename, variable_units) %>%
  filter(length(sitename) >= 15)

# Run summary statistics on sites with yearly data
yearly_var_summary <- summarySE(data = yearly_var.1, measurevar = "datavalue", groupvars = c("sitename", "variable_units"), na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)

# Trim data range - start at 2000
yearly_var_range <- subset(yearly_var.1, year >= 2000)

# Add sampling range (start and end years)
yearly_var_range.1 <- yearly_var_range %>%
  group_by(sitename, variable_units) %>%
  summarize(start_year = min(year), end_year = max(year))

# Subset only sites with at least one yearly data point between 2019 and 2023
yearly_var_range.2 <- subset(yearly_var_range.1, end_year >= 2019)

# Drop unused levels
yearly_var_range.2 <- droplevels(yearly_var_range.2)

# Add back in environmental data using only the chosen sites
yearly_var_range.3 <- left_join(yearly_var_range.2, yearly_var_range)

# Make year numeric
yearly_var_range.3$year <- as.numeric(yearly_var_range.3$year)

## Insert NAs for unsampled data ----

# Create a range of years
year <- c(2000:2022)

# Extract a list of site names
site <- yearly_var_range.2$sitename

# Extract a list of variables
var <- yearly_var_range.2$variable_units

# Repeat site rows for x years
site_rep <- data.frame(rep(site, each = length(year)))

# Repeat variable rows for x years
var_rep <- data.frame(rep(var, each = length(year)))

# Repeat year rows for x sites
year_rep <- data.frame(rep(year, times = length(site)))

# Combine site and year vectors into a data frame
site_var_year <- data.frame(site_rep, var_rep, year_rep)

# Re-name columns
colnames(site_var_year)[1] <- "sitename"
colnames(site_var_year)[2] <- "variable_units"
colnames(site_var_year)[3] <- "year"

# Merge inserted NAs with the original data
site_var_year.1 <- full_join(site_var_year, yearly_var_range.3)

# See how many years are missing for each site
site_var_year.sum <- site_var_year.1 %>%
  group_by(sitename, variable_units) %>%
  summarize(num_na = sum(is.na(datavalue))) %>%
  mutate(percent_na = (num_na / 22) * 100)

# Filter out sites/variable with more than 25% missing data
site_var_year.sum_25 <- site_var_year.sum %>%
  filter(percent_na <= 25)

# Merge to the main data frame to filter out sites/variables > 25% missing data
site_var_year.2 <- left_join(site_var_year.sum_25, site_var_year.1) %>% droplevels()

## Impute missing data ----

# Observe patterns for missing data
md.pattern(site_var_year.2)

# Function to apply mice imputation within each group with a specified method
impute_group <- function(data, method = "pmm") {
  imputed_data <- mice(data, method = method, m = 5, printFlag = FALSE)  # Apply the specified method
  complete(imputed_data)
}

# Wrapper function to apply imputation across groups
impute_by_site_variable <- function(data, method = "pmm") {
  data %>%
    group_by(sitename, variable_units) %>%
    group_split() %>%
    lapply(impute_group, method = method) %>%
    bind_rows()
}

# Example usage - Specify the method for imputation
mice_imputed_pmm <- impute_by_site_variable(site_var_year.2, method = "pmm")
mice_imputed_cart <- impute_by_site_variable(site_var_year.2, method = "cart")

# Plot different types of imputed data to see which one matches the original distribution the best
# PMM-imputed is the best fit
h1 <- ggplot(site_var_year.2, aes(x = datavalue)) +
  geom_histogram(fill = "#ad1538", color = "#000000", position = "identity") +
  ggtitle("Original distribution") +
  theme_classic()
h2 <- ggplot(mice_imputed_pmm, aes(x = datavalue)) +
  geom_histogram(fill = "#15ad4f", color = "#000000", position = "identity") +
  ggtitle("PMM-imputed distribution") +
  theme_classic()
h3 <- ggplot(mice_imputed_cart, aes(x = datavalue)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("CART-imputed distribution") +
  theme_classic()

plot_grid(h1, h2, h3, nrow = 2, ncol = 2)

# Clean up data and use PMM-imputed values
mice_imputed.1 <- mice_imputed_pmm %>%
  select(sitename, variable_units, year, datavalue)

# Export clean temperature data
write.csv(mice_imputed.1, "./data/Great Bay Environmental Data/clean_imputedPMM_Great Bay_environmental data.csv", row.names = F)
