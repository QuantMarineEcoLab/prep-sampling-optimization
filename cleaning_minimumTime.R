# Clean minimum time data

# Load plyr before dplyr to avoid package conflict issues
pacman::p_load(
  plyr,
  dplyr,
  here,
  ggplot2
  )

# Load data
data <- read.csv(here("data/Great Bay Environmental Data/clean_imputedPMM_Great Bay_environmental data.csv"))

# Source functions needed for this script
source(here("sample_optimization/scripts/minimumTime_functions.R"))

# Clean data ----

# Clean site names
data$sitename <- clean_vec(data$sitename)

# Clean variable_units
data$variable_units <- clean_vec(data$variable_units)

# Re-name columns
variable <- data %>%
  rename("response" = "datavalue") %>%
  rename("site" = "sitename")

# Set numerical outputs in decimal places, not in scientific notations.
options(scipen = 10)

# Average each variable across sites
variable.1 <- variable %>%
  group_by(variable_units, year) %>%
  summarize(response = mean(response))

# Organize each site data into one place by creating a nested data frame
variable.nest <- variable.1 %>%
  group_by(variable_units) %>%
  tidyr::nest()

# Remove unwanted environmental variables
variable.test <- variable %>%
  filter(!variable_units %in% c("nitrogen_ammonia_as_n_dissolved_mg_l", 
                                "nitrogen_dissolved_total_mg_l", 
                                "ph_none", 
                                "phosphorus_as_p_total_mg_l", 
                                "phosphorus_orthophosphate_as_p_dissolved_mg_l", 
                                "turbidity_ntu"))

# Plot results to visualize trends
ggplot(variable.test, aes(x = year, y = response, color = site)) +
  geom_point() + 
  geom_line() +
  facet_wrap(~variable_units)