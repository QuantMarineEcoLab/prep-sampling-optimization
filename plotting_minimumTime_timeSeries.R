# Plot time series across Great Bay sites for each environmental variable

pacman::p_load(
  plyr,
  dplyr,
  here,
  purrr,
  forcats,
  tidyr,
  ggplot2,
  cowplot
)

# Load cleaned data ----
source(here("sample_optimization/scripts/cleaning_minimumTime.R"))

# Format variable names for plotting ----
variable.nest$variable_units_formatted <- fct_recode(variable.nest$variable_units,
  "Dissolved Oxygen (mg/L)" = "dissolved_oxygen_mg_l",
  "Dissolved Oxygen Saturation (%)" = "dissolved_oxygen_saturation_percent",
  "Ammonia, dissolved (mg/L)" = "nitrogen_ammonia_as_n_dissolved_mg_l",
  "Nitrogen, dissolved (mg/L)" = "nitrogen_dissolved_total_mg_l",
  "Nitrite + Nitrate, dissolved (mg/L)" = "nitrogen_nitrite_no2_nitrate_no3_as_n_diss_mg_l",
  "pH" = "ph_none",
  "Phosphorus (mg/L)" = "phosphorus_as_p_total_mg_l",
  "Phosphorus [orthophosphate], dissolved (mg/L)" = "phosphorus_orthophosphate_as_p_dissolved_mg_l",
  "Suspended Solids (mg/L)" = "solids_suspended_total_mg_l",
  "Water Temperature (°C)" = "temperature_water_deg_c",
  "Turbidity (ntu)" = "turbidity_ntu"
)

# Create time series plots with slope and p-value annotations ----
variable.nest <- variable.nest %>%
  mutate(
    model_stats = map(data, function(df) {
      model <- lm(response ~ year, data = df)
      summary_model <- summary(model)
      slope <- coef(model)[["year"]]
      p_value <- coef(summary_model)[2, 4]
      list(slope = slope, p_value = p_value)
    }),
    time_series = map2(data, model_stats, function(data, stats) {
      ggplot(data = data, aes(x = year, y = response)) +
        geom_point(shape = 21, color = "hotpink2", fill = "pink", size = 3.5) +
        geom_smooth(method = "lm", se = TRUE, alpha = 0.1, color = "black") +
        labs(y = variable_units_formatted, x = "Year") +
        annotate("text",
          x = mean(range(data$year, na.rm = TRUE)), y = Inf,
          label = paste0(
            "slope = ", round(stats$slope, 3),
            "; p = ", round(stats$p_value, 3)
          ),
          hjust = 0.5, vjust = 1.5, size = 4, color = "black", fontface = "bold"
        ) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          legend.position = "none"
        )
    })
  )

# Save all variable plots individually ----
for (i in 1:length(variable.nest$time_series)) {
  ggsave(
    filename = paste0(variable.nest$variable_units[i], ".jpeg"),
    plot = variable.nest$time_series[[i]],
    path = here("sample_optimization/figures/minimum_time/time_series/across_sites")
  )
}

# Create publication-ready subset of plots ----

# Remove Rt. 108 Bridge Mill Pond
variable.108 <- variable %>%
  filter(site != "rte_108_bridge_mill_pond")

# Average each variable across sites by year
variable.108 <- variable.108 %>%
  group_by(variable_units, year) %>%
  summarize(response = mean(response), .groups = "drop")

# Nest the data
variable.nest <- variable.108 %>%
  group_by(variable_units) %>%
  nest()

# Re-format variable names
variable.nest$variable_units_formatted <- fct_recode(variable.nest$variable_units,
  "Dissolved Oxygen (mg/L)" = "dissolved_oxygen_mg_l",
  "Dissolved Oxygen Saturation (%)" = "dissolved_oxygen_saturation_percent",
  "Ammonia, dissolved (mg/L)" = "nitrogen_ammonia_as_n_dissolved_mg_l",
  "Nitrogen, dissolved (mg/L)" = "nitrogen_dissolved_total_mg_l",
  "Nitrite + Nitrate, dissolved (mg/L)" = "nitrogen_nitrite_no2_nitrate_no3_as_n_diss_mg_l",
  "pH" = "ph_none",
  "Phosphorus (mg/L)" = "phosphorus_as_p_total_mg_l",
  "Phosphorus [orthophosphate], dissolved (mg/L)" = "phosphorus_orthophosphate_as_p_dissolved_mg_l",
  "Suspended Solids (mg/L)" = "solids_suspended_total_mg_l",
  "Water Temperature (°C)" = "temperature_water_deg_c",
  "Turbidity (ntu)" = "turbidity_ntu"
)

# Create time series plots with slope and p-value for publication subset ----
variable.nest <- variable.nest %>%
  mutate(
    model_stats = map(data, function(df) {
      model <- lm(response ~ year, data = df)
      summary_model <- summary(model)
      slope <- coef(model)[["year"]]
      p_value <- coef(summary_model)[2, 4]
      list(slope = slope, p_value = p_value)
    }),
    time_series = map2(data, model_stats, function(data, stats) {
      ggplot(data = data, aes(x = year, y = response)) +
        geom_point(shape = 21, color = "hotpink2", fill = "pink", size = 3.5) +
        geom_smooth(method = "lm", se = TRUE, alpha = 0.1, color = "black") +
        labs(y = variable_units_formatted, x = "Year") +
        annotate("text",
          x = mean(range(data$year, na.rm = TRUE)), y = Inf,
          label = paste0(
            "slope = ", round(stats$slope, 3),
            "; p = ", round(stats$p_value, 3)
          ),
          hjust = 0.5, vjust = 1.5, size = 4, color = "black", fontface = "bold"
        ) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          legend.position = "none"
        )
    })
  )

# Select key variables for the manuscript ----
DO_mgL <- variable.nest %>% filter(variable_units == "dissolved_oxygen_mg_l")
DO_sat <- variable.nest %>% filter(variable_units == "dissolved_oxygen_saturation_percent")
nitrogen <- variable.nest %>% filter(variable_units == "nitrogen_nitrite_no2_nitrate_no3_as_n_diss_mg_l")
solids <- variable.nest %>% filter(variable_units == "solids_suspended_total_mg_l")
w_temp <- variable.nest %>% filter(variable_units == "temperature_water_deg_c")

# Arrange selected plots in a grid layout ----
grid.plot <- plot_grid(
  DO_mgL$time_series[[1]],
  DO_sat$time_series[[1]],
  nitrogen$time_series[[1]],
  solids$time_series[[1]],
  w_temp$time_series[[1]],
  align = "hv",
  labels = "AUTO",
  ncol = 2
)

# Save final multi-panel figure ----
ggsave(
  filename = "manuscript-variable-time-series.jpeg",
  plot = grid.plot,
  path = here("sample_optimization/figures/Minimum Effort, Maximum Insight Manuscript"),
  width = 9,
  height = 12
)
