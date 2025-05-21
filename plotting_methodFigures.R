# Method figures for the sample optimization paper 

# Load libraries
pacman::p_load(
  plyr,
  dplyr,
  here,
  ggplot2,
  cowplot
)

# Simulate declining variable
set.seed(42)
years <- 1:36
variable <- 50 - 0.3 * years + rnorm(length(years), 0, 1)
data <- data.frame(year = years, variable = pmax(variable, 20))

# Panel A with ±1 SD ribbon
main_model <- lm(variable ~ year, data = data)
predictions <- predict(main_model, newdata = data, se.fit = TRUE)

# Add predicted values and ±1 SD to data frame
data$fit <- predictions$fit
data$sd_upper <- data$fit + sd(data$variable)
data$sd_lower <- data$fit - sd(data$variable)

# Plot
main_plot <- ggplot(data, aes(x = year, y = variable)) +
  geom_ribbon(aes(ymin = sd_lower, ymax = sd_upper), fill = "red", alpha = 0.05) +
  geom_point(size = 2) +
  geom_line(aes(y = fit), color = "red", linewidth = 1) +
  labs(x = "Year", y = "Variable") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12)
  )

# Function for subsample plots (Panels B–D)
make_subsample_plot <- function(data, window_size) {
  n <- nrow(data)
  start_indices <- 1:(n - window_size + 1)
  base_plot <- ggplot(data, aes(x = year, y = variable)) + geom_point()
  for (i in start_indices) {
    sub_data <- data[i:(i + window_size - 1), ]
    base_plot <- base_plot +
      geom_smooth(data = sub_data, method = "lm", se = FALSE, color = "red", size = 0.6)
  }
  base_plot +
    labs(x = "Year", y = "Variable") +
    theme_minimal() +
    theme(axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 12))
}

# Panels B–D
plot_5yr <- make_subsample_plot(data, 5)
plot_10yr <- make_subsample_plot(data, 10)
plot_15yr <- make_subsample_plot(data, 15)

# Panel E: Subsample slope matching
window_length <- 3:22
percent_detected <- c(40, 55, 70, 80, 90, 95, rep(100, length(window_length) - 6))
first_80_year <- window_length[which(percent_detected >= 80)[1]]
first_100_year <- window_length[which(percent_detected >= 100)[1]]
df <- data.frame(Window = window_length, Percent = percent_detected)

plot_E <- ggplot(df, aes(x = Window, y = Percent)) +
  geom_line(color = "black", linewidth = 1) +
  geom_point(color = "black", size = 2) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(xintercept = first_80_year, linetype = "dotted", color = "blue", linewidth = 1) +
  labs(
    x = "Subsample Length (Year)",
    y = "% Matched Subsample Slopes"
  ) +
  ylim(30, 105) +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12))

# Assemble left column: A (main)
left_column <- plot_grid(main_plot, ncol = 1, labels = c("A"), label_size = 14)

# Assemble right column: B, C, D
middle_column <- plot_grid(plot_5yr, plot_10yr, plot_15yr, ncol = 1, labels = c("B", "C", "D"), label_size = 14)

# Assemble right column: E (slope match)
right_column <- plot_grid(plot_E, ncol = 1, labels = c("E"), label_size = 14)

# Final combined figure
final_plot <- plot_grid(left_column, middle_column, right_column, ncol = 3, rel_widths = c(1.3, 1, 1.3))
print(final_plot)

# Save the plot
ggsave(
  filename = paste0("method-figures.jpeg"),
  plot = final_plot,
  path = here("sample_optimization/figures/Minimum Effort, Maximum Insight Manuscript"),
  width = 10,
  height = 6
)
