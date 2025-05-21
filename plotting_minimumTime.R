# Create figures for minimum time across Great Bay sites for each environmental variable
pacman::p_load(
  plyr,
  dplyr,
  here,
  forcats,
  formattable,
  htmltools,
  webshot
)

# Minimum Time Table (80% and 100% correct) ----

# Load data
source(here("sample_optimization/scripts/modeling_minimumTime_trendMagnitude.R"))

# Extract minimum time (100% correct) and label threshold
variable.min.100 <- minTime_summary.100 %>%
  select(variable, `100% correct` = min_time_mag)

# Extract minimum time (80% correct) and label threshold
variable.min.80 <- minTime_summary.80 %>%
  select(variable, `80% correct` = min_time_mag)

# Combine both datasets into one long table
variable.min <- full_join(
  variable.min.80,
  variable.min.100
)

# Re-name threshold columns
variable.min.2 <- variable.min %>%
  select(
    Variable = variable,
    "Minimum Years (80% correct)" = `80% correct`,
    "Minimum Years (100% correct)" = `100% correct`
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
min_time.table <- formattable(variable.min.pub,
  align = c("l", "l", "l"),
  list(area(col = c("Minimum Years (80% correct)", "Minimum Years (100% correct)"))
  ~ normalize_bar("pink", 0.2))
)

# Function to export minimum time table
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

# Export table
export_formattable(min_time.table, here("sample_optimization/figures/Minimum Effort, Maximum Insight Manuscript/minimumTime_variable_acrossSites.png"))


# Minimum Time Variance Partitioning Table ----

# Function to export minimum time table
export_formattable_plain <- function(df, file, width = "100%", height = NULL,
                                     background = "white", delay = 0.2) {
  library(formattable)
  library(htmltools)
  library(webshot)
  library(htmlwidgets)

  # Create formattable object with "Minimum Years" rows bold and colored
  f <- formattable(df, list(
    area(
      row = which(final_varPart_Time_table$Metric %in% c("Minimum Years (80%)", "Minimum Years (100%)")),
      col = 1:ncol(final_varPart_Time_table)
    ) ~
    formatter("span", style = ~ style(
      font.weight = "bold",
      display = "inline-block",
      width = "100%"
    ))
  ))

  # Convert to HTML widget
  w <- as.htmlwidget(f, width = width, height = height)

  # Custom CSS: wider column 1
  style_tag <- tags$style(HTML("
    .formattable_widget table th:nth-child(1),
    .formattable_widget table td:nth-child(1) {
      min-width: 200px !important;
      white-space: nowrap !important;
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
    zoom = 5
  )
}

# Export table
export_formattable_plain(
  final_varPart_Time_table,
  here("sample_optimization/figures/Minimum Effort, Maximum Insight Manuscript/minimumTime_variancePartitioning_trendMagnitude.png")
)
