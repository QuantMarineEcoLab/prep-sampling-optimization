# Functions to calculate the minimum time required to detect trends

# Clean column names (lowercase, snake case, etc.)
clean_vec <- function(x, refactor = FALSE) {
  require(magrittr, quietly = TRUE)
  
  if (!(is.character(x) || is.factor(x))) {
    return(x)
  } else {
    old_names <- as.character(x)
    
    new_names <- old_names %>%
      gsub("'", "", .) %>%
      gsub("\"", "", .) %>%
      gsub("%", "percent", .) %>%
      gsub("^[ ]+", "", .) %>%
      make.names(.) %>%
      gsub("[.]+", "_", .) %>%
      gsub("[_]+", "_", .) %>%
      tolower(.) %>%
      gsub("_$", "", .)
    
    return(new_names)
  }
}

# Standardize response values
standardize <- function(data) {
  
  # Ensure year is numeric
  data$year <- as.numeric(as.character(data$year))
  
  # Standardize the response column
  data$stand.response <- (data$response - mean(data$response, na.rm = TRUE)) /
    sd(data$response, na.rm = TRUE)
  
  # Reorganize the columns based on what exists
  if ("site" %in% names(data)) {
    data <- data[, c("site", "year", "response", "stand.response")]
  } else {
    data <- data[, c("year", "response", "stand.response")]
  }
  
  return(data)
}

# Fit linear regression without site data
linefit_no_site <- function(data) {
  
  # Fit a simple linear regression model to standardized response over time
  model <- lm(stand.response ~ year, data = na.omit(data))
  
  # Extract model summary
  t.model <- coef(summary(model))
  slope_pval <- summary(model)$coefficients["year", "Pr(>|t|)"]
  
  # Package key statistics into a clean data frame
  output <- data.frame(
    start_year = min(data$year, na.rm = TRUE),              # First year in the window
    N_data = nrow(data),                                     # Number of observations
    N_years = length(unique(data$year)),                     # Number of unique years
    slope = t.model["year", "Estimate"],                     # Slope estimate
    slope_se = t.model["year", "Std. Error"],                # Slope standard error
    p_value = slope_pval,                                    # p-value for slope
    intercept = t.model["(Intercept)", "Estimate"],          # Intercept estimate
    intercept_se = t.model["(Intercept)", "Std. Error"],     # Intercept standard error
    intercept_p_value = t.model["(Intercept)", "Pr(>|t|)"]   # p-value for intercept
  )
  
  return(output)
}

# Fit linear regression with site data
linefit_with_site <- function(data) {
  
  # Load required libraries
  library(lme4)
  library(lmerTest)
  
  data <- na.omit(data)                     # Clean data
  n_sites <- length(unique(data$site))     # Count number of unique sites
  
  if (n_sites < 3) {
    # Fallback to simple linear model if not enough sites for random effects
    return(linefit_no_site(data))
  }
  
  # Fit a linear mixed-effects model with random intercept for site
  model <- lmer(stand.response ~ year + (1 | site), data = data)
  t.model <- coef(summary(model))
  slope_pval <- t.model["year", "Pr(>|t|)"]
  
  # Package fixed-effect statistics into a clean data frame
  output <- data.frame(
    start_year = min(data$year, na.rm = TRUE),
    N_data = nrow(data),
    N_years = length(unique(data$year)),
    slope = t.model["year", "Estimate"],
    slope_se = t.model["year", "Std. Error"],
    p_value = slope_pval,
    intercept = t.model["(Intercept)", "Estimate"],
    intercept_se = t.model["(Intercept)", "Std. Error"],
    intercept_p_value = t.model["(Intercept)", "Pr(>|t|)"]
  )
  
  return(output)
}

# Function to determine which linefit model to use
linefit_wrapper <- function(data) {
  
  # Decide which model to use based on presence of 'site' column
  if ("site" %in% names(data)) {
    result <- linefit_with_site(data)
    attr(result, "model_type") <- "lmer"   # Attach model type info
    return(result)
  } else {
    result <- linefit_no_site(data)
    attr(result, "model_type") <- "lm"
    return(result)
  }
}

# Subset data according to the window size
breakup <- function(data, window, verbose = FALSE) {
  remaining <- data  # Working copy of data
  
  # Initialize output data frame to store stats for each window
  output <- data.frame(
    start_year = numeric(0),
    N_data = numeric(0),
    N_years = numeric(0),
    slope = numeric(0),
    slope_se = numeric(0),
    p_value = numeric(0),
    intercept = numeric(0),
    intercept_se = numeric(0),
    intercept_p_value = numeric(0),
    model_type = character(0)  # Track which model was used (lm or lmer)
  )
  
  # Count total unique years in the dataset
  numyears <- length(unique(data$year))
  
  # Loop through the dataset, shifting the window one year at a time
  while (numyears > (window - 1)) {
    chunk <- subset(remaining, year < (min(year) + window))  # Subset a chunk
    
    # Only process if full window of unique years exists
    if (window == length(unique(chunk$year))) {
      out <- linefit_wrapper(chunk)  # Fit appropriate model
      out$model_type <- attr(out, "model_type")  # Add model type info
      
      # Optional logging to console
      if (verbose) {
        message("Fitting model from ", min(chunk$year), " to ", max(chunk$year),
                " | Model used: ", out$model_type)
      }
      
      output <- rbind(output, out)
    }
    
    # Move window forward by removing the earliest year
    remaining <- subset(remaining, year > min(year))
    numyears <- length(unique(remaining$year))
  }
  
  return(output)
}

# Subset data across all window sizes
multiple_breakups <- function(data, verbose = FALSE) {
  data1 <- standardize(data)       # Standardize response across all data
  count <- length(data1$year)      # Total number of data points
  output <- data.frame()           # Empty frame to collect all breakups
  
  # Loop through all possible window sizes (min 3 years)
  for (i in 3:count) {
    outeach <- breakup(data1, i, verbose = verbose)
    output <- rbind(output, outeach)
  }
  
  return(output)
}

# Calculate minimum time using trend magnitude metric
min_time_magnitude <- function(data, min_percent = 95, error_multiplier = 1) {
  
  # Step 1: Run all breakups for various window sizes
  test <- multiple_breakups(data)
  
  # Step 2: Use the slope of the longest window as the "true" slope
  count <- nrow(test)
  true_slope <- test[count, "slope"]
  true_error <- test[count, "slope_se"] * sqrt(test[count, "N_data"]) * error_multiplier
  
  # Step 3: Define the accepted slope range
  max_true <- true_slope + true_error
  min_true <- true_slope - true_error
  
  # Step 4: Get list of all window sizes
  windows <- unique(test$N_years)
  stability <- max(windows)
  prop.vec <- c()
  
  # Step 5: Loop through each window size and calculate % within true range
  for (i in 1:length(windows)) {
    window_length <- windows[i]
    test_subset <- test[test$N_years == window_length, ]
    number_of_windows <- nrow(test_subset)
    
    # Check if slope falls within the accepted range
    correct_subset <- test_subset[
      (test_subset$slope < max_true) & (test_subset$slope > min_true), ]
    number_of_correct <- nrow(correct_subset)
    
    percentage_correct <- 100 * number_of_correct / number_of_windows
    prop.vec <- c(prop.vec, percentage_correct)
    
    if (percentage_correct >= min_percent && window_length < stability) {
      stability <- window_length
    }
  }
  
  # Step 6: Return results as dataframe
  df <- data.frame(
    window_length = windows,
    proportion_correct = prop.vec,
    stability_time = stability
  )
  
  return(df)
}
