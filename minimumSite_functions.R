# Functions to calculate the minimum # of sites required to detect trends

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
  
  library(lme4)
  library(lmerTest)
  
  # Remove missing values
  data <- na.omit(data)
  
  # Count unique sites
  n_sites <- length(unique(data$site))
  
  # If there are fewer than 3 sites, fallback to simple linear regression
  if (n_sites < 3) {
    return(linefit_no_site(data))
  }
  
  # Fit a mixed-effects model with random intercepts for site
  model <- lmer(stand.response ~ year + (1 | site), data = data)
  t.model <- coef(summary(model))
  slope_pval <- t.model["year", "Pr(>|t|)"]
  
  # Return fixed effect estimates in a data frame
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

# Subset data with every combination of sites according to the window size
breakup_systematic <- function(data, window) {
  
  # Step 1: Initialize an empty output dataframe
  output <- data.frame(
    combo_num = numeric(0),
    N_data = numeric(0),
    N_years = numeric(0),
    N_sites = numeric(0),
    perc_site_removed = numeric(0),
    site_removed = character(0),
    year_slope = numeric(0),
    slope_sd = numeric(0),
    slope_se = numeric(0),
    slope_pval = numeric(0),
    intercept = numeric(0)
  )
  
  # Step 2: Calculate how many sites to remove for this window
  numsites <- length(unique(data$site))
  perc_site_removed <- round(numsites * (window / 100))
  
  # Step 3: Generate all combinations of sites to remove
  site_combinations <- combn(unique(data$site), perc_site_removed, simplify = FALSE)
  
  # Step 4: Loop through each site combination
  for (iterations in 1:length(site_combinations)) {
    
    chunk <- site_combinations[[iterations]]
    
    # Step 4a: Remove selected sites
    chunk.1 <- subset(data, !site %in% chunk)
    
    # Step 4b: Fit model to remaining data
    out <- linefit_with_site(chunk.1)
    
    # Step 4c: Add site-removal metadata
    out$site_removed <- paste(chunk, collapse = ", ")
    out$combo_num <- iterations
    out$perc_site_removed <- perc_site_removed
    out$N_sites <- length(unique(chunk.1$site))
    
    # Step 4d: Convert SE to SD
    out$slope_sd <- out$slope_se * sqrt(out$N_data)
    
    # Step 4e: Append to output
    output <- rbind(output, out)
  }
  
  # Step 5: Select final columns to return
  output <- dplyr::select(output,
                          combo_num, N_data, N_years, N_sites, perc_site_removed, site_removed,
                          year_slope = slope, slope_sd, slope_se, slope_pval = p_value, intercept
  )
  
  return(output)
}

# Subset data across all window sizes
multiple_breakups <- function(data) {
  
  # Step 1: Standardize the data
  data1 <- standardize(data)
  
  # Step 2: Initialize an empty dataframe
  output <- data.frame(
    combo_num = numeric(0),
    N_data = numeric(0),
    N_years = numeric(0),
    N_sites = character(0),
    perc_site_removed = numeric(0),
    site_removed = character(0),
    year_slope = numeric(0),
    slope_sd = numeric(0),
    slope_se = numeric(0),
    slope_pval = numeric(0),
    intercept = numeric(0)
  )
  
  # Step 3: Loop over different % site removal windows
  for (i in seq(1, 100, by = 10)) {
    outeach <- breakup_systematic(data1, i)
    output <- rbind(output, outeach)
  }
  
  return(output)
}

# Calculate minimum site using trend magnitude metric
min_site_magnitude <- function(data, min_percent = 95, error_multiplier = 1) {
  
  # Step 1: Run multiple site-removal breakups
  test <- multiple_breakups(data)
  
  # Step 2: Extract the 'true' slope and error (no sites removed)
  true_slope <- unique(subset(test, perc_site_removed == 0)$year_slope)
  true_error <- unique(subset(test, perc_site_removed == 0)$slope_se) * error_multiplier
  
  # Step 3: Define acceptable slope range
  max_true <- true_slope + true_error
  min_true <- true_slope - true_error
  
  # Step 4: Get all unique N_sites values
  windows <- unique(test$N_sites)
  stability <- max(windows)
  prop.vec <- c()
  
  # Step 5: Loop through each N_sites value
  for (i in 1:length(windows)) {
    num_site <- windows[i]
    test_subset <- test[test$N_sites == num_site, ]
    number_of_windows <- nrow(test_subset)
    
    # Step 5a: Count how many runs have a slope within the "true" range
    correct_subset <- test_subset[
      test_subset$year_slope < max_true &
        test_subset$year_slope > min_true, ]
    number_of_correct <- nrow(correct_subset)
    
    # Step 5b: Calculate % of runs considered "stable"
    percentage_correct <- 100 * number_of_correct / number_of_windows
    prop.vec <- c(prop.vec, percentage_correct)
    
    # Step 5c: Update stability threshold if criteria met
    if (percentage_correct >= min_percent && num_site < stability) {
      stability <- num_site
    }
  }
  
  # Step 6: Return summary dataframe
  df <- data.frame(
    num_site = windows,
    proportion_correct = prop.vec,
    minimum_site = stability
  )
  
  return(df)
}

# Identify which site provides the most/least unique information and minimum site combinations
usable_sites <- function(data, min_percent = 95, error_multiplier = 1) {
  
  # Step 1: Run all site-removal combinations and model fits
  test <- multiple_breakups(data)
  
  # Step 2: Get the "true" slope from the full dataset (no site removal)
  true_slope <- unique(subset(test, perc_site_removed == 0)$year_slope)
  
  # Step 3: Calculate absolute difference from true slope for each combination
  test$slope_diff <- abs(test$year_slope - true_slope)
  
  # Step 4: Determine the minimum number of sites required for stability
  site_stability <- min_site_magnitude(data, min_percent = min_percent, error_multiplier = error_multiplier)
  min_sites_required <- min(site_stability$minimum_site, na.rm = TRUE)
  
  # Step 5: Filter to combinations that used exactly the minimum number of sites
  eligible_combos <- test[test$N_sites == min_sites_required, ]
  
  # Step 6: Identify the best matching site combo (smallest slope diff)
  best_match_combo <- eligible_combos[which.min(eligible_combos$slope_diff), , drop = FALSE]
  
  # Step 7: Find the largest and smallest slope diffs for one-site removal only
  one_site_removed <- test[test$perc_site_removed == 1, ]
  
  largest_diff_one_site <- one_site_removed[which.max(one_site_removed$slope_diff), , drop = FALSE]
  smallest_diff_one_site <- one_site_removed[which.min(one_site_removed$slope_diff), , drop = FALSE]
  
  # Step 8: Generate full sorted list of slope differences
  all_differences <- test %>%
    dplyr::select(site_removed, N_sites, year_slope, slope_diff) %>%
    dplyr::arrange(desc(slope_diff))
  
  # Step 9: Return all results in a nested tibble
  result <- tibble::tibble(
    minimum_sites_required = min_sites_required,
    smallest_difference_min_site = list(best_match_combo),
    largest_difference_one_site_removed = list(largest_diff_one_site),
    smallest_difference_one_site_removed = list(smallest_diff_one_site),
    all_slope_differences = list(all_differences)
  )
  
  return(result)
}
