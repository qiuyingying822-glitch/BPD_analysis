# Load necessary packages
setwd("/SparCC")
library(readxl)
library(MASS)
library(dplyr)
library(ggplot2)
library(reshape2)

# Set parameters
output_dir <- "model_results"
if (!dir.exists(output_dir)) dir.create(output_dir)
cat("=== BPD Severity Ordinal Logistic Regression Model ===\n")

# Read network topology metrics data
network_data <- read_excel("Network_topology_index.xlsx")
colnames(network_data)[1] <- "SampleID"
cat("Network topology data: ", nrow(network_data), "rows ¡Á", ncol(network_data), "columns\n")

# Read clinical metadata
metadata <- read.csv("metadata_mice.csv", sep = ",", row.names = 1)
metadata$SampleID <- row.names(metadata)
cat("Clinical metadata: ", nrow(metadata), "rows ¡Á", ncol(metadata), "columns\n")

# Data preprocessing

# Network metrics standardization
network_vars <- c(
  "Vertex.number", "Edge.number", "Diameter", "Average.path.length",
  "Average.nearest.neighbor.degree", "Betweenness.centralization",
  "Density", "Degree.centralization", "Degree.assortativity", "Transitivity"
)
str(network_data)

convert_to_numeric <- function(x) {
  if(is.character(x)) {
    # Remove possible spaces and special characters
    x_clean <- gsub("[^0-9.-]", "", x)
    # Convert to numeric
    result <- as.numeric(x_clean)
    # Check if conversion was successful
    if(all(is.na(result))) {
      cat("Warning: Column conversion failed, keeping original\n")
      return(x)
    } else {
      cat("Successfully converted to numeric\n")
      return(result)
    }
  } else {
    return(x)
  }
}

# Apply conversion (exclude Sample ID column)
columns_to_convert <- setdiff(names(network_data), c("SampleID", "group"))
for(col in columns_to_convert) {
  if(is.character(network_data[[col]])) {
    cat("Converting column:", col, "...")
    network_data[[col]] <- convert_to_numeric(network_data[[col]])
  }
}

str(network_data)

# Merge data
merged_data <- network_data %>%
  inner_join(metadata, by = c("SampleID" = "SampleID")) %>%
  filter(!is.na(Group))  # Remove missing BPD groups

cat("Merged data: ", nrow(merged_data), "samples\n")
str(merged_data)
merged_data <- merged_data[,-2]

# Variable renaming and transformation
df <- merged_data %>%
  mutate(
    # Convert BPD severity to ordered factor
    Group = ordered(Group, 
                    levels = c(0, 1, 2, 3),
                    labels = c("None", "Mild", "Moderate", "Severe")),
    
    # Categorical variable transformation
    Gender = factor(Gender, levels = c("Male", "Female"), labels = c("Male", "Female")),
    Antimicrobial_type = as.factor(Antimicrobial_type),
    Delivery_mode = factor(Delivery_mode, levels = c(0, 1), labels = c("C-section", "Vaginal")),
    Chorioamnionitis = factor(Chorioamnionitis, levels = c(0, 1)),
    ACS = factor(ACS, levels = c(0, 1, 2), ordered = TRUE),
    MSF = factor(MSF, levels = c(0, 1, 2, 3), ordered = TRUE),
    Intrauterine = factor(Intrauterine, levels = c(0, 1)),
    EUGR = factor(EUGR, levels = c(0, 1)),
    hospital_stay = as.numeric(hospital_stay)
  )

str(df)
print(table(merged_data$Gender, useNA = "always"))

# Standardize network metrics
for (var in network_vars) {
  if (var %in% names(df)) {
    # Check if numeric
    if (is.numeric(df[[var]])) {
      # Remove missing values
      valid_values <- df[[var]][!is.na(df[[var]])]
      if (length(valid_values) > 0 && sd(valid_values) > 0) {
        df[[paste0(var, "_std")]] <- scale(df[[var]])
        cat("Successfully standardized:", var, "\n")
      } 
    }
  }
}

# Define variable priority (based on clinical importance)
variable_priority <- c(
  # Network metrics (core independent variables)
  "Vertex.number_std", "Edge.number_std", "Density_std", 
  "Transitivity_std", "Degree.centralization_std",
  "Betweenness.centralization_std", "Average.path.length_std",
  
  # Most important clinical confounders
  "Gestational_age", "Birth_Weight", "Gender",
  
  # Other clinical factors
  "Delivery_mode", "ACS", "Chorioamnionitis",
  "MSF", "Intrauterine_distress", "EUGR", "hospital_stay"
)

# Forward stepwise variable selection function
forward_stepwise_selection <- function(df, outcome_var, candidate_vars) {
  current_vars <- character(0)
  stable_models <- list()
  
  for(var in candidate_vars) {
    if(!var %in% names(df)) {
      cat("Variable does not exist:", var, "\n")
      next
    }
    
    # Check if variable has sufficient variation
    if(is.numeric(df[[var]])) {
      if(sd(df[[var]], na.rm = TRUE) == 0) {
        cat("Variable has no variation, skipping:", var, "\n")
        next
      }
    }
    
    # Build current model formula
    if(length(current_vars) == 0) {
      formula_str <- paste(outcome_var, "~", var)
    } else {
      formula_str <- paste(outcome_var, "~", paste(c(current_vars, var), collapse = " + "))
    }
    
    formula_obj <- as.formula(formula_str)
    
    # Try to fit model
    model <- tryCatch({
      polr(formula_obj, data = df, Hess = TRUE, method = "logistic")
    }, error = function(e) {
      cat("Model convergence issue, skipping variable:", var, "Error:", e$message, "\n")
      return(NULL)
    })
    
    # Check Hessian matrix singularity
    if(!is.null(model)) {
      if(!any(is.na(coef(model))) && !any(diag(vcov(model)) < 0)) {
        current_vars <- c(current_vars, var)
        stable_models[[length(stable_models) + 1]] <- list(
          variables = current_vars,
          model = model,
          formula = formula_str
        )
        cat("Successfully added variable:", var, "\n")
      } else {
        cat("Numerical instability, skipping variable:", var, "\n")
      }
    }
  }
  
  return(stable_models)
}

# Execute forward stepwise selection
cat("3. Performing forward stepwise variable selection...\n")

candidate_vars <- variable_priority[variable_priority %in% names(df)]
stable_models <- forward_stepwise_selection(df, "Group", candidate_vars)

# Select final model
if(length(stable_models) > 0) {
  final_model <- stable_models[[length(stable_models)]]$model
  final_formula <- stable_models[[length(stable_models)]]$formula
  cat("Final model formula:", final_formula, "\n")
} else {
  stop("No stable model found")
}

# Model diagnostics and result extraction
cat("4. Model diagnostics and result extraction...\n")

# Calculate odds ratios and confidence intervals
model_summary <- summary(final_model)
coef_table <- coef(model_summary)
odds_ratios <- exp(coef_table[, "Value"])

cat("Attempting bootstrap method for confidence intervals...\n")
n_boot <- 1000
boot_coefs <- matrix(NA, nrow = n_boot, ncol = length(coef(final_model)))
colnames(boot_coefs) <- names(coef(final_model))

for(i in 1:n_boot) {
  tryCatch({
    # Bootstrap resampling
    boot_indices <- sample(1:nrow(df), replace = TRUE)
    boot_data <- df[boot_indices, ]
    boot_model <- polr(final_model$terms, data = boot_data, Hess = FALSE, method = "logistic")
    boot_coefs[i, ] <- coef(boot_model)
  }, error = function(e) {
    # Ignore single bootstrap failure
  })
}

# Calculate bootstrap confidence intervals
ci_boot <- apply(boot_coefs, 2, function(x) {
  quantile(na.omit(x), c(0.025, 0.975))
})
ci <- t(exp(ci_boot))
cat("Bootstrap confidence intervals calculated successfully\n")

# Exclude thresholds (cutpoints) from coef_table
# Threshold names usually contain "|" symbol
threshold_indices <- grep("\\|", rownames(coef_table))
coefficient_indices <- setdiff(1:nrow(coef_table), threshold_indices)

cat("Debug information:\n")
cat("Total rows in coef_table:", nrow(coef_table), "\n")
cat("Threshold rows:", length(threshold_indices), "\n")
cat("Coefficient rows:", length(coefficient_indices), "\n")
cat("Coefficient names:", rownames(coef_table)[coefficient_indices], "\n")

# Keep only coefficient part (exclude thresholds)
coef_table_coefs <- coef_table[coefficient_indices, ]

# Extract data based on actual column names
coefficient_values <- coef_table_coefs[, 1]  # Coefficient values
std_errors <- coef_table_coefs[, 2]          # Standard errors
t_values_coef <- coef_table_coefs[, 3]       # t-values (directly from model)

# Manually calculate p-values (based on t-test)
p_values <- 2 * (1 - pnorm(abs(t_values_coef)))

odds_ratios <- exp(coefficient_values)

# Create results table
results_table <- data.frame(
  Variable = rownames(coef_table_coefs),
  Coefficient = coefficient_values,
  Std_Error = std_errors,
  t_value = t_values_coef,
  OR = odds_ratios,
  CI_lower = ci[, 1],
  CI_upper = ci[, 2],
  p_value = p_values,
  Significance = ifelse(p_values < 0.001, "***",
                        ifelse(p_values < 0.01, "**",
                               ifelse(p_values < 0.05, "*", 
                                      ifelse(p_values < 0.1, ".", " ")))),
  stringsAsFactors = FALSE
)

# Save results
write.csv(results_table, file.path(output_dir, "ordinal_regression_results.csv"), row.names = FALSE)

cat("Final model results:\n")
print(results_table)