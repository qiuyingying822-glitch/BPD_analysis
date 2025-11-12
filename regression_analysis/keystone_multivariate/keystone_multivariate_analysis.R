#!/usr/bin/env Rscript
#
# Keystone Taxa Multivariate Analysis
# Association between keystone ASVs and BPD severity with clinical adjustments
#
# This script performs multivariate ordinal regression to assess the relationship
# between keystone microbial taxa and BPD severity, adjusting for clinical covariates.

cat("=============================================\n")
cat("Keystone Taxa Multivariate Analysis\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================\n")

# Load required packages
cat("1. Loading required packages...\n")
required_packages <- c("missForest", "randomForest", "ordinal", "dplyr", 
                       "ggplot2", "readr", "tidyr", "MASS", "car", "readxl",
                       "patchwork", "tibble")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
cat("✅ All packages loaded successfully\n")

# Data preprocessing and imputation
cat("2. Data preprocessing and imputation...\n")

# Read metadata
metadata <- read_excel("metadata.xlsx", sheet = "Sheet1", 
                       na = c("", "NA", "N/A", "NULL", " "))

cat("📊 Metadata dimensions:", dim(metadata), "\n")

# Check missing values
cat("🔍 Missing value summary:\n")
missing_summary <- sapply(metadata, function(x) sum(is.na(x)))
print(missing_summary[missing_summary > 0])

# Antibiotic type encoding
cat("3. Encoding antibiotic types...\n")
antibiotic_mapping <- data.frame(
  original = c("β-lactams", "β-lactams+Glycopeptides", "β-lactams+Macrolides",
               "β-lactams+Triazoles", "β-lactams+Triazoles+Glycopeptides",
               "β-lactams+Triazoles+Imidazoles+Glycopeptides",
               "β-lactams+Triazole+Macrolides", "β-lactams+Glycopeptides+Macrolides",
               "β-lactams+Triazoles+Glycopeptides+Macrolides"),
  numeric_code = 1:9
)

cat("📋 Antibiotic encoding mapping:\n")
print(antibiotic_mapping)

# Apply encoding
metadata_encoded <- metadata %>%
  mutate(
    Antimicrobial_numeric = case_when(
      Antimicrobial.type.Lable. == "β-lactams" ~ 1,
      Antimicrobial.type.Lable. == "β-lactams+Glycopeptides" ~ 2,
      Antimicrobial.type.Lable. == "β-lactams+Macrolides" ~ 3,
      Antimicrobial.type.Lable. == "β-lactams+Triazoles" ~ 4,
      Antimicrobial.type.Lable. == "β-lactams+Triazoles+Glycopeptides" ~ 5,
      Antimicrobial.type.Lable. == "β-lactams+Triazoles+Imidazoles+Glycopeptides" ~ 6,
      Antimicrobial.type.Lable. == "β-lactams+Triazole+Macrolides" ~ 7,
      Antimicrobial.type.Lable. == "β-lactams+Glycopeptides+Macrolides" ~ 8,
      Antimicrobial.type.Lable. == "β-lactams+Triazoles+Glycopeptides+Macrolides" ~ 9,
      TRUE ~ NA_real_
    )
  )

# Data cleaning
cat("4. Cleaning and preparing data...\n")
metadata_clean <- metadata_encoded %>%
  rename(
    Sample_ID = `Sample ID`,
    Delivery_mode = `Deliver way`,
    Antimicrobial_type = Antimicrobial.type.Lable.,
    Intrauterine = `Intrauterine Distress`,
    hospital_stay = `Hospital Stay(d)`
  ) %>%
  mutate(
    # Ordered categorical variables
    Group = factor(Group, levels = c(0, 1, 2, 3), ordered = TRUE),
    ACS = factor(ACS, levels = c(0, 1, 2), ordered = TRUE),
    MSF = factor(MSF, levels = c(0, 1, 2, 3), ordered = TRUE),
    
    # Unordered categorical variables
    Gender = factor(Gender, levels = c(1, 2), labels = c("Male", "Female")),
    Antimicrobial_type = as.factor(Antimicrobial_type),
    Chorioamnionitis = as.factor(Chorioamnionitis),
    Delivery_mode = as.factor(Delivery_mode),
    
    # Continuous variables
    Gestational_age = as.numeric(`Gestational.age(w)`),
    Birth_Weight = as.numeric(Birth.Weight),
    hospital_stay = as.numeric(hospital_stay)
  ) %>%
  select(Sample_ID, Group, Gender, Gestational_age, Birth_Weight, Delivery_mode, 
         Antimicrobial_type, ACS, Chorioamnionitis, MSF, Intrauterine,
         EUGR, hospital_stay)

cat("✅ Data cleaning completed\n")

# Multiple imputation
cat("5. Performing multiple imputation...\n")
if (!require("mice")) {
  install.packages("mice")
  library(mice)
}

# Prepare data for imputation
impute_data_for_mice <- metadata_clean %>%
  select(-Sample_ID) %>%
  mutate(
    # Convert ordered factors to unordered for imputation
    Group = factor(Group, ordered = FALSE),
    ACS = factor(ACS, ordered = FALSE),
    MSF = factor(MSF, ordered = FALSE),
    across(where(is.factor), ~ factor(as.character(.))),
    across(where(is.numeric), as.numeric)
  )

# Perform imputation
set.seed(123)
mice_imputed <- mice(impute_data_for_mice, 
                     m = 5,
                     maxit = 10,
                     method = "pmm")

# Extract first imputed dataset
final_imputed_data <- complete(mice_imputed, 1)
final_imputed_data <- cbind(metadata$`Sample ID`, final_imputed_data)
colnames(final_imputed_data)[1] <- "Sample_ID"

# Convert back to appropriate types
final_data <- final_imputed_data %>%
  mutate(
    Group = factor(Group, levels = c(0, 1, 2, 3), ordered = TRUE),
    ACS = factor(ACS, levels = c(0, 1, 2), ordered = TRUE),
    MSF = factor(MSF, levels = c(0, 1, 2, 3), ordered = TRUE)
  )

cat("✅ Multiple imputation completed\n")

# Save imputed data
write.csv(final_data, "metadata_mice.csv", row.names = FALSE)
cat("💾 Imputed data saved to: metadata_mice.csv\n")

# Load keystone taxa data
cat("6. Loading keystone taxa data...\n")
keystone_data <- read_delim("keystone_taxa_detailed.txt", delim = "\t", show_col_types = FALSE)
keystone_data1 <- keystone_data[, -c(4:12)]

# Select top keystone ASVs
top_keystone_asvs <- keystone_data1 %>%
  group_by(group) %>%
  slice_max(keystone_score, n = 2) %>%
  ungroup()

cat("🔑 Top keystone ASVs:\n")
print(top_keystone_asvs[, c("group", "node", "keystone_score", "keystone_rank")])

# Load ASV abundance data
cat("7. Loading ASV abundance data...\n")
asv_abundance <- read.csv("rerafy_rel_asvtable.csv", sep = ",", row.names = 1)
asv_abundance <- as.data.frame(t(asv_abundance))
asv_abundance <- rownames_to_column(asv_abundance, "Sample_ID")

# Merge datasets
merged_data <- final_data %>%
  left_join(asv_abundance, by = "Sample_ID") %>%
  filter(!is.na(Group))

cat("📊 Merged data dimensions:", dim(merged_data), "\n")

# Prepare for multivariate analysis
cat("8. Preparing for multivariate analysis...\n")
multivariate_data <- merged_data %>%
  mutate(
    Group = factor(Group, levels = c(0, 1, 2, 3), 
                   labels = c("none", "mild", "moderate", "severe"),
                   ordered = TRUE),
    across(c(Gender, Delivery_mode, Antimicrobial_type, ACS, 
             Chorioamnionitis, MSF, Intrauterine, EUGR), as.factor)
  ) %>%
  rename(
    Intrauterine_Distress = Intrauterine,
    Hospital_Stay = hospital_stay
  )

# Standardize variables
cat("9. Standardizing variables...\n")
multivariate_data_std <- multivariate_data %>%
  mutate(
    # Standardize ASV abundances
    ASV1303_std = scale(log1p(ASV1303)),
    ASV2215_std = scale(log1p(ASV2215)),
    ASV1702_std = scale(log1p(ASV1702)),
    ASV708_std = scale(log1p(ASV708)),
    
    # Standardize continuous covariates
    Gestational_age_std = scale(Gestational_age),
    Birth_Weight_std = scale(Birth_Weight),
    Hospital_Stay_std = scale(Hospital_Stay)
  )

cat("✅ Variable standardization completed\n")

# Stepwise multivariate analysis
cat("10. Performing stepwise multivariate analysis...\n")

asv_list <- c("ASV1303", "ASV2215", "ASV1702", "ASV708")
covariates_std <- c("Gender", "Gestational_age_std", "Birth_Weight_std", 
                    "Delivery_mode", "Antimicrobial_type", "ACS", 
                    "Chorioamnionitis", "MSF", "Intrauterine_Distress", 
                    "EUGR", "Hospital_Stay_std")

variable_order <- c(
  "Gender", "Gestational_age_std", "Birth_Weight_std",
  "Delivery_mode", "ACS", "Chorioamnionitis", "MSF",
  "Antimicrobial_type", "Intrauterine_Distress", "EUGR", "Hospital_Stay_std"
)

stepwise_results <- list()

for (asv in asv_list) {
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Analyzing ASV:", asv, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  asv_std <- paste0(asv, "_std")
  stepwise_results[[asv]] <- list()
  successful_vars <- c()
  current_vars <- asv_std
  
  # Step 1: ASV only
  formula_step1 <- paste("Group ~", asv_std)
  model_step1 <- clm(as.formula(formula_step1), data = multivariate_data_std)
  status1 <- ifelse(all(!is.na(vcov(model_step1))), "✅ Success", "❌ Failed")
  cat("1. ASV only:", status1, "\n")
  
  # Add variables stepwise
  for (i in 1:length(variable_order)) {
    var <- variable_order[i]
    current_vars <- c(current_vars, var)
    formula_current <- paste("Group ~", paste(current_vars, collapse = " + "))
    
    model_current <- clm(as.formula(formula_current), data = multivariate_data_std)
    
    if (all(!is.na(vcov(model_current)))) {
      status <- "✅ Success"
      successful_vars <- c(successful_vars, var)
      cat(sprintf("%d. + %-20s %s\n", i+1, var, status))
    } else {
      status <- "❌ Failed"
      cat(sprintf("%d. + %-20s %s\n", i+1, var, status))
      current_vars <- current_vars[-length(current_vars)]
    }
  }
  
  # Fit final model with successful variables
  if (length(successful_vars) > 0) {
    final_formula <- paste("Group ~", paste(c(asv_std, successful_vars), collapse = " + "))
    final_model <- clm(as.formula(final_formula), data = multivariate_data_std)
    
    # Extract results
    model_summary <- summary(final_model)
    asv_coef <- model_summary$coefficients[asv_std, ]
    conf_int <- exp(confint(final_model)[asv_std, ])
    
    stepwise_results[[asv]][["final_result"]] <- data.frame(
      ASV = asv,
      Estimate = asv_coef["Estimate"],
      Std_Error = asv_coef["Std. Error"],
      z_value = asv_coef["z value"],
      p_value = asv_coef["Pr(>|z|)"],
      OR = exp(asv_coef["Estimate"]),
      CI_lower = conf_int[1],
      CI_upper = conf_int[2],
      Successful_Vars = length(successful_vars),
      stringsAsFactors = FALSE
    )
    
    cat("📊 Final results for", asv, ":\n")
    cat("   OR:", round(exp(asv_coef["Estimate"]), 4), "\n")
    cat("   95% CI: [", round(conf_int[1], 4), ", ", round(conf_int[2], 4), "]\n", sep = "")
    cat("   p-value:", round(asv_coef["Pr(>|z|)"], 4), "\n")
  }
}

# Create summary table
cat("\n", rep("=", 80), "\n", sep = "")
cat("SUMMARY OF ALL ASV ANALYSES\n")
cat(rep("=", 80), "\n", sep = "")

final_summary <- data.frame()
for (asv in asv_list) {
  if (!is.null(stepwise_results[[asv]][["final_result"]])) {
    final_summary <- rbind(final_summary, stepwise_results[[asv]][["final_result"]])
  }
}

# Add significance and interpretation
if (nrow(final_summary) > 0) {
  final_summary <- final_summary %>%
    mutate(
      Significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        p_value < 0.1 ~ ".",
        TRUE ~ "ns"
      ),
      Direction = ifelse(OR > 1, "Risk Increase", "Risk Decrease"),
      Interpretation = paste0(ifelse(OR > 1, "+", "-"), 
                             round(abs(OR-1)*100, 1), "%")
    )
  
  print(final_summary[, c("ASV", "OR", "CI_lower", "CI_upper", "p_value", "Significance", 
                          "Successful_Vars")])
  
  # Save results
  write.csv(final_summary, "keystone_multivariate_results.csv", row.names = FALSE)
  cat("💾 Results saved to: keystone_multivariate_results.csv\n")
}

# Create publication-quality plot
if (nrow(final_summary) > 0 && require(ggplot2) && require(patchwork)) {
  cat("11. Creating publication-quality visualization...\n")
  
  enhanced_plot <- ggplot(final_summary, aes(x = OR, y = reorder(ASV, OR))) +
    geom_rect(aes(xmin = CI_lower, xmax = CI_upper, 
                  ymin = as.numeric(reorder(ASV, OR)) - 0.3, 
                  ymax = as.numeric(reorder(ASV, OR)) + 0.3),
              fill = "gray90", alpha = 0.5) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), 
                   height = 0.2, color = "gray40", size = 0.8) +
    geom_point(aes(fill = Significance, size = abs(log(OR))), 
               shape = 21, color = "black", stroke = 0.5) +
    geom_text(data = subset(final_summary, Significance == "ns"),
              aes(label = "ns"), 
              vjust = -1.5, size = 3.5, color = "#33A02C", fontface = "bold") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 0.8) +
    scale_fill_manual(
      values = c("***" = "#E31A1C", "**" = "#FF7F00", "*" = "#1F78B4", 
                 "." = "#6A3D9A", "ns" = "#33A02C")
    ) +
    scale_size_continuous(range = c(3, 8), guide = "none") +
    labs(
      title = "Association of Keystone Taxa with BPD Severity",
      subtitle = "Multivariate ordinal logistic regression",
      x = "Odds Ratio per Standard Deviation Increase in Abundance",
      y = "Keystone Amplicon Sequence Variants (ASVs)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, color = "gray40", hjust = 0.5)
    ) +
    scale_x_continuous(
      trans = "log10",
      breaks = c(0.5, 0.67, 1, 1.5, 2, 3),
      labels = c("0.5", "0.67", "1.0", "1.5", "2.0", "3.0")
    )
  
  ggsave("keystone_association_plot.png", enhanced_plot, 
         width = 12, height = 8, dpi = 300, bg = "white")
  cat("🖼️  Visualization saved to: keystone_association_plot.png\n")
}

cat("\n", "="^50, "\n")
cat("KEYSTONE MULTIVARIATE ANALYSIS COMPLETED\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("="^50, "\n")
