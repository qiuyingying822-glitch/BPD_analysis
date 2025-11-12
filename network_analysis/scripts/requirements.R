
# Required R Packages for BPD Network Analysis
# Install these packages before running the analysis

required_packages <- c(
  "dplyr",      # Data manipulation
  "ggplot2",    # Data visualization
  "igraph",     # Network analysis
  "tidyr",      # Data tidying
  "reshape2",   # Data reshaping
  "readxl",     # Excel file reading
  "MASS"        # Statistical functions
)

# Function to check and install packages
install_if_missing <- function(packages) {
  for(pkg in packages) {
    if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing package:", pkg, "\n")
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    } else {
      cat("Package", pkg, "is already installed\n")
    }
  }
}

# Install required packages
cat("=== Installing Required R Packages ===\n")
install_if_missing(required_packages)

cat("\n=== All required packages are installed ===\n")
cat("You can now run the BPD network analysis scripts.\n")
