#!/usr/bin/env Rscript

cat("Installing SensIAT package dependencies...\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Core dependencies from DESCRIPTION
imports <- c(
  "assertthat",
  "dplyr",
  "generics",
  "ggplot2",
  "glue",
  "KernSmooth",
  "MASS",
  "MAVE",
  "methods",
  "orthogonalsplinebasis",
  "pracma",
  "purrr",
  "Rcpp",
  "rlang",
  "splines",
  "stats",
  "survival",
  "tibble",
  "tidyr",
  "utils"
)

suggests <- c(
  "dfoptim",
  "furrr",
  "future",
  "inline",
  "knitr",
  "ManifoldOptim",
  "metR",
  "progress",
  "rmarkdown",
  "spelling",
  "testthat",
  "tidyverse"
)

# Additional development tools
dev_tools <- c(
  "BB"  # Used in fit functions
)

cat("\nInstalling Imports packages...\n")
install.packages(imports, dependencies = TRUE)

cat("\nInstalling Suggests packages...\n")
install.packages(suggests, dependencies = TRUE)

cat("\nInstalling additional development dependencies...\n")
install.packages(dev_tools, dependencies = TRUE)

cat("\nAll dependencies installed successfully!\n")
