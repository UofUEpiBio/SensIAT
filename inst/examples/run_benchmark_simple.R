# Run this from your R console:
# setwd("c:/Users/u0092104/OneDrive - University of Utah/Projects/Dan S - PCORI/pcoriRPackage")
# source("inst/examples/run_benchmark_simple.R")

cat("Loading required packages...\n")
library(SensIAT)
library(dplyr)
library(purrr)

if (!requireNamespace("tictoc", quietly = TRUE)) {
    install.packages("tictoc")
}
library(tictoc)

cat("Sourcing benchmark script...\n")
source("inst/examples/benchmark_term2_caching.R")
