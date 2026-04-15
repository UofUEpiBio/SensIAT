# Quick Benchmark Runner
# Run this to quickly compare original vs cached term2 computation
# with minimal output

library(SensIAT)

# Load and run benchmark
cat("Loading benchmark script...\n")
source("inst/examples/benchmark_term2_caching.R")
