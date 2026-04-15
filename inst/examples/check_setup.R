# Quick syntax check and dry run
# Run from R console after devtools::load_all()

cat("Step 1: Load term2_fast.R\n")
source("R/term2_fast.R")
cat("✓ term2_fast.R loaded\n\n")

cat("Step 2: Check make_term2_integrand_fast exists\n")
if (exists("make_term2_integrand_fast")) {
    cat("✓ make_term2_integrand_fast function found\n")
    cat("  Args:", paste(names(formals(make_term2_integrand_fast)), collapse = ", "), "\n\n")
} else {
    stop("make_term2_integrand_fast not found!")
}

cat("Step 3: Try loading benchmark dependencies\n")
library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
if (!requireNamespace("tictoc", quietly = TRUE)) {
    cat("Installing tictoc...\n")
    install.packages("tictoc", repos = "https://cloud.r-project.org")
}
library(tictoc, warn.conflicts = FALSE)
cat("✓ Dependencies loaded\n\n")

cat("Step 4: Check for SensIAT functions\n")
required_funs <- c(
    "fit_SensIAT_single_index_fixed_coef_model",
    "pcoriaccel_evaluate_basis",
    "pcoriaccel_estimate_pmf",
    "pcoriaccel_integrate_simp",
    "compute_SensIAT_expected_values",
    "SplineBasis",
    "GramMatrix",
    "extrapolate_from_last_observation"
)
for (fn in required_funs) {
    if (exists(fn)) {
        cat("  ✓", fn, "\n")
    } else {
        cat("  ✗", fn, "NOT FOUND\n")
    }
}

cat("\nAll checks passed! Ready to run benchmark.\n")
cat("Run: source('inst/examples/benchmark_term2_caching.R')\n")
