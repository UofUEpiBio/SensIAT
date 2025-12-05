# Minimal benchmark test - can paste directly into R console
# Assumes you've run devtools::load_all() or library(SensIAT)

cat(strrep("=", 60), "\n")
cat("MINIMAL BENCHMARK TEST\n")
cat(strrep("=", 60), "\n\n")

# Load packages
suppressPackageStartupMessages({
    # Prefer loading local source (devtools::load_all) when in repo to pick up latest changes
    if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
        devtools::load_all(quiet = TRUE)
    } else if (requireNamespace("SensIAT", quietly = TRUE)) {
        library(SensIAT)
    } else {
        stop("Cannot load SensIAT: install the package or run within the repo with devtools available.")
    }
    library(dplyr)
    library(purrr)
})

# Ensure required dependent packages are available when running stand-alone
if (!requireNamespace("orthogonalsplinebasis", quietly = TRUE)) {
    cat("Installing orthogonalsplinebasis package...\\n")
    install.packages("orthogonalsplinebasis")
}
if (!requireNamespace("survival", quietly = TRUE)) {
    cat("Installing survival package...\\n")
    install.packages("survival")
}

# Setup
cat("Setting up data and models...\n")
data_with_lags <- SensIAT_example_data |>
    group_by(Subject_ID) |>
    mutate(
        ..prev_outcome.. = lag(Outcome, default = NA_real_, order_by = Time),
        ..prev_time.. = lag(Time, default = 0, order_by = Time),
        ..delta_time.. = Time - lag(Time, default = NA_real_, order_by = Time)
    ) |>
    ungroup()

intensity.model <- survival::coxph(
    survival::Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
    data = data_with_lags |> filter(Time > 0)
)

outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
    Outcome ~ splines::ns(..prev_outcome.., df = 3) + ..delta_time.. - 1,
    id = Subject_ID,
    data = data_with_lags |> filter(Time > 0)
)

alpha <- 0
knots <- c(60, 260, 460)
spline.degree <- 3L
knots_extended <- c(rep(head(knots, 1), spline.degree), knots, rep(tail(knots, 1), spline.degree))
base <- orthogonalsplinebasis::SplineBasis(knots_extended, order = spline.degree + 1L)
tmin <- base@knots[base@order]
tmax <- base@knots[length(base@knots) - base@order + 1]

V <- orthogonalsplinebasis::GramMatrix(base)
V.inv <- solve(V)

set.seed(123)
observed_outcomes <- data_with_lags$Outcome[!is.na(data_with_lags$Outcome)]
mean_outcome <- mean(observed_outcomes)
beta_test <- rep(log(mean_outcome) / ncol(base), ncol(base)) + rnorm(ncol(base), 0, 0.1)

test_id <- unique(data_with_lags$Subject_ID)[1]
cat("✓ Setup complete\n\n")

# Test FAST integrand
cat("Testing FAST integrand builder...\n")
data_i <- data_with_lags[data_with_lags$Subject_ID == test_id, ]

tryCatch({
    # Access the fast integrand builder from the SensIAT namespace if not exported
    f_builder <- tryCatch(
        get("make_term2_integrand_fast", envir = asNamespace("SensIAT")),
        error = function(e) get("make_term2_integrand_fast")
    )

    integrand_fast <- f_builder(
        outcome.model = outcome.model,
        base = base,
        alpha = alpha,
        patient_times = data_i$Time,
        patient_outcomes = data_i$Outcome,
        marginal_beta = beta_test,
        V_inv = V.inv
    )
    cat("✓ Fast integrand created\n")
    
    # Test at a single point
    test_val <- integrand_fast(100)
    cat("✓ Integrand evaluated at t=100\n")
    cat("  Result length:", length(test_val), "\n")
    cat("  First few values:", head(test_val, 3), "\n\n")
    
    # Test integration
    cat("Testing integration...\n")
    result <- pcoriaccel_integrate_simp(integrand_fast, tmin, tmax)
    cat("✓ Integration successful\n")
    cat("  Q:", result$Q[1:3], "...\n")
    cat("  fcnt:", result$fcnt, "\n")
    cat("  estim.prec:", result$estim.prec, "\n\n")
    
    cat(strrep("=", 60), "\n")
    cat("SUCCESS! Fast integrand works correctly.\n")
    cat("You can now run the full benchmark:\n")
    cat("  source('inst/examples/benchmark_term2_caching.R')\n")
    cat(strrep("=", 60), "\n")
    
}, error = function(e) {
    cat("✗ ERROR:\n")
    cat(conditionMessage(e), "\n\n")
    cat("Traceback:\n")
    print(sys.calls())
})
