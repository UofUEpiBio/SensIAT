# Benchmark: Original vs Cached Term 2 Influence Computation
#
# This script compares the performance of the original method for computing
# term 2 influence values against an optimized cached version.

library(dplyr)
library(purrr)

# Prefer loading local source to pick up latest code when running from repo
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(quiet = TRUE)
} else {
    library(SensIAT)
}

# Check for tictoc, install if needed
if (!requireNamespace("tictoc", quietly = TRUE)) {
    cat("Installing tictoc package...\n")
    install.packages("tictoc")
}
library(tictoc)

# Setup test data ----
# Using the example data from the package
data_with_lags <- SensIAT_example_data |>
    dplyr::group_by(Subject_ID) |>
    dplyr::mutate(
        ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
        ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
        ..delta_time.. = Time - dplyr::lag(.data$Time, default = NA_real_, order_by = Time)
    ) |>
    dplyr::ungroup()

# Create models ----
intensity.model <-
    survival::coxph(
        survival::Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
        data = data_with_lags |> dplyr::filter(.data$Time > 0)
    )

outcome.model <-
    fit_SensIAT_single_index_fixed_coef_model(
        Outcome ~ splines::ns(..prev_outcome.., df = 3) + ..delta_time.. - 1,
        id = Subject_ID,
        data = data_with_lags |> filter(Time > 0)
    )

# Setup parameters ----
alpha <- 0
knots <- c(60, 260, 460)
spline.degree <- 3L

knots_extended <- c(
    rep(head(knots, 1), spline.degree),
    knots,
    rep(tail(knots, 1), spline.degree)
)
base <- orthogonalsplinebasis::SplineBasis(knots_extended, order = spline.degree + 1L)

tmin <- base@knots[base@order]
tmax <- base@knots[length(base@knots) - base@order + 1]

# Link function setup (log link)
link.fun <- log
inv.link <- exp
d1.inv.link <- exp

# Weight function for lp_mse with log link
V <- orthogonalsplinebasis::GramMatrix(base)
V.inv <- solve(V)

W <- function(t, beta) {
    B <- as.vector(pcoriaccel_evaluate_basis(base, t))
    mu <- sum(B * beta)  # More robust than crossprod for this case
    (V.inv %*% B) / mu
}

# Imputation function
impute_data_fn <- function(t, df) {
    data_wl <- df |>
        mutate(
            ..prev_time.. = .data$Time,
            ..prev_outcome.. = .data$Outcome,
            ..delta_time.. = 0
        )
    extrapolate_from_last_observation(t, data_wl, 'Time', slopes = c('..delta_time..' = 1))
}

# Test beta vector - initialize to give reasonable mu values
# For log link: log(mu) = B*beta, initialize to constant log(mean_outcome)
set.seed(123)
observed_outcomes <- data_with_lags$Outcome[!is.na(data_with_lags$Outcome)]
mean_outcome <- mean(observed_outcomes)
beta_test <- rep(log(mean_outcome) / ncol(base), ncol(base))
beta_test <- beta_test + rnorm(ncol(base), 0, 0.1)  # Small random perturbation

# Select a subset of patients for testing
test_ids <- unique(data_with_lags$Subject_ID)[1:5]  # Use first 5 patients
cat("Testing with", length(test_ids), "patients\n\n")


# ORIGINAL METHOD ----
compute_term2_for_patient_original <- function(patient_id, data, id_var, beta) {
    data_i <- data[id_var == patient_id, ]

    term2_integrand <- function(t) {
        weight <- W(t, beta)
        expected <- compute_SensIAT_expected_values(
            model = outcome.model,
            alpha = alpha,
            new.data = impute_data_fn(t, data_i)
        )
        B <- pcoriaccel_evaluate_basis(base, t)
        weight * as.numeric(
            expected$E_Yexp_alphaY / expected$E_exp_alphaY -
            inv.link(crossprod(B, beta))
        )
    }

    rslt <- pcoriaccel_integrate_simp(
        term2_integrand,
        tmin,
        tmax
    )
    rslt$Q
}


# CACHED METHOD ----
compute_term2_for_patient_cached <- function(patient_id, data, id_var, beta) {
    data_i <- data[id_var == patient_id, ]

    # Cache for storing computed expectations
    expected_cache <- new.env(hash = TRUE)
    # Note: R CMD check may warn about unused variable, but it IS used in the closure below

    term2_integrand <- function(t) {
        # Create cache key (round to avoid floating point issues)
        cache_key <- as.character(round(t, 10))

        if (!exists(cache_key, envir = expected_cache)) {
            weight <- W(t, beta)
            imputed <- impute_data_fn(t, data_i)
            expected <- compute_SensIAT_expected_values(
                model = outcome.model,
                alpha = alpha,
                new.data = imputed
            )
            B <- pcoriaccel_evaluate_basis(base, t)

            # Store in cache
            expected_cache[[cache_key]] <- list(
                weight = weight,
                expected = expected,
                B = B
            )
        }

        cached <- expected_cache[[cache_key]]
        cached$weight * as.numeric(
            cached$expected$E_Yexp_alphaY / cached$expected$E_exp_alphaY -
            inv.link(crossprod(cached$B, beta))
        )
    }

    rslt <- pcoriaccel_integrate_simp(
        term2_integrand,
        tmin,
        tmax
    )
    rslt$Q
}


# BENCHMARKING ----
cat(strrep("=", 60), "\n")
cat("BENCHMARK: Term 2 Influence Computation\n")
cat(strrep("=", 60), "\n\n")

cat("Setup:\n")
cat("  - Number of patients:", length(test_ids), "\n")
cat("  - Alpha:", alpha, "\n")
cat("  - Integration range: [", tmin, ",", tmax, "]\n")
cat("  - Spline basis dimension:", ncol(base), "\n\n")

# Warm-up run (to compile any lazy functions)
cat("Warming up...\n")
tryCatch({
    invisible(compute_term2_for_patient_original(test_ids[1], data_with_lags, data_with_lags$Subject_ID, beta_test))
    invisible(compute_term2_for_patient_cached(test_ids[1], data_with_lags, data_with_lags$Subject_ID, beta_test))
    cat("Warm-up complete.\n\n")
}, error = function(e) {
    cat("ERROR during warm-up:\n")
    cat(conditionMessage(e), "\n")
    cat("\nDebug info:\n")
    cat("  - test_ids[1]:", test_ids[1], "\n")
    cat("  - beta length:", length(beta_test), "\n")
    cat("  - base ncol:", ncol(base), "\n")
    stop("Warm-up failed. Fix errors before benchmarking.", call. = FALSE)
})

# Run original method
cat("Running ORIGINAL method...\n")
tic("Original method")
results_original <- map(
    test_ids,
    compute_term2_for_patient_original,
    data = data_with_lags,
    id_var = data_with_lags$Subject_ID,
    beta = beta_test
)
time_original <- toc(quiet = TRUE)
cat("Original method completed in", round(time_original$toc - time_original$tic, 2), "seconds\n\n")

# Run cached method
cat("Running CACHED method...\n")
cached_error <- NULL
results_cached <- NULL
time_cached <- NULL
tryCatch({
    tic("Cached method")
    results_cached <- map(
        test_ids,
        compute_term2_for_patient_cached,
        data = data_with_lags,
        id_var = data_with_lags$Subject_ID,
        beta = beta_test
    )
    time_cached <- toc(quiet = TRUE)
    cat("Cached method completed in", round(time_cached$toc - time_cached$tic, 2), "seconds\n\n")
}, error = function(e) {
    cached_error <<- conditionMessage(e)
    cat("ERROR in CACHED method:\n", cached_error, "\n\n", sep = "")
})

# Run FAST method (closure-based, low allocation)
cat("Running FAST method...\n")
make_fast_integrand_for_patient <- function(patient_id) {
    data_i <- data_with_lags[data_with_lags$Subject_ID == patient_id, ]
    # Access the fast integrand builder from the SensIAT namespace if not exported
    f_builder <- tryCatch(
        get("make_term2_integrand_fast", envir = asNamespace("SensIAT")),
        error = function(e) get("make_term2_integrand_fast")
    )
    f_builder(
        outcome.model = outcome.model,
        base = base,
        alpha = alpha,
        patient_times = data_i$Time,
        patient_outcomes = data_i$Outcome,
        marginal_beta = beta_test,
        V_inv = V.inv
    )
}

fast_error <- NULL
results_fast <- NULL
time_fast <- NULL
tryCatch({
    tic("Fast method")
    results_fast <- map(
        test_ids,
        ~{
            f <- make_fast_integrand_for_patient(.x)
            rslt <- pcoriaccel_integrate_simp(f, tmin, tmax)
            rslt$Q
        }
    )
    time_fast <- toc(quiet = TRUE)
    cat("FAST method completed in", round(time_fast$toc - time_fast$tic, 2), "seconds\n\n")
}, error = function(e) {
    fast_error <<- conditionMessage(e)
    cat("ERROR in FAST method:\n", fast_error, "\n\n", sep = "")
})

# VALIDATION ----
cat(strrep("=", 60), "\n")
cat("VALIDATION: Comparing Results\n")
cat(strrep("=", 60), "\n\n")

# Check that results match (element-wise)
if (!is.null(results_cached)) {
    results_match <- all.equal(unlist(results_original), unlist(results_cached), tolerance = 1e-10)
    cat("Results match:", isTRUE(results_match), "\n")
    if (!isTRUE(results_match)) {
        cat("Differences:\n")
        print(results_match)
    }
} else {
    cat("Results match: skipped (cached method failed)\n")
}

# Compare individual patient results
orig_vec <- unlist(results_original)
if (!is.null(results_cached)) {
    cached_vec <- unlist(results_cached)
    comparison_df <- tibble(
        patient_id = test_ids,
        original = orig_vec,
        cached = cached_vec,
        difference = abs(original - cached),
        rel_diff = ifelse(abs(original) > 0, abs(original - cached) / abs(original), NA_real_)
    )
    cat("\nPer-patient comparison:\n")
    print(comparison_df, n = Inf)
} else {
    cat("\nPer-patient comparison: skipped (cached method failed)\n")
}

# SUMMARY ----
cat("\n")
cat(strrep("=", 60), "\n")
cat("PERFORMANCE SUMMARY\n")
cat(strrep("=", 60), "\n\n")

elapsed_original <- time_original$toc - time_original$tic
elapsed_cached <- if (!is.null(time_cached)) time_cached$toc - time_cached$tic else NA_real_
elapsed_fast <- if (!is.null(time_fast)) time_fast$toc - time_fast$tic else NA_real_
speedup_cached <- if (!is.na(elapsed_cached)) elapsed_original / elapsed_cached else NA_real_
speedup_fast <- if (!is.na(elapsed_fast)) elapsed_original / elapsed_fast else NA_real_

cat("Original method: ", round(elapsed_original, 2), " seconds\n", sep = "")
if (!is.na(elapsed_cached)) {
    cat("Cached method:   ", round(elapsed_cached, 2), " seconds\n", sep = "")
} else {
    cat("Cached method:   FAILED (", cached_error, ")\n", sep = "")
}
if (!is.na(elapsed_fast)) {
    cat("FAST method:     ", round(elapsed_fast, 2), " seconds\n", sep = "")
} else {
    cat("FAST method:     FAILED (", fast_error, ")\n", sep = "")
}
if (!is.na(speedup_cached)) {
    cat("Speedup (Cached vs Original): ", round(speedup_cached, 2), "x\n", sep = "")
} else {
    cat("Speedup (Cached vs Original):  n/a (cached failed)\n", sep = "")
}
if (!is.na(speedup_fast)) {
    cat("Speedup (FAST vs Original):   ", round(speedup_fast, 2), "x\n", sep = "")
    cat("Time saved (FAST):            ", round(elapsed_original - elapsed_fast, 2), " seconds (",
        round(100 * (1 - elapsed_fast/elapsed_original), 1), "%)\n", sep = "")
} else {
    cat("Speedup (FAST vs Original):    n/a (fast failed)\n", sep = "")
}

cat("\nPer-patient average:\n")
cat("Original: ", round(elapsed_original / length(test_ids), 3), " seconds/patient\n", sep = "")
cat("Cached:   ", round(elapsed_cached / length(test_ids), 3), " seconds/patient\n", sep = "")
cat("FAST:     ", round(elapsed_fast / length(test_ids), 3), " seconds/patient\n", sep = "")

# Estimate full dataset time
total_patients <- length(unique(data_with_lags$Subject_ID))
cat("\nEstimated time for all", total_patients, "patients:\n")
cat("Original: ", round(elapsed_original / length(test_ids) * total_patients / 60, 1), " minutes\n", sep = "")
if (!is.na(elapsed_cached)) {
    cat("Cached:   ", round(elapsed_cached / length(test_ids) * total_patients / 60, 1), " minutes\n", sep = "")
} else {
    cat("Cached:   n/a (failed)\n")
}
if (!is.na(elapsed_fast)) {
    cat("FAST:     ", round(elapsed_fast / length(test_ids) * total_patients / 60, 1), " minutes\n", sep = "")
} else {
    cat("FAST:     n/a (failed)\n")
}

cat("\n")
cat(strrep("=", 60), "\n")
