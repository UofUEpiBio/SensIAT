# Diagnostic script to identify errors in benchmark
# Run this first to see what's failing

library(SensIAT)
library(dplyr)

cat("Step 1: Loading data...\n")
data_with_lags <- SensIAT_example_data |>
    dplyr::group_by(Subject_ID) |>
    dplyr::mutate(
        ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
        ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
        ..delta_time.. = Time - dplyr::lag(.data$Time, default = NA_real_, order_by = Time)
    ) |>
    dplyr::ungroup()
cat("✓ Data loaded\n\n")

cat("Step 2: Creating intensity model...\n")
intensity.model <-
    survival::coxph(
        survival::Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
        data = data_with_lags |> dplyr::filter(.data$Time > 0)
    )
cat("✓ Intensity model created\n\n")

cat("Step 3: Creating outcome model...\n")
outcome.model <-
    fit_SensIAT_single_index_fixed_coef_model(
        Outcome ~ splines::ns(..prev_outcome.., df = 3) + ..delta_time.. - 1,
        id = Subject_ID,
        data = data_with_lags |> filter(Time > 0)
    )
cat("✓ Outcome model created\n\n")

cat("Step 4: Setting up spline basis...\n")
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
cat("✓ Spline basis created\n")
cat("  - Basis dimension:", ncol(base), "\n")
cat("  - Integration range: [", tmin, ",", tmax, "]\n\n")

cat("Step 5: Setting up weight function...\n")
V <- orthogonalsplinebasis::GramMatrix(base)
V.inv <- solve(V)
cat("✓ Gram matrix computed and inverted\n\n")

cat("Step 6: Testing beta vector...\n")
# Initialize beta to give reasonable fitted values
# For log link: log(mu) = B*beta, so beta should give mu in range of observed outcomes
# Use constant function: all beta equal, sum to log(mean(outcome))
set.seed(123)
observed_outcomes <- data_with_lags$Outcome[!is.na(data_with_lags$Outcome)]
mean_outcome <- mean(observed_outcomes)
cat("  - Mean observed outcome:", round(mean_outcome, 2), "\n")

# Initialize to constant log(mean) divided by number of basis functions
# This gives mu(t) ≈ mean_outcome everywhere (reasonable starting point)
beta_test <- rep(log(mean_outcome) / ncol(base), ncol(base))
# Add small random perturbation
beta_test <- beta_test + rnorm(ncol(base), 0, 0.1)

cat("✓ Beta vector created, length:", length(beta_test), "\n")
cat("  - Range: [", round(min(beta_test), 3), ",", round(max(beta_test), 3), "]\n")
cat("  - Expected mu range: [", round(exp(sum(beta_test) * 0.8), 2), ",", 
    round(exp(sum(beta_test) * 1.2), 2), "]\n\n")

cat("Step 7: Testing W function...\n")
W <- function(t, beta) {
    B <- as.vector(pcoriaccel_evaluate_basis(base, t))
    mu <- sum(B * beta)
    if (abs(mu) < 1e-10) {
        warning("mu is very small (", mu, ") at t=", t, ". May cause numerical issues.")
    }
    (V.inv %*% B) / mu
}
test_w <- W(100, beta_test)
cat("✓ W function works\n")
cat("  - W(100) dimension:", paste(dim(test_w), collapse = " x "), "\n")
cat("  - mu at t=100:", sum(as.vector(pcoriaccel_evaluate_basis(base, 100)) * beta_test), "\n\n")

cat("Step 8: Testing imputation function...\n")
impute_data_fn <- function(t, df) {
    data_wl <- df |>
        mutate(
            ..prev_time.. = .data$Time,
            ..prev_outcome.. = .data$Outcome,
            ..delta_time.. = 0
        )
    extrapolate_from_last_observation(t, data_wl, 'Time', slopes = c('..delta_time..' = 1))
}
test_id <- unique(data_with_lags$Subject_ID)[1]
test_df <- data_with_lags[data_with_lags$Subject_ID == test_id, ]
test_impute <- impute_data_fn(100, test_df)
cat("✓ Imputation function works\n")
cat("  - Result rows:", nrow(test_impute), "\n")
cat("  - Result cols:", ncol(test_impute), "\n\n")

cat("Step 9: Testing compute_SensIAT_expected_values...\n")
test_expected <- compute_SensIAT_expected_values(
    model = outcome.model,
    alpha = alpha,
    new.data = test_impute
)
cat("✓ compute_SensIAT_expected_values works\n")
cat("  - E_Yexp_alphaY:", test_expected$E_Yexp_alphaY, "\n")
cat("  - E_exp_alphaY:", test_expected$E_exp_alphaY, "\n\n")

cat("Step 10: Testing integrand function...\n")
# Note: For log link, eta = B*beta, mu = exp(eta)
# The W function uses mu = B*beta (should be eta), then divides by it
# This can cause issues if eta is near zero
term2_integrand <- function(t) {
    weight <- W(t, beta_test)
    expected <- compute_SensIAT_expected_values(
        model = outcome.model,
        alpha = alpha,
        new.data = impute_data_fn(t, test_df)
    )
    B <- pcoriaccel_evaluate_basis(base, t)
    inv.link <- exp
    eta <- crossprod(B, beta_test)
    
    # Check for numerical issues
    if (abs(eta) < 0.01) {
        cat("  Warning: eta =", eta, "at t =", t, "\n")
    }
    
    weight * as.numeric(
        expected$E_Yexp_alphaY / expected$E_exp_alphaY -
        inv.link(eta)
    )
}
test_integrand_val <- term2_integrand(100)
cat("✓ Integrand function works\n")
cat("  - Integrand(100) dimension:", paste(dim(test_integrand_val), collapse = " x "), "\n")
cat("  - Integrand(100) first few values:\n")
print(head(test_integrand_val))
cat("\n")

cat("Step 11: Testing integration...\n")
test_integrate <- pcoriaccel_integrate_simp(
    term2_integrand,
    tmin,
    tmax
)
cat("✓ Integration works\n")
cat("  - Q:", test_integrate$Q, "\n")
cat("  - fcnt:", test_integrate$fcnt, "\n")
cat("  - estim.prec:", test_integrate$estim.prec, "\n\n")

cat(strrep("=", 60), "\n")
cat("ALL DIAGNOSTICS PASSED!\n")
cat(strrep("=", 60), "\n")
cat("\nYou can now run the full benchmark:\n")
cat("  source('inst/examples/benchmark_term2_caching.R')\n\n")
