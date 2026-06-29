#!/usr/bin/env Rscript
## Benchmark: Just measure model.matrix computation overhead
## This isolates Optimization 1 impact: pre-computing X_orig and lp0

library(SensIAT)
library(dplyr)
library(splines)

cat("\n========== MODEL MATRIX OVERHEAD BENCHMARK ==========\n\n")

# Load or fit models
if (file.exists("model_log.RData")) {
    cat("Loading pre-saved model_log.RData...\n")
    load("model_log.RData")
} else {
    cat("Fitting models (first time only)...\n")
    set.seed(42)
    n_subjects <- 25
    end_time <- 365
    
    beta_prev <- -0.35
    beta_dt   <-  0.003
    beta_time <-  3e-4
    intercept <-  1.2
    sigma_eps <-  0.25
    
    outcome_simulator <- function(prev_outcome, time, delta_time, ...) {
        log_prev <- log(max(prev_outcome, 1e-6))
        eta <- intercept +
            beta_prev * log_prev +
            beta_dt   * sqrt(delta_time) +
            beta_time * time
        exp(stats::rnorm(1, mean = eta, sd = sigma_eps))
    }
    
    baseline_outcome_fn <- function() stats::rlnorm(1, meanlog = log(3), sdlog = 0.4)
    
    wb_shape <- 1.35
    wb_scale <- 220
    gamma_y    <- 0.20
    y_ref      <- 3
    lambda_min <- 1 / 240
    lambda_max <- 1 / 35
    
    weibull_baseline_hazard <- function(t, shape = wb_shape, scale = wb_scale) {
        tt <- max(t, 1e-6)
        (shape / scale) * (tt / scale)^(shape - 1)
    }
    
    intensity_fn <- function(t, prev_outcome, visit_num) {
        h0 <- weibull_baseline_hazard(t)
        outcome_mult <- exp(gamma_y * (log(max(prev_outcome, 1e-6)) - log(y_ref)))
        rate <- h0 * outcome_mult
        min(max(rate, lambda_min), lambda_max)
    }
    intensity_bound <- lambda_max
    
    sim_data <- simulate_SensIAT_data(
        n_subjects          = n_subjects,
        End                 = end_time,
        baseline_outcome_fn = baseline_outcome_fn,
        outcome_simulator   = outcome_simulator,
        intensity_fn        = intensity_fn,
        intensity_bound     = intensity_bound,
        seed                = 42,
        max_visits          = 30
    )
    
    knots <- c(90, 180, 270)
    model_log <- fit_SensIAT_within_group_model(
        group.data = sim_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        id = "Subject_ID",
        outcome = "Outcome",
        time = "Time",
        knots = knots,
        alpha = 0,
        End = end_time,
        outcome.args = list(
            model = ~ns(..prev_outcome.., df = 3) + scale(..delta_time..) - 1
        ),
        influence.args = list(
            term2_grid_n = 20,
            BBsolve.control = list(maxit = 200, tol = 1e-4)
        ),
        loss = "lp_mse",
        link = "log",
        term2_method = "fixed_grid"
    )
    
    save(model_log, file = "model_log.RData")
    cat("Models fitted and saved\n")
}

outcome_model <- model_log$models$outcome
sampled_coef <- outcome_model$coef

cat("\nOutcome model loaded. Benchmarking model.matrix computation...\n\n")

# ============================================================================
# Extract model info
# ============================================================================
model_terms <- delete.response(terms(outcome_model))
data_orig <- outcome_model$data
default_names <- c(prev_outcome = "..prev_outcome..", time = "..time..", delta_time = "..delta_time..")

# ============================================================================
# VERSION 1: RECOMPUTE X_orig EVERY CALL
# ============================================================================
v1_simulator <- function(prev_outcome, time, delta_time) {
    # Create newdata frame
    nd <- data.frame(
        prev_outcome = prev_outcome,
        time = time,
        delta_time = delta_time,
        stringsAsFactors = FALSE
    )
    names(nd) <- unname(default_names)
    
    # BOTTLENECK: Recompute both X_new and X_orig every call
    X_new <- model.matrix(model_terms, data = nd)
    X_orig <- model.matrix(model_terms, data = data_orig)
    lp0 <- as.vector(X_orig %*% as.numeric(sampled_coef[colnames(X_orig)]))
    xb <- as.vector(X_new %*% as.numeric(sampled_coef[colnames(X_new)]))
    
    # Return dummy value for benchmarking
    return(list(lp0 = lp0, xb = xb))
}

# ============================================================================
# VERSION 2: PRE-COMPUTE X_orig ONCE (OPTIMIZATION 1)
# ============================================================================
# Set up: pre-compute X_orig and lp0 once
X_orig <- model.matrix(model_terms, data = data_orig)
lp0 <- as.vector(X_orig %*% as.numeric(sampled_coef[colnames(X_orig)]))

v2_simulator <- function(prev_outcome, time, delta_time) {
    # Create newdata frame
    nd <- data.frame(
        prev_outcome = prev_outcome,
        time = time,
        delta_time = delta_time,
        stringsAsFactors = FALSE
    )
    names(nd) <- unname(default_names)
    
    # OPTIMIZED: Only compute X_new (lp0 is pre-computed)
    X_new <- model.matrix(model_terms, data = nd)
    xb <- as.vector(X_new %*% as.numeric(sampled_coef[colnames(X_new)]))
    
    # Return using pre-computed lp0
    return(list(lp0 = lp0, xb = xb))
}

# ============================================================================
# BENCHMARK FUNCTION
# ============================================================================
benchmark_version <- function(simulator, n_obs = 200, version_name = "Version") {
    cat(sprintf("Testing %s (%d observations)...\n", version_name, n_obs))
    
    # Generate random inputs
    set.seed(123)
    prev_outcomes <- rlnorm(n_obs, meanlog = log(3), sdlog = 0.4)
    times <- runif(n_obs, 0, 365)
    delta_times <- pmax(1, rexp(n_obs, 0.01))
    
    # Warm up (JIT compilation)
    try({
        for (i in 1:5) simulator(prev_outcomes[i], times[i], delta_times[i])
    }, silent = TRUE)
    
    # Benchmark
    t_start <- proc.time()
    results <- list()
    for (i in seq_len(n_obs)) {
        results[[i]] <- simulator(prev_outcomes[i], times[i], delta_times[i])
    }
    elapsed <- (proc.time() - t_start)[[3]]
    
    per_obs_ms <- (elapsed / n_obs) * 1000
    
    cat(sprintf("  Time: %.4f seconds\n", elapsed))
    cat(sprintf("  Per observation: %.4f ms\n", per_obs_ms))
    
    invisible(list(elapsed = elapsed, per_obs_ms = per_obs_ms))
}

# ============================================================================
# RUN BENCHMARKS
# ============================================================================
cat("\n=== Benchmarking model.matrix computation overhead ===\n\n")

# Pre-compute setup time
cat("v2 (pre-computation) setup:\n")
t_setup <- proc.time()
X_orig_setup <- model.matrix(model_terms, data = data_orig)
lp0_setup <- as.vector(X_orig_setup %*% as.numeric(sampled_coef[colnames(X_orig_setup)]))
elapsed_setup <- (proc.time() - t_setup)[[3]]
cat(sprintf("  X_orig computed in %.4f s\n\n", elapsed_setup))

# Run benchmarks at increasing observation counts
for (n_obs in c(50, 200, 500)) {
    cat(sprintf("\n>>> Benchmark with n_observations = %d <<<\n\n", n_obs))
    
    result_v1 <- benchmark_version(v1_simulator, n_obs, "v1_recompute")
    result_v2 <- benchmark_version(v2_simulator, n_obs, "v2_precompute")
    
    speedup <- result_v1$elapsed / result_v2$elapsed
    
    cat(sprintf("\nSpeedup (v2 vs v1): %.2f x\n", speedup))
    cat(sprintf("Per-observation improvement: %.4f ms -> %.4f ms (%.2f ms saved per call)\n",
                result_v1$per_obs_ms, result_v2$per_obs_ms, 
                result_v1$per_obs_ms - result_v2$per_obs_ms))
}

cat("\n========== BENCHMARK COMPLETE ==========\n")
cat("\nInterpretation:\n")
cat("- v1 recomputes X_orig (entire training data) at every call\n")
cat("- v2 pre-computes X_orig once in the closure\n")
cat("- Speedup shows reduction in per-observation time\n")
cat("- With pcoriaccel_NW added, v2 will be even more dominant\n")
