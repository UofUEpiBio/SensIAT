#!/usr/bin/env Rscript
## Quick benchmark: just load package and run simulator benchmark
## No bootstrap, no full model fitting - just timing comparisons

library(SensIAT)
library(dplyr)
library(splines)

cat("\n========== QUICK BENCHMARK TEST ==========\n\n")

# Load the model_log from the saved RData if it exists
if (file.exists("model_log.RData")) {
    cat("Loading pre-saved model_log.RData...\n")
    load("model_log.RData")
} else {
    # Fit the models quickly
    cat("Fitting models (first time only)...\n")
    
    # Same data generation as in main script
    set.seed(42)
    n_subjects <- 25
    end_time <- 365
    
    # --- Outcome simulator (log link) ---
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
    
    # --- Intensity function ---
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
    
    # Simulate data
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
    
    # Fit models
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
    
    # Save for future runs
    save(model_log, file = "model_log.RData")
    cat("Models fitted and saved to model_log.RData\n")
}

outcome_model <- model_log$models$outcome
# Extract coefficients directly
sampled_coef <- outcome_model$coef

cat("\nOutcome model loaded. Running benchmark...\n\n")

# ============================================================================
# VERSION 1: ORIGINAL (NON-OPTIMIZED) SIMULATOR
# ============================================================================
make_parametric_single_index_simulator_v1_original <- function(outcome_model, sampled_coef, covariate_mapping = NULL) {
    if (!inherits(outcome_model, "SensIAT::Single-index-outcome-model")) 
        stop("outcome_model must be a fitted single-index outcome model")
    
    bandwidth <- if (!is.null(outcome_model$bandwidth)) outcome_model$bandwidth else stop("bandwidth missing in outcome_model")
    data_orig <- outcome_model$data
    default_names <- c(prev_outcome = "..prev_outcome..", time = "..time..", delta_time = "..delta_time..")
    if (!is.null(covariate_mapping)) {
        if (!is.character(covariate_mapping) || is.null(names(covariate_mapping))) 
            stop("covariate_mapping must be a named character vector")
        default_names[names(covariate_mapping)] <- covariate_mapping
    }

    model_terms <- delete.response(terms(outcome_model))
    y_resp <- model.response(model.frame(outcome_model))
    y_seq <- sort(unique(y_resp))

    function(prev_outcome, time, delta_time, newdata = NULL) {
        if (!is.null(newdata)) {
            if (!is.data.frame(newdata) || nrow(newdata) < 1) stop("newdata must be a data.frame with at least one row")
            nd <- newdata[1, , drop = FALSE]
        } else {
            nd <- data.frame(prev_outcome = prev_outcome, time = time, delta_time = delta_time, stringsAsFactors = FALSE)
            names(nd) <- unname(default_names)
        }

        # RECOMPUTES EVERY CALL (BOTTLENECK)
        X_new <- tryCatch(
            model.matrix(model_terms, data = nd), 
            error = function(e) stop("Failed to build model matrix for outcome newdata: ", conditionMessage(e))
        )
        X_orig <- model.matrix(model_terms, data = data_orig)
        lp0 <- as.vector(X_orig %*% as.numeric(sampled_coef[colnames(X_orig)]))
        xb <- as.vector(X_new %*% as.numeric(sampled_coef[colnames(X_new)]))

        Fhat <- pcoriaccel_NW(Xb = lp0, Y = y_resp, xb = xb[1], y_seq = y_seq, h = bandwidth, kernel = attr(outcome_model, "kernel"))
        Fhat_vec <- as.vector(Fhat)
        pmf <- pmax(0, diff(c(0, Fhat_vec)))
        if (sum(pmf) <= 0) return(mean(y_resp, na.rm = TRUE))
        pmf <- pmf / sum(pmf)
        sampled <- sample(y_seq, size = 1, prob = pmf)
        return(sampled)
    }
}

# ============================================================================
# VERSION 2: OPTIMIZED (Option 1: Pre-compute lp0)
# ============================================================================
make_parametric_single_index_simulator_v2_precompute <- function(outcome_model, sampled_coef, covariate_mapping = NULL) {
    if (!inherits(outcome_model, "SensIAT::Single-index-outcome-model")) 
        stop("outcome_model must be a fitted single-index outcome model")
    
    bandwidth <- if (!is.null(outcome_model$bandwidth)) outcome_model$bandwidth else stop("bandwidth missing in outcome_model")
    data_orig <- outcome_model$data
    default_names <- c(prev_outcome = "..prev_outcome..", time = "..time..", delta_time = "..delta_time..")
    if (!is.null(covariate_mapping)) {
        if (!is.character(covariate_mapping) || is.null(names(covariate_mapping))) 
            stop("covariate_mapping must be a named character vector")
        default_names[names(covariate_mapping)] <- covariate_mapping
    }

    model_terms <- delete.response(terms(outcome_model))
    y_resp <- model.response(model.frame(outcome_model))
    y_seq <- sort(unique(y_resp))
    kernel_type <- attr(outcome_model, "kernel")
    
    # PRE-COMPUTE ONCE (OPTIMIZATION 1)
    X_orig <- model.matrix(model_terms, data = data_orig)
    lp0 <- as.vector(X_orig %*% as.numeric(sampled_coef[colnames(X_orig)]))
    
    function(prev_outcome, time, delta_time, newdata = NULL) {
        if (!is.null(newdata)) {
            if (!is.data.frame(newdata) || nrow(newdata) < 1) stop("newdata must be a data.frame with at least one row")
            nd <- newdata[1, , drop = FALSE]
        } else {
            nd <- data.frame(prev_outcome = prev_outcome, time = time, delta_time = delta_time, stringsAsFactors = FALSE)
            names(nd) <- unname(default_names)
        }

        X_new <- tryCatch(
            model.matrix(model_terms, data = nd), 
            error = function(e) stop("Failed to build model matrix for outcome newdata: ", conditionMessage(e))
        )
        xb <- as.vector(X_new %*% as.numeric(sampled_coef[colnames(X_new)]))

        # lp0 is pre-computed, only xb computed per call
        Fhat <- pcoriaccel_NW(Xb = lp0, Y = y_resp, xb = xb[1], y_seq = y_seq, h = bandwidth, kernel = kernel_type)
        Fhat_vec <- as.vector(Fhat)
        pmf <- pmax(0, diff(c(0, Fhat_vec)))
        if (sum(pmf) <= 0) return(mean(y_resp, na.rm = TRUE))
        pmf <- pmf / sum(pmf)
        sampled <- sample(y_seq, size = 1, prob = pmf)
        return(sampled)
    }
}

# ============================================================================
# BENCHMARK FUNCTION
# ============================================================================
benchmark_simulator <- function(simulator, n_obs = 200, version_name = "Version") {
    cat(sprintf("Testing %s (%d observations)...\n", version_name, n_obs))
    
    # Generate random input sequences
    set.seed(123)
    prev_outcomes <- rlnorm(n_obs, meanlog = log(3), sdlog = 0.4)
    times <- runif(n_obs, 0, 365)
    delta_times <- pmax(1, rexp(n_obs, 0.01))  # inter-visit gaps
    
    # Warm up
    try({
        for (i in 1:5) simulator(prev_outcomes[i], times[i], delta_times[i])
    }, silent = TRUE)
    
    # Time the simulation
    t_start <- proc.time()
    
    # Scalar version: call per observation
    results <- numeric(n_obs)
    for (i in seq_len(n_obs)) {
        results[i] <- simulator(prev_outcomes[i], times[i], delta_times[i])
    }
    
    elapsed <- (proc.time() - t_start)[[3]]
    
    cat(sprintf("  Time: %.3f seconds\n", elapsed))
    cat(sprintf("  Per observation: %.4f ms\n", (elapsed / n_obs) * 1000))
    
    invisible(list(elapsed = elapsed, results = results, per_obs_ms = (elapsed / n_obs) * 1000))
}

# ============================================================================
# RUN BENCHMARKS
# ============================================================================
cat("\n=== Building simulators ===\n\n")

t_build <- proc.time()
sim_v1 <- make_parametric_single_index_simulator_v1_original(outcome_model, sampled_coef)
elapsed_v1_build <- (proc.time() - t_build)[[3]]
cat(sprintf("v1_original built in %.4f s\n", elapsed_v1_build))

t_build <- proc.time()
sim_v2 <- make_parametric_single_index_simulator_v2_precompute(outcome_model, sampled_coef)
elapsed_v2_build <- (proc.time() - t_build)[[3]]
cat(sprintf("v2_precompute built in %.4f s\n\n", elapsed_v2_build))

# Run benchmarks with increasing observation counts
for (n_obs in c(50, 200, 500)) {
    cat(sprintf("\n>>> Benchmark with n_observations = %d <<<\n\n", n_obs))
    
    result_v1 <- benchmark_simulator(sim_v1, n_obs, "v1_original")
    result_v2 <- benchmark_simulator(sim_v2, n_obs, "v2_precompute")
    
    speedup_v2 <- result_v1$elapsed / result_v2$elapsed
    
    cat(sprintf("\nSpeedup (v2 vs v1): %.2f x\n", speedup_v2))
}

cat("\n========== BENCHMARK COMPLETE ==========\n")
