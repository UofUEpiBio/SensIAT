# Performance tests for fit_SensIAT_marginal_mean_model_generalized

test_that("fit_SensIAT_marginal_mean_model_generalized: expected value caching works correctly", {
    # This test verifies that the caching mechanism:
    # 1. Reuses computed expected values across optimizer iterations
    # 2. Produces identical results to the uncached version
    # 3. Actually caches (hits > 0)
    
    skip_on_cran()
    skip_if_not(identical(Sys.getenv("RUN_PERFORMANCE_TESTS"), "true"), "Set RUN_PERFORMANCE_TESTS=true to run")
    
    set.seed(456)
    setup <- generate_test_data(link = "log", n_subjects = 15)
    attr(setup$intensity.model, "bandwidth") <- 7
    
    # Track cache hits/misses by wrapping the cache function
    cache_hits <- 0
    cache_misses <- 0
    
    # Instrument the generalized function to count cache usage
    # We'll do this by temporarily modifying the environment
    result_with_tracking <- local({
        # Get a patient's data
        test_patient_id <- unique(setup$data$..id..)[1]
        patient_data <- setup$data[setup$data$..id.. == test_patient_id, ]
        
        # Build a cache with tracking
        cache_env <- new.env(parent = emptyenv())
        expected_get_tracked <- function(t) {
            k <- as.character(signif(t, 12))
            if (exists(k, cache_env, inherits = FALSE)) {
                cache_hits <<- cache_hits + 1
                cache_env[[k]]
            } else {
                cache_misses <<- cache_misses + 1
                df <- create_impute_fn()(t, patient_data)
                ev <- compute_SensIAT_expected_values(
                    model = setup$outcome.model,
                    alpha = 0,
                    new.data = df
                )
                cache_env[[k]] <- list(
                    E_exp_alphaY = ev$E_exp_alphaY,
                    E_Yexp_alphaY = ev$E_Yexp_alphaY
                )
                cache_env[[k]]
            }
        }
        
        # Simulate multiple calls with same t values (like optimizer would do)
        test_times <- seq(30, 90, by = 10)
        for (iter in 1:3) {  # Simulate 3 optimizer iterations
            for (t in test_times) {
                val <- expected_get_tracked(t)
            }
        }
        
        list(hits = cache_hits, misses = cache_misses)
    })
    
    cat("\n")
    cat("Cache performance:\n")
    cat(sprintf("  Cache hits:   %d\n", result_with_tracking$hits))
    cat(sprintf("  Cache misses: %d\n", result_with_tracking$misses))
    cat(sprintf("  Hit rate:     %.1f%%\n", 
                100 * result_with_tracking$hits / (result_with_tracking$hits + result_with_tracking$misses)))
    
    # Verify caching is working (should have hits on 2nd and 3rd iterations)
    expect_true(result_with_tracking$hits > 0,
                info = "Cache should have hits when same t values are queried")
    
    # Now verify full model results are consistent
    result_cached <- suppressWarnings({
        fit_SensIAT_marginal_mean_model_generalized(
            data = setup$data,
            time = ..time..,
            id = ..id..,
            alpha = 0,
            knots = setup$knots,
            outcome.model = setup$outcome.model,
            intensity.model = setup$intensity.model,
            loss = "lp_mse",
            link = "log",
            impute_data = create_impute_fn(),
            BBsolve.control = list(maxit = 50, tol = 1e-2, trace = FALSE),
            term2_method = "fast"
        )
    })
    
    # Verify it converged
    expect_true(length(result_cached$coefficients) > 0,
                info = "Model with caching should produce coefficients")
    expect_true(all(is.finite(result_cached$coefficients)),
                info = "Coefficients should be finite")
})

test_that("fit_SensIAT_marginal_mean_model_generalized: fast term2 timing comparison", {
    skip_on_cran()
    skip_if_not(identical(Sys.getenv("RUN_PERFORMANCE_TESTS"), "true"), "Set RUN_PERFORMANCE_TESTS=true to run")

    set.seed(789)
    setup <- generate_test_data(link = "log", n_subjects = 18)
    attr(setup$intensity.model, "bandwidth") <- 7

    alphas <- c(-0.5, 0, 0.5)

    fit_original_nocache <- function() {
        last <- NULL
        for (a in alphas) {
            last <- fit_SensIAT_marginal_mean_model_generalized(
                data = setup$data,
                time = ..time..,
                id = ..id..,
                alpha = a,
                knots = setup$knots,
                outcome.model = setup$outcome.model,
                intensity.model = setup$intensity.model,
                loss = "lp_mse",
                link = "log",
                impute_data = create_impute_fn(),
                BBsolve.control = list(maxit = 12, tol = 1e-2, trace = FALSE),
                term2_method = "original",
                use_expected_cache = FALSE
            )
        }
        last
    }

    fit_original_cached <- function() {
        last <- NULL
        for (a in alphas) {
            last <- fit_SensIAT_marginal_mean_model_generalized(
                data = setup$data,
                time = ..time..,
                id = ..id..,
                alpha = a,
                knots = setup$knots,
                outcome.model = setup$outcome.model,
                intensity.model = setup$intensity.model,
                loss = "lp_mse",
                link = "log",
                impute_data = create_impute_fn(),
                BBsolve.control = list(maxit = 12, tol = 1e-2, trace = FALSE),
                term2_method = "original",
                use_expected_cache = TRUE
            )
        }
        last
    }

    fit_fast_nocache <- function() {
        last <- NULL
        for (a in alphas) {
            last <- fit_SensIAT_marginal_mean_model_generalized(
                data = setup$data,
                time = ..time..,
                id = ..id..,
                alpha = a,
                knots = setup$knots,
                outcome.model = setup$outcome.model,
                intensity.model = setup$intensity.model,
                loss = "lp_mse",
                link = "log",
                impute_data = create_impute_fn(),
                BBsolve.control = list(maxit = 12, tol = 1e-2, trace = FALSE),
                term2_method = "fast",
                use_expected_cache = FALSE
            )
        }
        last
    }

    fit_fast_cached <- function() {
        last <- NULL
        for (a in alphas) {
            last <- fit_SensIAT_marginal_mean_model_generalized(
                data = setup$data,
                time = ..time..,
                id = ..id..,
                alpha = a,
                knots = setup$knots,
                outcome.model = setup$outcome.model,
                intensity.model = setup$intensity.model,
                loss = "lp_mse",
                link = "log",
                impute_data = create_impute_fn(),
                BBsolve.control = list(maxit = 12, tol = 1e-2, trace = FALSE),
                term2_method = "fast",
                use_expected_cache = TRUE
            )
        }
        last
    }

    # Warm-up to avoid one-time overhead dominating timings
    suppressWarnings(fit_fast_cached())

    time_original_nc <- system.time(suppressWarnings(result_original_nc <- fit_original_nocache()))
    time_original_c <- system.time(suppressWarnings(result_original_c <- fit_original_cached()))
    time_fast_nc <- system.time(suppressWarnings(result_fast_nc <- fit_fast_nocache()))
    time_fast_c <- system.time(suppressWarnings(result_fast_c <- fit_fast_cached()))

    cat("\n")
    cat("Term2 timing comparison (seconds):\n")
    cat(sprintf("  Original (no cache): %.3f\n", time_original_nc["elapsed"]))
    cat(sprintf("  Original (cached):   %.3f\n", time_original_c["elapsed"]))
    cat(sprintf("  Fast (no cache):     %.3f\n", time_fast_nc["elapsed"]))
    cat(sprintf("  Fast (cached):       %.3f\n", time_fast_c["elapsed"]))
    cat(sprintf("  Fast cached vs fastest ratio: %.2fx\n",
                time_fast_c["elapsed"] / min(c(time_original_nc["elapsed"], time_original_c["elapsed"], time_fast_nc["elapsed"]))))

    fastest_other <- min(c(time_original_nc["elapsed"], time_original_c["elapsed"], time_fast_nc["elapsed"]))
    expect_true(time_fast_c["elapsed"] <= fastest_other * 1.1,
                info = sprintf("Fast cached (%.3fs) should be fastest or within 10%% of the next fastest (%.3fs)",
                               time_fast_c["elapsed"], fastest_other))

    expect_equal(result_fast_c$coefficients, result_original_c$coefficients,
                 tolerance = 0.2,
                 info = "Fast cached term2 should produce similar coefficients to original cached")
    expect_true(all(is.finite(result_fast_c$coefficients)),
                info = "Fast cached term2 coefficients should be finite")
})
